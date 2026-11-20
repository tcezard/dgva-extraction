"""Off-hours checks (sprint 3, tickets 10 and 13): what the full scans tell about the database, as files and a summary.

Run with `dgva-extraction checks --output-dir DIR` outside busy hours: several steps read whole tables. Each step
writes a TSV and adds a section to summary.md; a step that fails is reported and the next one runs. `--study`
restricts the per-study steps to one study, for a short trial during the day.

This module holds the analysis, independent of the database; dgva_extraction.cli runs the queries.
"""
import itertools
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Tuple

# Prefix of a merged experiment description's ID, per experiment type (ticket 10, point 5).
DESCRIPTION_PREFIXES = {'Discovery': 'D', 'Validation': 'V', 'Genotyping': 'G'}
DESCRIPTION_PARTS = ('method_type', 'analysis_type', 'detection_method', 'platform_name')
# Experiment types whose calls make calls analyses (validation experiments validate calls found by others).
CALLING_TYPES = ('Discovery', 'Genotyping')
# Submitted variants counted in sprint 2 (ticket 6), to check list-studies against.
SPRINT2_SUBMITTED_VARIANTS = 44156767


@dataclass
class Experiment:
    experiment_id: int
    study: str
    submitter_id: int
    experiment_type: str
    has_calls: bool
    parts: Dict[str, set] = field(default_factory=lambda: {part: set() for part in DESCRIPTION_PARTS})

    @property
    def signature(self):
        return (self.experiment_type,) + tuple(tuple(sorted(self.parts[part])) for part in DESCRIPTION_PARTS)


@dataclass
class Descriptions:
    """Experiments merged per study when their whole description is identical (ticket 10, point 5).

    Each description gets an ID per study and experiment type, in the order of the lowest submitter experiment ID it
    covers: D1, D2... for discovery, V1... for validation, G1... for genotyping.
    """
    experiments: Dict[int, Experiment]
    description_of: Dict[int, str]

    @classmethod
    def from_rows(cls, rows: Iterable[dict]):
        experiments: Dict[int, Experiment] = {}
        for row in rows:
            experiment = experiments.get(row['experiment_id'])
            if experiment is None:
                experiment = experiments[row['experiment_id']] = Experiment(
                    row['experiment_id'], row['study_accession'], row['submitter_experiment_id'],
                    row['experiment_type'], bool(row['has_calls']))
            for part in DESCRIPTION_PARTS:
                if row[part] is not None:
                    experiment.parts[part].add(row[part])
        description_of = {}
        by_study_type = defaultdict(list)
        for experiment in experiments.values():
            by_study_type[(experiment.study, experiment.experiment_type)].append(experiment)
        for (study, experiment_type), group in by_study_type.items():
            signatures = {}
            for experiment in sorted(group, key=lambda item: (item.submitter_id, item.experiment_id)):
                if experiment.signature not in signatures:
                    prefix = DESCRIPTION_PREFIXES.get(experiment_type, 'X')
                    signatures[experiment.signature] = f'{prefix}{len(signatures) + 1}'
                description_of[experiment.experiment_id] = signatures[experiment.signature]
        return cls(experiments, description_of)

    def rows(self):
        for experiment in sorted(self.experiments.values(), key=lambda item: (item.study, item.experiment_type,
                                                                             item.submitter_id)):
            yield {
                'study_accession': experiment.study, 'experiment_type': experiment.experiment_type,
                'description_id': self.description_of[experiment.experiment_id],
                'submitter_experiment_id': experiment.submitter_id, 'experiment_id': experiment.experiment_id,
                'has_calls': int(experiment.has_calls),
                **{part: '; '.join(sorted(experiment.parts[part])) for part in DESCRIPTION_PARTS},
            }

    def studies(self):
        return sorted({experiment.study for experiment in self.experiments.values()})

    def summary(self):
        lines = []
        for experiment_type in sorted({item.experiment_type for item in self.experiments.values()}):
            group = [item for item in self.experiments.values() if item.experiment_type == experiment_type]
            per_study = defaultdict(set)
            for item in group:
                per_study[item.study].add(self.description_of[item.experiment_id])
            several = sum(1 for descriptions in per_study.values() if len(descriptions) > 1)
            lines.append(f'- {experiment_type}: {len(group)} experiments ({sum(item.has_calls for item in group)} '
                         f'with calls) in {len(per_study)} studies; {sum(len(v) for v in per_study.values())} '
                         f'merged descriptions; {several} studies with more than one.')
        return lines


def calls_analyses(rows: Iterable[dict], descriptions: Descriptions) -> Dict[Tuple[str, str, str], int]:
    """Calls per (study, submitted assembly, discovery or genotyping description): one EVA calls analysis each.

    A call is counted once per assembly; the experiments of one description add up.
    """
    analyses: Dict[Tuple[str, str, str], int] = Counter()
    for row in rows:
        experiment = descriptions.experiments.get(row['experiment_id'])
        if experiment is None:
            description = f'unknown experiment {row["experiment_id"]}'
        elif experiment.experiment_type not in CALLING_TYPES:
            description = f'{experiment.experiment_type} experiment {descriptions.description_of[row["experiment_id"]]}'
        else:
            description = descriptions.description_of[row['experiment_id']]
        analyses[(row['study_accession'], row['assembly_accession'], description)] += int(row['calls'])
    return analyses


def calls_analyses_summary(analyses: Dict[Tuple[str, str, str], int], studies_read: int):
    per_study_assembly = Counter((study, assembly) for study, assembly, _ in analyses)
    per_study = Counter(study for study, _, _ in analyses)
    largest = per_study.most_common(5)
    odd = sorted({description for _, _, description in analyses if not description[:1] in ('D', 'G')
                  or ' ' in description})
    lines = [
        f'- {studies_read} studies read; {len(per_study)} with calls on a submitted assembly.',
        f'- {len(per_study_assembly)} study-assemblies with calls: today, one calls analysis each.',
        f'- **{len(analyses)} calls analyses** with one per discovery description and submitted assembly; '
        f'{sum(1 for count in per_study.values() if count > 1)} studies have more than one.',
        f'- Largest: {", ".join(f"{study} ({count})" for study, count in largest)}.',
        f'- Studies on several submitted assemblies: '
        f'{sum(1 for study in per_study if len({a for s, a, _ in analyses if s == study}) > 1)}.',
    ]
    if odd:
        lines.append(f'- **Calls linked to an experiment that is not discovery or genotyping, or unknown:** '
                     f'{", ".join(odd[:10])}.')
    return lines


@dataclass
class ValidationSummary:
    calls: int = 0
    several_descriptions: int = 0
    mixed_results: int = 0
    mixed_within_one_description: int = 0
    unknown_experiments: int = 0
    results: Counter = field(default_factory=Counter)
    examples: List[dict] = field(default_factory=list)

    def lines(self):
        return [
            f'- {self.calls} calls validated more than once ({sum(self.results.values())} results: '
            f'{", ".join(f"{key} {value}" for key, value in sorted(self.results.items()))}).',
            f'- **{self.several_descriptions} calls validated by more than one validation description** (method).',
            f'- **{self.mixed_results} calls with mixed results** (for example Pass and Fail), '
            f'{self.mixed_within_one_description} of them within one description.',
            f'- Results whose experiment is unknown: {self.unknown_experiments}.',
        ]


def validation_summary(rows: Iterable[dict], descriptions: Descriptions, max_examples=1000) -> ValidationSummary:
    """Rows of check_multiple_validations.sql, ordered by call accession."""
    summary = ValidationSummary()
    for accession, group in itertools.groupby(rows, key=lambda row: row['accession']):
        group = list(group)
        summary.calls += 1
        per_description = defaultdict(set)
        for row in group:
            summary.results[row['validation_result']] += 1
            description = descriptions.description_of.get(row['experiment_id'])
            if description is None:
                summary.unknown_experiments += 1
                description = f'unknown experiment {row["experiment_id"]}'
            per_description[description].add(row['validation_result'])
        results = set().union(*per_description.values())
        several = len(per_description) > 1
        mixed = len(results) > 1
        summary.several_descriptions += several
        summary.mixed_results += mixed
        summary.mixed_within_one_description += any(len(values) > 1 for values in per_description.values())
        if (several or mixed) and len(summary.examples) < max_examples:
            study = descriptions.experiments.get(group[0]['experiment_id'])
            summary.examples.append({
                'accession': accession, 'study_accession': study.study if study else '',
                'validations': ', '.join(f'{description}:{"/".join(sorted(values))}'
                                         for description, values in sorted(per_description.items())),
                'several_descriptions': int(several), 'mixed_results': int(mixed)})
    return summary
