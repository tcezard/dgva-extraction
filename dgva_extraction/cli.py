"""dgva-extraction command line.

Every command that reads the database runs only with the SELECT-only account: the connection is refused if the
session holds any write privilege (see dgva_extraction.database).
"""
import argparse
import logging
import os
import resource
import sys
import time

from ebi_eva_common_pyutils.logger import logging_config as log_cfg

from dgva_extraction import __version__
from dgva_extraction.config import database_config, load_config
from dgva_extraction.database import DgvaDatabase
from dgva_extraction.keys import SymbolicAlleles
from dgva_extraction.reader import StudyReader, placement_class, study_side_list, write_side_list
from dgva_extraction.reference import FastaReference, Genome
from dgva_extraction.vcf_writer import StudyVcfWriter, VcfLineBuilder

logger = log_cfg.get_logger(__name__)


def peak_memory_mb():
    """Peak resident memory of this process (ru_maxrss is in bytes on macOS, kilobytes on Linux)."""
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return peak / (1024 * 1024) if sys.platform == 'darwin' else peak / 1024


def write_rows(path, columns, rows):
    with open(path, 'w') as handle:
        handle.write('\t'.join(columns) + '\n')
        for row in rows:
            handle.write('\t'.join('' if row.get(column) is None else str(row[column]) for column in columns) + '\n')


def connect(args):
    return DgvaDatabase.connect(database_config())


def check_privileges(args):
    with connect(args):
        print('The session holds no write privilege.')
    return 0


def list_studies(args):
    with connect(args) as database:
        rows = list(database.stream('study_assembly_list', {'study': args.study}))
    if not rows:
        logger.error('No study found')
        return 1
    write_rows(args.output, list(rows[0]), rows)
    submitted = sum(row['submitted_variants'] or 0 for row in rows)
    to_convert = sum(1 for row in rows if (row['submitted_placements'] or 0) > 0)
    print(f'{len({row["study_accession"] for row in rows})} studies, {to_convert} submitted study-assemblies, '
          f'{submitted} submitted variants. Written to {args.output}')
    return 0


def check_types(args):
    alleles = SymbolicAlleles.load()
    with connect(args) as database:
        rows = list(database.stream('so_terms', {}))
    missing = 0
    for row in rows:
        entry = alleles.for_term(row['so_term'] or '')
        allele = 'NOT IN TABLE' if entry is None else (entry.allele or 'none (refused)')
        missing += entry is None and row['variants'] > 0
        print(f'{row["region_or_call"]}\t{row["exchange_type"]}\t{row["so_term"]}\t{row["so_id"]}\t'
              f'{row["variants"]}\t{allele}')
    return 1 if missing else 0


def side_list(args):
    with connect(args) as database:
        entries = study_side_list(database, args.study)
    write_side_list(args.output, entries)
    print(f'{len(entries)} variants of {args.study} are in no VCF. Written to {args.output}')
    return 0


DUPLICATE_COLUMNS = ('study_accession', 'accession', 'variant_type', 'assembly_accession', 'placement_id', 'placement_class',
                     'placement_method', 'placement_rank', 'placements_per_assembly', 'assembly_unit', 'alt_status',
                     'breakpoint_order',
                     'chromosome', 'start_pos', 'stop_pos', 'outer_start_pos', 'inner_start_pos', 'inner_stop_pos',
                     'outer_stop_pos')


def _start_stop(row):
    # The reader's order: start, then outer, then inner coordinates.
    start = next((row.get(key) for key in ('start_pos', 'outer_start_pos', 'inner_start_pos')
                  if row.get(key) is not None), None)
    stop = next((row.get(key) for key in ('stop_pos', 'outer_stop_pos', 'inner_stop_pos')
                 if row.get(key) is not None), None)
    return start, stop


def describe_duplicates(placements):
    """How the submitted placements of one variant on one assembly differ, in a few words."""
    places = [(row['chromosome'], *_start_stop(row)) for row in placements]
    if len(set(places)) == 1:
        return 'identical positions'
    if len({chromosome for chromosome, _, _ in places}) > 1:
        return 'different chromosomes'
    spans = sorted((start, stop) for _, start, stop in places if start is not None and stop is not None)
    overlapping = all(later[0] <= earlier[1] for earlier, later in zip(spans, spans[1:]))
    return 'same chromosome, overlapping' if overlapping else 'same chromosome, apart'


def duplicate_placements(args):
    """Real examples of variants with several submitted placements on one assembly (ticket 10, open point 4)."""
    seeds = {}
    per_study = {}
    variants = set()
    with connect(args) as database:
        # Read every seed to count the variants; keep only a few as examples.
        for row in database.stream('duplicate_placement_seeds', {'study': args.study}):
            key = row['variant_study_id']
            variants.add(key)
            if key in seeds or per_study.get(row['study_accession'], 0) >= args.per_study \
                    or len(seeds) >= args.limit:
                continue
            seeds[key] = row
            per_study[row['study_accession']] = per_study.get(row['study_accession'], 0) + 1
        summary = list(database.stream('duplicate_placement_summary', {'study': args.study}))
        breakpoint_orders = {row['breakpoint_order']: row['description']
                             for row in database.stream('breakpoint_order_cv', {})}
        rows = []
        for key, seed in seeds.items():
            for placement in database.stream('variant_placements', {'variant_study_id': key}):
                rows.append(dict(placement, study_accession=seed['study_accession'], accession=seed['accession'],
                                 placement_class=placement_class(placement['placement_method'],
                                                                 placement['has_remap_score'])))
    write_rows(args.output, DUPLICATE_COLUMNS, rows)
    by_variant = {}
    for row in rows:
        if row['placement_class'] == 'submitted':
            by_variant.setdefault((row['study_accession'], row['accession'], row['assembly_accession']), []).append(row)
    for (study, accession, assembly), placements in by_variant.items():
        ranks = ', '.join(f'{row["placement_rank"]}/{row["breakpoint_order"]}' for row in placements)
        print(f'{study}\t{accession}\t{placements[0]["variant_type"]}\t{assembly}\t'
              f'{len(placements)} submitted placements (rank/breakpoint order: {ranks})\t'
              f'{describe_duplicates(placements)}')
    print('Breakpoint orders: ' + '; '.join(f'{value} = {description}'
                                            for value, description in breakpoint_orders.items()))
    print('Variants with two or more submitted placements on one assembly, whatever their rank:')
    print('type\tplacements\twith breakpoint order\tlowest\thighest\tvariants')
    for row in summary:
        print('\t'.join(str(row[column]) for column in ('variant_type', 'placements', 'with_breakpoint_order',
                                                        'lowest_breakpoint_order', 'highest_breakpoint_order',
                                                        'variants')))
    print(f'{len(variants)} variants of {args.study} have a NotBest submitted placement; {len(seeds)} shown. '
          f'Every placement of those is written to {args.output}')
    return 0


# Regions and calls go to separate VCFs, for separate EVA analyses: the regions are DGVa's clusters of the submitted
# calls, not the raw submitted data (docs/vcf-representation.md).
RECORD_FILES = ('calls', 'regions')


def extract(args):
    started = time.monotonic()
    os.makedirs(args.output_dir, exist_ok=True)
    prefix = os.path.join(args.output_dir, f'{args.study}_{args.assembly}')
    fasta = FastaReference(args.fasta)
    genome = Genome(fasta, args.assembly_report, args.assembly)
    builder = VcfLineBuilder(genome)
    writers = {
        kind: StudyVcfWriter(f'{prefix}.{kind}.vcf', builder, reference_name=args.assembly,
                             meta_lines=[('DGVa_study', args.study), ('DGVa_assembly', args.assembly),
                                         ('DGVa_records', kind)],
                             source=f'dgva-extraction {__version__}')
        for kind in RECORD_FILES
    }
    status = 0
    try:
        with connect(args) as database:
            reader = StudyReader(database, args.study, args.assembly)
            try:
                for record in reader.records():
                    writers['regions' if record.kind.value == 'region' else 'calls'].write(record)
                refused = sum(sum(writer.statistics.refused.values()) for writer in writers.values())
                if refused:
                    # One refused record fails the whole study: neither VCF is written.
                    status = 1
                    logger.error('%s %s failed: %s record(s) refused; no VCF was written. See the statistics.',
                                 args.study, args.assembly, refused)
                    for writer in writers.values():
                        writer.abort()
                else:
                    for kind, writer in writers.items():
                        if writer.statistics.written == 0:
                            writer.abort()
                            continue
                        if reader.counts['null_method_placements']:
                            writer.meta_lines.append(('DGVa_placement_method', 'unknown'))
                        writer.finish()
            finally:
                write_side_list(prefix + '.side_list.tsv', reader.side_list)
                with open(prefix + '.stats.tsv', 'w') as handle:
                    handle.write('file\tcategory\tvalue\tcount\n')
                    for kind, writer in writers.items():
                        for category, value, count in writer.statistics.rows():
                            handle.write(f'{kind}\t{category}\t{value}\t{count}\n')
                        for reason, examples in sorted(writer.statistics.refused_examples.items()):
                            for example in examples:
                                handle.write(f'{kind}\trefused_example\t{reason}\t{example}\n')
                    for name, count in sorted(reader.counts.items()):
                        handle.write(f'all\treader\t{name}\t{count}\n')
                    handle.write(f'all\trun\tseconds\t{time.monotonic() - started:.1f}\n')
                    handle.write(f'all\trun\tpeak_memory_mb\t{peak_memory_mb():.1f}\n')
    finally:
        for writer in writers.values():
            writer.abort()  # no-op for a finished writer; removes temporary files otherwise
        fasta.close()
    written = ', '.join(f'{writers[kind].statistics.written} {kind}' for kind in RECORD_FILES)
    print(f'{args.study} {args.assembly}: {written} written, {len(reader.side_list)} on the side list. '
          f'Statistics in {prefix}.stats.tsv')
    return status


def build_parser():
    parser = argparse.ArgumentParser(prog='dgva-extraction', description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version=__version__)
    parser.add_argument('--config', help='YAML configuration with a DGVA section, kept outside the repository. '
                                         'Default: $DGVA_EXTRACTION_CONFIG, then ~/.dgva_extraction_config.yml')
    parser.add_argument('--debug', action='store_true', default=False,
                        help='Set the logging to debug level')
    commands = parser.add_subparsers(dest='command', required=True)

    commands.add_parser('check-privileges', help='Connect and check that the session cannot write')

    command = commands.add_parser('list-studies', help='Studies and assemblies with their counts per class')
    command.add_argument('--study', default='ALL', help='One study accession, or ALL (default; slow)')
    command.add_argument('--output', required=True)

    commands.add_parser('check-types', help='Variant types in use, and their symbolic allele')

    command = commands.add_parser('side-list', help='Variants of a study that are in no VCF')
    command.add_argument('--study', required=True)
    command.add_argument('--output', required=True)

    command = commands.add_parser('duplicate-placements',
                                  help='Examples of variants with several submitted placements on one assembly')
    command.add_argument('--study', required=True, help='Study accession')
    command.add_argument('--limit', type=int, default=20, help='Number of variants to show (default 20)')
    command.add_argument('--per-study', type=int, default=3, help='At most this many variants per study (default 3)')
    command.add_argument('--output', required=True)

    command = commands.add_parser('extract', help='Write the VCF of one study on one submitted assembly')
    command.add_argument('--study', required=True)
    command.add_argument('--assembly', required=True, help='Assembly accession, as in VARIANT_PLACEMENT')
    command.add_argument('--fasta', required=True, help='Uncompressed FASTA of the assembly (indexed if needed)')
    command.add_argument('--assembly-report', help='NCBI assembly report: contig order and name aliases')
    command.add_argument('--output-dir', required=True)
    return parser


COMMANDS = {'check-privileges': check_privileges, 'list-studies': list_studies, 'check-types': check_types,
            'side-list': side_list, 'duplicate-placements': duplicate_placements, 'extract': extract}


def main(argv=None):
    log_cfg.add_stdout_handler()
    args = build_parser().parse_args(argv)
    if args.debug:
        log_cfg.set_log_level(logging.DEBUG)
    load_config(args.config)
    return COMMANDS[args.command](args)


if __name__ == '__main__':
    sys.exit(main())
