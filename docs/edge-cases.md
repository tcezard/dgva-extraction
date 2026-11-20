# Edge cases not covered by the extraction

The list of DGVa data that dgva-extraction does not yet write to a VCF correctly. Each entry says what happens today
and what it needs.

**Keep it up to date.**
- Add an entry when a run, a check or a review finds a case the code does not handle.
- Remove an entry when the code handles the case and a test covers it, and say so in the commit.
- `tests/test_edge_cases_list.py` fails if a refusal or side-list reason in the code is missing from this file.

Last updated: 2026-11-20, after implementing ticket 10's points 4 (variants with several placements) and 5 (key names).

How the code treats a case today:
- **Refused:** the writer refuses the record, and the whole study fails with no VCF. The statistics list each refusal.
- **Written wrong:** the VCF is valid, but says something DGVa doesn't. **The most urgent group**: nothing flags it.
- **Approximated:** the record is written, but with a default or a less precise allele than DGVa could support.
- **Study fails:** the reader stops the study. None of these has been seen in the data; they are guards.
- **Side list:** the record is not written, and is listed in the study's side list (`*.side_list.tsv`). Sprint 4 turns
  the side list into the DGVa JSON.
- **Not read:** the reader does not read the data, so nothing is written and nothing fails. **These are silent
  losses**, the most important group to close.

## 1. Refused: needs a representation decision (ticket 10)

| Reason code | Case | Where it occurs | What it needs |
|---|---|---|---|
| `no_symbolic_allele` | A type whose allele is `null` in `etc/symbolic_alleles.yaml` | None since 2026-11-20: every type in the table has an allele | A guard, for a type added without a decision |
| `unknown_type` | A type not in `etc/symbolic_alleles.yaml` | `complex chromosomal rearrangement` was missing until 2026-11-20 | Run `dgva-extraction check-types` over the whole database, and add every type in use |
| `missing_repeat_unit` | `<CNV:TR>` (short tandem repeat) with no `RUS` or `RUL` | Any tandem repeat until the repeat fields are read (sprint 4). `nstd128` has motifs on its regions only. | Read the motifs (sprint 4). Then decide what to do with calls that have a count but no motif. |
| `repeat_key_without_tandem_repeat` | A repeat count on a record that isn't `<CNV:TR>` | `nstd128` calls (4.4 M) if their type isn't a tandem repeat: type not checked | Check `nstd128`'s types when the repeat fields are read (sprint 4) |
| `sequence_not_allele` | A submitted sequence with characters other than A, C, G, T, N | Sequences aren't read yet (sprint 4) | Write it in `SUBMITTED_SEQUENCE` instead |
| `alt_equals_ref` | A submitted sequence equal to the reference | Sequences aren't read yet | Decide: drop the sequence, or write it in `SUBMITTED_SEQUENCE` |
| `imprecise_with_sequence` | A sequence-resolved allele with non-zero confidence intervals | Sequences aren't read yet | Decide: symbolic allele with `SUBMITTED_SEQUENCE`, or sequence without the intervals |
| `unknown_chromosome` | A chromosome or contig name that matches neither the FASTA nor the assembly report | Not seen yet with a full FASTA. Includes `Submitted artifact` placements on unplaced contigs (45,964), whose names are untested. | The assembly FASTA and report for each assembly (pipeline, plan step 5.1); check contig names on an artifact study |
| `outside_contig` | A position beyond the end of its sequence | Not seen yet. Possible with the GRCh37.p9/p13 naming mix-up. | Check against the right assembly version |
| `invalid_coordinates` | End before start (not an insertion), or start below 1 | Not seen yet | Look at the case |
| `invalid_confidence_interval` | A confidence interval that does not contain 0 | Not seen yet | Look at the case |
| `invalid_reference_base` | A FASTA character that isn't a base or an IUPAC code | Not seen yet | Look at the FASTA |
| `unmapped_value`, `duplicate_key` | A value with no key in `etc/vcf_keys.yaml`, or a key the writer sets itself | Only when new fields are added | Add the key (ticket 10 naming rules) |
| `invalid_breakend` | A breakend record that spans more than one base, or carries a sequence | A guard: the reader makes breakends only from single-base placements | Look at the case |
| `unknown_event_type` | An `EVENTTYPE` that VCF 4.4 doesn't reserve and dgva-extraction doesn't declare | A guard: the reader writes `TRA` and `COMPLEX` only | Declare it in `DECLARED_EVENT_TYPES` (`vcf_writer.py`), which writes its `##EVENTTYPE` line |

## 1b. Approximated: written, with less than DGVa could give

Decided in ticket 10 on 2026-11-20. Each is valid VCF 4.4; the exact type stays in `VC`.

| Case | Written as | Where it occurs | What would make it exact |
|---|---|---|---|
| `SVCLAIM` on `<DEL>`/`<DUP>` from `deletion`, `duplication`, `tandem duplication` or the mobile element deletions, which don't say whether the call claims abundance or adjacency | `SVCLAIM=D` (abundance), the default: most DGVa technologies (arrays, read depth) measure abundance | `nstd83`: 15 deletions; most studies | Derive it from the study metadata, for example the detection method (sequencing calls may be `J` or `DJ`). A later ticket, if the gain is worth the cost. |
| Translocations with one placement (the partner breakpoint isn't a placement of the variant) | `<TRA>`, defined in the file, with `SVLEN=.` | `estd236` (74 of its 89 rows), `nstd83` (20 call placements and 10 regions), and other studies | Breakend notation where DGVa has the mate (`VARIANT_MATE`, not read: section 4c). Translocations given as two placements are written as breakend pairs since 2026-11-20. |
| Translocations given as two breakpoints with no breakpoint order | A breakend pair, `From` taken as the lower position (natural chromosome order, then position), as dbVar does. The orientation depends on that choice when the two strands differ. | `nstd115` (29,142 variants), `nstd107` (1,383), `nstd99` (1,077), and the `translocation` regions of `nstd106` (2,862) | The real breakpoint order; none is known. Decided 2026-11-20 (ticket 10, point 4). |
| Translocations given as two breakpoints with strand `not reported` or missing | Two `<TRA>` records linked by `EVENT`, `EVENTTYPE=TRA`: breakend notation can't say the orientation is unknown | About 2,150 intrachromosomal translocations: `nstd106` (1,211), `nstd107` (461), `nstd99` (475); a few dozen with one strand missing | The strands; none is known |
| Inversions and CNVs given as two single-base breakpoints on one chromosome | One record from the first breakpoint to the second; the second placement is on the side list (`merged_breakpoint`) | 140 variants in `nstd133` | Nothing, if the breakpoints are the variant's ends, as assumed |
| Complex chromosomal rearrangements given as one placement per breakpoint | One single breakend per breakpoint (`T.` or `.T` from `From`/`To` and the strand), linked by `EVENT`, `EVENTTYPE=COMPLEX` (declared in a `##EVENTTYPE` line). Which breakpoints are joined is not given. | 245 variants in `nstd133`, 2 to 114 breakpoints each | The pairing of the breakpoints; none is known |
| `complex chromosomal rearrangement`, `complex substitution`, `complex structural alteration` | `<COMPLEX_CHROMOSOMAL_REARRANGEMENT>`, `<COMPLEX_SUBSTITUTION>`, `<COMPLEX_STRUCTURAL_ALTERATION>`, defined in the file | `estd236`: 7 regions; others not counted | Nothing, unless a sequence is submitted (then written from it) |
| `sequence alteration` | `<SEQUENCE_ALTERATION>`, defined in the file | `estd1`: 25,193 calls | Nothing, unless a sequence is submitted |
| `indel`, `delins` without a sequence | `<INDEL>`, `<DELINS>`, defined in the file | not counted | The sequence, where DGVa has one (sprint 4) |
| `copy-neutral loss of heterozygosity`, `loss of heterozygosity` | `<COPY_NUMBER_NEUTRAL_LOSS_OF_HETEROZYGOSITY>`, `<LOSS_OF_HETEROZYGOSITY>`, defined in the file | not counted | Nothing |

Their names follow the SO term (ticket 10, decided 2026-11-20).

## 1c. Written wrong: valid VCF that misstates the variant

None known. The variants with several submitted placements on one assembly (71,548 in 10 studies), found on
2026-11-20, were written as several records with the same ID, each a whole variant. Since ticket 10's point 4 they
are written as breakend pairs, single breakends, one spanning record, or records linked by `EVENT`; the remaining
approximations are in section 1b. Not yet run on the database: see section 5.

## 2. Study fails in the reader (guards; none seen in the data)

| Case | Measured in sprint 2 | Check in the code |
|---|---|---|
| Placement method `conflict` (submitted method with a remap score) or an unknown method | 0 | `placement_class()` |
| W2: a type with no SO term, a variant missing from `VARIANT_REGION`/`VARIANT_CALL`, a region with no assertion method, an assembly missing from `GENOME_ASSEMBLY` | 0 | `_check_w2()` |
| `CI_START`/`CI_STOP` with the wrong sign for the study (signed in `estd219`, unsigned elsewhere) | 0 | `_interval_from_*()` |
| Both CI sets (`CIPOS*`/`CIEND*` and `CI_START`/`CI_STOP`) on one placement | 0 | `_confidence_intervals()` |
| Duplicated rows that disagree. The CI and allele tables have no primary key. | not counted | `_single_placement_row()` |
| A parent region in another study | 0 | `_parents()` |

## 3. Not written; on the side list

| Reason code | Case | Where it occurs | What it needs |
|---|---|---|---|
| `no_start` | W3: submitted placement with no start | `estd236`: `esv4011099`, which has another placement that is written | Nothing: the DGVa JSON lists it (sprint 4) |
| `no_chromosome` | W3: a start but no chromosome or contig name | `estd240`: both variants | Find the chromosome, for example in the Exchange XML. Until then the study has an empty VCF. |
| `no_end` | No stop, outer stop or inner stop, on a type that isn't an insertion | not counted | Look at the cases |
| `duplicate_placement` | Identical placements of one accession, written once | `estd208`: 169 variants × 5 placements | Find which assembly they are on (ticket 13, not done) |
| `positioned_from_calls` | A translocation region with no placement of its own, positioned from its first supporting call. It is listed so that the derived position is visible. | `nstd83`: 10 regions | Written as `<TRA>` (section 1b); `BND` will need both calls' positions |
| `merged_breakpoint` | The second breakpoint of an inversion or CNV given as two single-base placements on one chromosome: written as the end of the record of the first | 140 variants in `nstd133` | Nothing: listed so that the merge is visible |
| (study side list) `withdrawn`, `remapped_only`, `no_position_only`, `no_placement` | Variants in no VCF | 101,575 withdrawn; 976 remapped-only (`nstd35`); 751 with no position only; 112 with no placement (`estd180`, `nstd83`) | Nothing: the DGVa JSON lists them (sprint 4). The remapped placements go to the remappings file (D3). |

## 4. Not read: silent losses until they are added

### 4a. Gap fields, planned for sprint 4

The gap table of `reference/dgva-oracle-schema.md`, in the planning repository:
- HGVS names, placement rank and placements per assembly;
- sequences and samplesets (validation: section 4b2);
- repeat counts and motifs, cytobands and mutation order;
- alt status and submitted types;
- `estd214`/`estd219` copy numbers and allele frequencies, from the load tables.

Keys are proposed in `etc/vcf_keys.yaml` for each of these.

### 4b. Samples and genotypes, decided on 2026-11-19 and planned for sprint 4

Reading `VARIANT_CALL_SAMPLE`, `VARIANT_CALL_ZYGOSITY` and `VARIANT_CALL_COPY_NUMBER`, then:
- one VCF per sample;
- sampleset and `<study>_unspecified_samples` pseudo-samples;
- a call with several samples written in each of their files.

See `docs/vcf-representation.md`. `genotype.py` and `etc/ploidy.yaml` exist, but nothing calls them yet.

eva-sub-cli (2026-11-20) requires a sample on every analysis: the regions analysis, and the calls of a study with no
samples or samplesets, must list samples without a VCF column. A study with no sample at all needs a study-level
pseudo-sample in the metadata (open, ticket 14). Pseudo-samples need `collection date` and `geographic location`
(`not provided` is accepted).

### 4b2. Experiments, decided on 2026-11-20 and planned for sprint 4

Not read today, so every call of a study and assembly goes to one calls analysis. Decided (ticket 10, point 5):
- experiments with an identical description (type, method type, analysis type, detection method, platforms) are
  merged; each merged description has an ID and is described once in the metadata;
- one calls analysis per discovery (or genotyping) description and submitted assembly: 363 instead of 206 before
  the multi-assembly studies; 51 studies split, the largest `nstd37` (30);
- `VALIDATION_RESULT=<description ID>:<result>` for the 2,992,979 validation results.

Needs `VARIANT_CALL_EXPERIMENT` (one experiment per call), `VARIANT_CALL_VALIDATION`, and the experiment tables.

### 4c. Fields the GVF already carries, which no ticket covers

The gap table lists only what the GVF lacks. These columns reach the GVF (reference document, "Columns that reach
the GVF"), so a VCF without them loses data that the GVF has. None has a key in `etc/vcf_keys.yaml` yet, and
**no sprint 3 or sprint 4 ticket reads them**:

| Data | Tables | GVF attribute | Legacy VCF key |
|---|---|---|---|
| Submitter variant ID | `VARIANT_STUDY.SUBMITTER_VARIANT_ID` | `Alias`, `submitter_variant_*_id` | `NAME` |
| Description | `VARIANT_DESCRIPTION` | `variant_*_description` | `DESC` |
| Event ID | `VARIANT_EVENT` | `event_id` | |
| Variant attributes: log2 value, support count, loss of heterozygosity, insertion length | `VARIANT_ATTRIBUTE` | same names | `LOG2VAL`, `SUPPNUM`, `LOH` |
| Clinical significance, mode of inheritance, origin | `VARIANT_CALL_CLIN_SIGNIF`, `_MODE_OF_INHERIT`, `_ORIGIN` | `clinical_significance`, `mode_of_inheritance`, `var_origin` | `CLINSIG`, `MOI`, `ORIGIN` |
| Phenotype links and descriptions | `VARIANT_CALL_PHENOTYPE`, `VARIANT_CALL_PHENO_DESCRIPTION` | `phenotype_link`, `phenotype_description` | `PHENOLINK`, `PHENODESC` |
| External links | `VARIANT_CALL_EXT_LINK` → `EXTERNAL_LINK` | `Dbxref` | `DBXREF` |
| Supporting and evidence sequence links | `VAR_CALL_SUPPORT_EXT_LINK` (455,310 calls, 15 studies) | `<type>_link` | `SUPPLINK`, `EVIDLINK`, `3INLINK`, … |
| Breakpoint order | `LINKED_VAR_BREAKPOINT_ORDER`, `VAR_PLAC_BREAKPOINT_ORDER` | `mutation_order`, `breakpoint_order` | `BNDORDER` |
| Strand (on breakpoints it gives the breakend join orientation: section 1c) | `VARIANT_PLACEMENT_STRAND` | column 7 | |
| Allele count and number | `ALLELE_COUNT`, `ALLELE_NUMBER` | `allele_count`, `allele_number` | `AC`, `AN` |
| Allele frequency outside the load tables | `ALLELE_FREQUENCY` | `allele_frequency` | `AF` (key exists, not read) |
| Mates | `VARIANT_MATE` (never read by the exporter, bug B9) | wrong `mate_id` | |

Needed: count which of these hold submitted data, add keys in ticket 10's table, and add them to sprint 4's reader
work alongside the gap fields. Metadata that reaches the GVF (study, subject, sample and publication pragmas) is for
the metadata work (plan step 4).

## 5. Assumptions not yet checked on real data

| Assumption | Risk if wrong | How to check |
|---|---|---|
| The multi-placement plan (ticket 10, point 4) on real data: `multi_placement_variants.sql`, the breakend pairs and the merged breakpoints | A wrong orientation or pairing, or a study failing on a case not foreseen | Extract `nstd106`, `nstd133` and `nstd107` and compare their breakends with dbVar's VCFs (they matched the strand rule exactly) |
| An insertion's sequence goes before `start` (`POS` = start − 1), as convertGVFtoVCF assumes | Insertions off by one | Compare with the Exchange XML or dbVar for a few insertions |
| The database sorts chromosome names like Python does (binary collation) | The translocation regions merged into the stream arrive out of order, and the writer stops (`OrderError`). Nothing is written wrong. | A run on a study with translocation regions on several chromosomes |
| The reader's ordering keeps memory flat on large studies | Out-of-memory error, or slow runs | `estd214`: run time and peak memory (ticket 13.6). Note the unexplained 552 s run of `estd236`. |
| Contig names of `Submitted artifact` placements resolve through the assembly report | `unknown_chromosome` refusals | An artifact study |
| Assembly names and versions: the FTP has GRCh37.p9 where the database has GRCh37.p13 | Wrong FASTA picked | The pipeline takes the assembly accession from the database (plan step 5.1) |
| Withdrawn accessions are dropped whichever study withdrew them, as the exporters do (W4) | A variant withdrawn from one study and valid in another is lost from both | Count accessions withdrawn by a study other than their own |

## 6. Sprint decisions that change what is extracted

| Decision | Status | Effect on the code |
|---|---|---|
| Sprint 3 decision 3: the 32,192 variants only in `VARIANT_STUDY_NEW` | Open; recommended: leave out | Not read, so nothing to change if the recommendation holds |
| Sprint 3 decision 4: the 2 live `estd214` chrY regions (`esv3818028`, `esv3817941`) | Open; recommended: treat as withdrawn | **Not implemented: the reader writes them today.** Needs an exclusion list if the recommendation is taken. |
| Ticket 10 points 1–3: VCF 4.4; `SVCLAIM` defaults to `D`; types with no VCF allele get one the file defines | Decided 2026-11-20 | Implemented; see section 1b |
| Ticket 10 points 4–5: variants with several placements (best placement keeps the accession, `EVENT` links the records, breakends for breakpoints); key names | Decided 2026-11-20 | Implemented, unit-tested; not yet run on the database (section 5). Allele names decided the same day. |
