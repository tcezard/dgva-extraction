# VCF representation of the DGVa data (sprint 3, ticket 10)

**Status: proposed, for review by a second developer.** The writer implements these rules, so they can be tested and
checked with vcf-validator, but nothing here is agreed until the review. Changing a rule after the studies are
submitted means resubmitting them.

The keys themselves are in [`dgva_extraction/etc/vcf_keys.yaml`](../dgva_extraction/etc/vcf_keys.yaml), and the
symbolic allele for each variant type in [`dgva_extraction/etc/symbolic_alleles.yaml`](../dgva_extraction/etc/symbolic_alleles.yaml).
Those two files are the single source: the writer builds the header from them, checks every value against them, and
refuses any value they don't list. Each key has its `Number`, `Type`, `Description`, database source, whether VCF
4.4 reserves it, its status and the sprint whose reader fills it.

The representation is checked with vcf-validator 0.10.2 on VCFs written by the writer (`tests/test_vcf_writer.py`,
`tests/test_cli.py`), and passed **eva-sub-cli validation** on 2026-11-20 (see "eva-sub-cli validation" below).

## eva-sub-cli validation (ticket 10, step 4; 2026-11-20)

eva-sub-cli 0.7.2.dev4 (local clone at `2cd9731`), run natively (Nextflow, vcf-validator, vcf-assembly-checker,
biovalidator), on a test submission written by dgva-extraction's writer on GRCh37 chromosomes 21 and 22:
- a regions VCF with no samples (`AF=.`): `<DEL>` with `SVCLAIM` and confidence intervals, `<DUP:TANDEM>`, `<CNV>`,
  `<INV>`, `<INS:ME>`, the alleles the file defines (`<COMPLEX_CHROMOSOMAL_REARRANGEMENT>`, `<SEQUENCE_ALTERATION>`,
  `<LOSS_OF_HETEROZYGOSITY>`, `<TRA>`), alternative locations linked by `EVENT` with `PLACEMENT_RANK`, a `<TRA>` pair
  linked by `EVENT`, region-to-region `PARENT_REGION_ACCESSION`;
- a sampleset pseudo-sample calls VCF (`GT` 0/1, 1/1, 1/.): `<DEL>`, `<DUP>`, a sequence-resolved insertion, a
  breakend pair with `MATEID`, single breakends with `EVENTTYPE=COMPLEX` and its `##EVENTTYPE` line,
  `VALIDATION_RESULT=V1:Pass,V2:Fail`.

Result: **SUCCESS**. vcf-validator (inside eva-sub-cli): no error in either file; assembly check: 13/13 and 7/7 REF
alleles match; evidence types detected as designed (`allele_frequency` for the regions, `genotype` for the calls);
sample names concordant; FASTA sequences INSDC, compatible with `GCA_000001405.1`.

The first run failed on the **metadata** only, which says what the sprint 4 metadata must do:
- **Every analysis needs at least one sample**, even with `allele_frequency` evidence and no sample column ("No
  sample found for the analysis"). The regions analysis, and the calls analysis of a study with no samples or
  samplesets, must list samples without `sampleInVCF`: the study's samples or pseudo-samples. A study with no sample
  at all needs a study-level pseudo-sample in the metadata (open).
- **A pseudo-sample's `bioSampleObject` needs the ERC000011 fields** `collection date` and `geographic location
  (country and/or sea)`; the missing-value term `not provided` is accepted.
- The Data Access & Sharing statement must be acknowledged (`statements`), as for any submission.
- `fileSize` given in the metadata JSON came out empty after eva-sub-cli's formatting step and failed the schema
  ("must be number"); leaving out `fileSize` and `md5` lets eva-sub-cli compute them. A tool quirk to report.

Not covered: the sprint 4 keys (no values yet), `experimentType` (a fixed list: the DGVa method types need a mapping)
and `platform` (one string per analysis).

## Keys

Names decided on 2026-11-20 (point 5 below): a VCF 4.4 reserved key first, then a key EVA already writes, then a new
long, descriptive name; no key for a value encoded elsewhere. `etc/vcf_keys.yaml` has every key with its source.

| Key | Field | Number | Type | Reserved (4.4) | Source | Sprint |
|---|---|---|---|---|---|---|
| `END` | INFO | 1 | Integer | yes | stop (or outer/inner stop); `POS` for insertions | 3 |
| `SVLEN` | INFO | A | Integer | yes | stop − start + 1; `.` for insertions and the alleles the file defines | 3 |
| `IMPRECISE` | INFO | 0 | Flag | yes | any CI bound not 0 | 3 |
| `CIPOS`, `CIEND` | INFO | 2 | Integer | yes | `CIPOS*`/`CIEND*`, else `CI_START`/`CI_STOP`, else inner/outer | 3 |
| `SVCLAIM` | INFO | A | String | yes | `D` by default, which copy number gain/loss also imply (point 2 below) | 3 |
| `VC` | INFO | 1 | String | no; EVA release key | SO accession of the type, as in the EVA release (`VC=SO:0000159`) | 3 |
| `SUBMITTED_VARIANT_TYPE` | INFO | 1 | String | no | submitter type, where it differs beyond case and separators | 4 |
| `PARENT_REGION_ACCESSION` | INFO | . | String | no | `LINKED_VARIANT`: call → region and region → region | 3 |
| `REGION_ASSERTION_METHOD` | INFO | 1 | String | no | region assertion method; also in the metadata when the same for the whole study | 3 |
| `PLACEMENT_METHOD` | INFO | 1 | String | no | placement method, when not `Submitted genomic`; the study default in the metadata | 3 |
| `EVENT` | INFO | A | String | yes | the accession, on every record of a variant written as several records (point 4) | 3 |
| `EVENTTYPE` | INFO | A | String | yes | `TRA`, or `COMPLEX` (declared in a `##EVENTTYPE` line) (point 4) | 3 |
| `MATEID` | INFO | A | String | yes | the other breakend of a translocation given as two placements (point 4) | 3 |
| `PLACEMENT_RANK` | INFO | 1 | String | no | `VAR_PLAC_PLACEMENT_RANK`; on variants written as several records from sprint 3 | 3 and 4 |
| `CN` | INFO | A | Float | yes | copy number (incl. the `estd214`/`estd219` load tables) | 4 |
| `AF` | INFO | A | Float | yes | allele frequency (incl. the `estd219` load table); `.` on every line of a VCF with no samples | 3 |
| `RUC`, `RUS` | INFO | A, . | Float, String | yes | repeat counts and motifs (`nstd128`) | 4 |
| `HGVS`, `CYTOBAND` | INFO | ., 1 | String | no | `VARIANT_PLACEMENT_HGVS`, `VARIANT_PLACEMENT_CYTOBAND` | 4 |
| `MUTATION_ORDER` | INFO | . | Integer | no | `VAR_PLAC_MUTATION_ORDER` (369 rows) | 4 |
| `PLACEMENT_ALT_STATUS` | INFO | 1 | String | no | `VARIANT_PLACEMENT_ALT_STATUS` (`nstd19` only); not a repeat of the rank (point 5) | 4 |
| `VALIDATION_RESULT` | INFO | . | String | no | `<experiment description ID>:<result>` (`Pass`, `Fail`, `Inconclusive`), one per validation; replaces the reserved `VALIDATED` (point 5) | 4 |
| `SUBMITTED_SEQUENCE` | INFO | 1 | String | no | a sequence that is not an allele (gap row 89) | 4 |
| `SUBMITTER_VARIANT_ID`, `DESCRIPTION`, `PHENOTYPE`, `EXTERNAL_LINKS` | INFO | 1, 1, ., . | String | no | `edge-cases.md` section 4c | 4 |
| `CLNSIG` | INFO | . | String | no | clinical significance; the ClinVar/dbVar name, kept so that tools recognise it | 4 |
| `AC`, `AN`, `SOMATIC` | INFO | A, 1, 0 | Integer, Integer, Flag | yes | `edge-cases.md` section 4c | 4 |
| `GT` | FORMAT | 1 | String | yes | zygosity of the call's sample and the species' ploidy (rule 4) | 4 |
| `CN` | FORMAT | 1 | Float | yes | copy number of the call's sample | 4 |

Dropped: `PLACEMENTS` (encoded by `EVENT` and the record IDs), `SAMPLESET` (encoded by the sample column and the
metadata), the validation experiment (metadata). Old names: `TYPESOID` → `VC`, `ASSERTM` → `REGION_ASSERTION_METHOD`,
`PLACEMETHOD` → `PLACEMENT_METHOD`, `PARENT` → `PARENT_REGION_ACCESSION`, `PLACRANK` → `PLACEMENT_RANK`, `ALTSTATUS` →
`PLACEMENT_ALT_STATUS`, `MUTORDER` → `MUTATION_ORDER`, `SVSEQ` → `SUBMITTED_SEQUENCE`, `SUBMITTED_TYPE` →
`SUBMITTED_VARIANT_TYPE`, `VALIDMETH`/`VALIDRES` → `VALIDATION_RESULT` (the reserved `VALIDATED` flag is not used).

Header lines and INFO keys are written in the order of `vcf_keys.yaml`, and only the keys that a VCF uses are in its
header.

## Record rules

1. **One VCF record per DGVa accession and placement**, with the accession (`esv`/`nsv`/`essv`/`nssv`) in the ID
   column. Records at the same position are never merged. A variant with several distinct submitted placements on the
   assembly is the exception (point 4 below): its best placement keeps the accession as ID, the others get
   `<accession>_2`, `_3`..., and all carry `EVENT=<accession>`.
2. **Regions and calls are both written, in separate VCFs for separate EVA analyses** (agreed 2026-11-19). The regions
   are DGVa's calls made on top of the submitted calls, not the raw submitted data, so they get their own analysis.
   A call points at its region(s) with `PARENT_REGION_ACCESSION`, and so does a region nested in another region (`estd1`, `estd20`,
   `estd59`). A call with no region (W1: `estd236`, `estd237`) has no `PARENT`. A withdrawn parent is left out of
   `PARENT`. See "File layout" below.
3. **A value with no target fails the study.** The writer refuses a record whose type isn't in the allele table, or that
   carries a value with no key in the table, or a value that doesn't match its key's `Number` and `Type`. A study with
   a refused record produces no VCF; its statistics list every refusal with its reason.
4. **Genotypes come from the zygosity and the species' ploidy** (agreed 2026-11-19; `dgva_extraction/genotype.py`,
   ploidy map in `dgva_extraction/etc/ploidy.yaml`, keyed by taxonomy ID). DGVa has no genotype table: a call belongs
   to one sample, and `VARIANT_CALL_ZYGOSITY` may give `Heterozygous`, `Homozygous` or `Hemizygous`.

   | Zygosity | Diploid | Haploid sequence (mitochondrion, chloroplast) |
   |---|---|---|
   | `Heterozygous` | `0/1` | refused |
   | `Homozygous` | `1/1` | `1` |
   | `Hemizygous` | `1` | `1` |
   | none recorded | `1/.`: the sample carries the call on one copy, and the other is unknown | `1` |

   A species missing from the ploidy map fails the study when a genotype has to be written. Sex chromosomes aren't in
   the map: their ploidy depends on the sample's sex, which `Hemizygous` gives when the submitter recorded it.
   Sprint 3's files have no sample columns; reading the samples and zygosities is sprint 4.
5. **Placement method.** `PLACEMENT_METHOD` is written on each record whose method is not `Submitted genomic`
   (`Submitted cytogenetic`, `Submitted artifact`). For the 57 studies with no method at all, the header has
   `##DGVa_placement_method=unknown`.
6. **Coordinates.**
   - Start is the first affected base, from start, else outer start, else inner start; end likewise from stop, outer
     stop, inner stop.
   - A symbolic allele is padded with the preceding base: `POS` = start − 1, and `REF` is that base. At start 1, `POS`
     stays 1 and `REF` is the base at 1.
   - An insertion is placed before `start`: `POS` = start − 1, `END` = `POS`, and `SVLEN` = `.` unless the sequence is
     known.
   - A sequence-resolved allele (sprint 4) is padded with the preceding base, or with the following base at position 1.
7. **Confidence intervals.**
   - The sources, in order: `CIPOS*`/`CIEND*` (unsigned distances); then `CI_START`/`CI_STOP` (signed in `estd219`,
     unsigned elsewhere; a value with the wrong sign fails the study); then the inner and outer coordinates.
   - A missing bound is written `.`.
   - `IMPRECISE` is set when any bound is not 0.
   - `0,0` intervals are written, without `IMPRECISE`.
8. **Reference bases.** `REF` is upper case. An IUPAC code is replaced by the first base it stands for, in alphabetical
   order (as convertGVFtoVCF does).
9. **CHROM.** `CHROM` is the sequence name in the assembly FASTA. A database name is matched as it is, then through the
   assembly report (sequence name, GenBank and RefSeq accessions, UCSC name, `chr` prefix). Contigs are written in
   assembly-report order, and only the ones used are in the header.
10. **Determinism.** Lines are sorted by contig in assembly-report order, then position, accession and line text.
   Running the writer twice on the same records gives byte-identical files.

## File layout and evidence type (agreed 2026-11-19)

eva-sub-cli (`add-sequence-collections`) needs one evidence type per analysis (`evidence_type_checker.py:30-40`):
- `genotype`: sample columns, with `GT` on every sample of every line;
- `allele_frequency`: no samples, with `AF` (or `AC` and `AN`) on every line.

Only the key's presence is checked (`file_utils.py:91-105`). Its own detection, run on test VCFs, classes a file with
`AF=.` on every line as `allele_frequency`, and files with `GT=.` or `GT=1/.` as `genotype`. vcf-validator 0.10.2
accepts all three.

Per study and submitted assembly. The calls are further split by discovery description (point 5 of the decisions:
one calls analysis per merged experiment description, 2026-11-20); the files below are per calls analysis, with the
description ID in their name (for example `<study>_<assembly>_<description>.calls.vcf`) when the study has several.

| VCF | Records | Evidence type | EVA analysis |
|---|---|---|---|
| `<study>_<assembly>.regions.vcf` | the regions | `allele_frequency` (`AF=.`) | its own analysis |
| `<study>_<assembly>.calls.vcf` | the calls, when the study has neither samples nor samplesets on its calls | `allele_frequency` (`AF`, or `.`) | the calls analysis |
| one VCF per sample (sprint 4) | the calls of that sample, when the study has samples; a call with several samples is in each of their files | `genotype` | the calls analysis |
| one VCF per sampleset pseudo-sample (sprint 4) | the calls linked only to that sampleset, whether or not the study also has individual samples | `genotype` | the calls analysis |
| `<study>_unspecified_samples` pseudo-sample VCF (sprint 4) | the calls with neither a sample nor a sampleset, in a study whose other calls have one | `genotype` | the calls analysis |

- **`AF` is never invented.** A VCF with no sample columns gets `AF=.` on every line where DGVa has no frequency
  (the writer adds it). A mock `AF=1` was rejected: it would state that every copy carries the variant.
- **One VCF per sample** avoids millions of empty sample cells (`estd214`: about 6.9 M calls × 2,500 samples), since
  each call belongs to one sample. `extract` writes the regions and calls files today; the split by sample comes with
  the sample reader (sprint 4).
- A study with no region writes no regions file, and a refused record in either file fails the whole study: neither
  file is written.

Calls with no individual sample, and calls with several samples (agreed 2026-11-19):
- **A call linked only to a sampleset gets a pseudo-sample named after the sampleset.** These are the 7.5 M calls in
  `VARIANT_CALL_SAMPLESET` with no `VARIANT_CALL_SAMPLE` row. The pseudo-sample stands for the sampleset (a pool), so
  its calls go in one per-sample VCF in the calls analysis, with `genotype` evidence. This applies to studies that
  also have individual samples and to sampleset-only studies alike. Only a study whose calls have neither samples nor
  samplesets gets a single calls VCF with `AF`.
  - **Name:** the sampleset's submitter ID (`DGVA_SAMPLESET.SUBMITTER_SAMPLESET_ID`), as `sampleInVCF` and the VCF
    column name. A pseudo-sample name equal to a real sample name in the same study fails the study, rather than
    being renamed silently.
  - **`GT`:** no zygosity is recorded for a pool, so rule 4 gives `1/.` (diploid): at least one copy in the pool
    carries the call.
  - **Only calls with no individual sample** are written for a pseudo-sample. A call linked to samples and to a
    sampleset is written in its samples' files only, and the metadata keeps the link to the sampleset (no INFO key:
    point 5).
  - **The sample metadata must say what it is.** In the eva-sub-cli JSON, each pseudo-sample is a `bioSampleObject`
    (there is no BioSample accession for a sampleset). Its characteristics must say that it is not an individual:
    - a `description` of the form "DGVa sampleset <submitter ID> of study <study>: a pool of <SAMPLESET_SIZE>
      samples, not an individual. Its calls were submitted for the sampleset without an individual sample.";
    - the sampleset's name, size, population, sex and phenotypes, from `DGVA_SAMPLESET` (gap rows 15–21, 57–59).

    The DGVa JSON lists the same pseudo-samples, with their sampleset ID, so that the VCF columns can be traced
    back. This is part of the metadata work in sprint 4 (plan step 4.1).
- **A call with neither a sample nor a sampleset, in a study whose other calls have one, gets a study-level
  pseudo-sample** (agreed 2026-11-19). DGVa's rule that a call has one or the other is not enforced. These calls go in
  the same calls analysis, in one more per-sample VCF with `genotype` evidence, and `GT` follows rule 4 (`1/.` when
  diploid).
  - **Name:** `<study>_unspecified_samples` (for example `estd1_unspecified_samples`). A real sample or sampleset with
    that name fails the study.
  - **Sample metadata:** its `description` says "DGVa study <study>: calls submitted without a sample or a sampleset.
    The samples they come from are not recorded." The DGVa JSON lists it too.
  - A study where no call has a sample or a sampleset keeps the single calls VCF with `AF=.`.
- **A call linked to several samples is written in each sample's VCF**, under the same accession. It can be merged
  later if needed. The writer's one-record-per-accession rule applies within a file.

Still open with this layout:
- **Number of files:** about 2,500 per assembly for `estd214` and 2,750 for `nstd37`. eva-sub-cli validates each file,
  and eva-submission brokers each to ENA. Count the samples per study, and check with the core team (ticket 14).
- **eva-sub-cli's manual fallback parser** (`file_utils.py:128`, used only when pysam fails) tests `find('AF=')` as
  a boolean: true when `AF=` is missing, false when `AF` is the first INFO key. Report it to the eva-sub-cli
  maintainers.

## Choices from the gap table (ticket 10, step 2)

| Choice | Proposal |
|---|---|
| Sequences that are not an allele (row 89) | INFO `SUBMITTED_SEQUENCE` |
| `PLACEMENTS` always, or only when not `Single` (row 116) | **Dropped** (point 5): `EVENT` and the record IDs show the variants with several placements. |
| `SUBMITTED_VARIANT_TYPE` (rows 86, 108) | Only where it differs from the Exchange type beyond case and `_`/`-` |
| Placement method | Per line for non-default methods; header line for the null-method studies (rule 5) |

## Decisions on the representation

[`edge-cases.md`](edge-cases.md) lists every case the extraction does not yet handle, including those that
depend on the decisions below.

These came out of checking the writer's output with vcf-validator, or have no answer in the gap table. Points 1–3 were
decided on 2026-11-20 and are implemented, as are points 4 and 5 and the allele names of point 3, decided the same
day. Point 6 is open.

1. **VCF version: 4.4 (decided).** 4.4 reserves the SV keys the gap table uses (`RUC`, `RUS`, `SVCLAIM`, `CNV:TR`).
   Its rules, which vcf-validator enforces and the writer follows:
   - `SVLEN` on every symbolic allele: `.` for `<INS>` and for the alleles the file defines (point 3);
   - `SVCLAIM` on every `<DEL>` and `<DUP>` (point 2);
   - `RUS` or `RUL` on every `<CNV:TR>`;
   - the repeat keys only on `<CNV:TR>`.
2. **`SVCLAIM` for `<DEL>`/`<DUP>`: abundance (`D`) by default (decided).** The technologies behind most DGVa calls
   (arrays, read depth) measure abundance and rarely resolve an adjacency (`J`). A copy number gain or loss is an
   abundance claim by definition; `deletion`, `duplication`, `tandem duplication` and the mobile element deletions get
   the default. A submitted `SVCLAIM` would be kept.

   A later ticket can refine it from the study metadata, if the gain is worth the cost. For example, from the analysis
   type of the call's experiment (`VARIANT_CALL_EXPERIMENT` → `EXPERIMENT_ANALYSIS` → `DGVA_ANALYSIS.ANALYSIS_TYPE`):
   - `Paired-end mapping`, `Split read mapping`, `Split read and paired-end mapping`, `Sequence alignment` and the
     assembly-based analyses would give `J`;
   - `Read depth and paired-end mapping` would give `DJ`;
   - the others keep `D`. Regions have no experiment of their own; their claim could come from their calls.
3. **Types with no VCF 4.4 symbolic allele: a symbolic allele the file defines (decided).** A record with a sequence
   is written from it, and a type that fits a VCF 4.4 allele uses that allele. Otherwise the type gets an allele defined
   in the `##ALT` header of each VCF that uses it, as VCF 4.4 allows. These can't have colon-separated subtypes (only
   `DEL`, `INS`, `DUP`, `INV` and `CNV` can), and are written with `SVLEN=.`. The exact type stays in `VC`.

   **Names (decided 2026-11-20): the SO term name in upper case, one allele per term**, so that the allele alone says
   the type; except the translocations, which share `<TRA>`, the name VCF 4.4 reserves as an `EVENTTYPE`. The rule
   also avoids short names other tools give a narrower meaning (GATK-SV's `<CPX>`).

   | Types | Allele |
   |---|---|
   | translocation, interchromosomal translocation, intrachromosomal translocation | `<TRA>` |
   | complex structural alteration | `<COMPLEX_STRUCTURAL_ALTERATION>` |
   | complex chromosomal rearrangement | `<COMPLEX_CHROMOSOMAL_REARRANGEMENT>` |
   | complex substitution | `<COMPLEX_SUBSTITUTION>` |
   | sequence alteration (`estd1`: 25,193 calls) | `<SEQUENCE_ALTERATION>` |
   | `indel`, `delins` without a sequence | `<INDEL>`, `<DELINS>` |
   | copy-number neutral loss of heterozygosity | `<COPY_NUMBER_NEUTRAL_LOSS_OF_HETEROZYGOSITY>` |
   | loss of heterozygosity | `<LOSS_OF_HETEROZYGOSITY>` |

   Translocations given as two placements are written in breakend notation (point 4); `<TRA>` stays for the others,
   such as the `nstd83` regions that the reader positions from their supporting calls, until `VARIANT_MATE` is read.
4. **Variants with several placements (decided 2026-11-20; option C).** A variant with two submitted placements on one
   assembly gives two records with the same ID. vcf-validator 0.10.2 accepts this, but VCF 4.4 says no identifier
   should be present in more than one data record.

   **What the cases are** (`dgva-extraction duplicate-placements`, 2026-11-20). Variants with a `NotBest` submitted
   placement are in 10 studies only: 71,548 variants. All have exactly 2 submitted placements on the assembly, except
   245 `nstd133` complex chromosomal rearrangements, which have 2 to 114 (one per breakpoint). "Best" and "not best"
   mean different things depending on the case:

   | Case | Studies (variants) | Example | What the two placements are |
   |---|---|---|---|
   | **a. The two breakpoints of one rearrangement** | `nstd115` (58,284), `nstd106` (5,752), `nstd107` (2,766), `nstd99` (2,154), `nstd133` (1,167), `estd210` (2) | `nstd99` `nssv3764364` interchromosomal translocation: chr17:38,378,792 and chr5:79,718,097. `nstd133` `nssv13659711` copy number loss: chr14:31,269,856 and chr14:39,513,878. | Two single-base placements, one per breakpoint: translocations (most), inversions and large CNVs (`nstd133`). Which one is "best" is arbitrary. |
   | **b. Alternative locations of the same sequence** | `nstd116` (1,286), `nstd102` (4) | `nstd116` `nssv8640080` CNV: chr8:7,286,416–7,366,833 or chr8:7,669,242–7,740,180, the 8p23.1 segmental duplications. `nstd102` `nsv1398014` CNV: X:60,814–537,179 or Y:10,814–487,179, in PAR1. | The same variant, which can't be told apart between two copies of a sequence. |
   | **c. One variant in pieces** | `nstd10` (131), `nstd40` (2) | `nstd10` `nsv429273` CNV on CanFam2 (`GCF_000002285.1`): AAEX02001078.1:1,121–11,033 and AAEX02001079.1:1–321, two consecutive WGS contigs. `nstd40` `nsv482237` CNV, cytogenetic: 6:164.5–171.1 Mb and 11:130.8–135.0 Mb. | Parts of one variant: split across contigs, or two segments of one rearrangement. |

   **Case a was written wrong, not just with a repeated ID**, until this decision: each breakpoint became its own
   record of the whole variant (`nstd133`'s 8 Mb deletion `nssv13659711` gave two 1 bp `<DEL>` records, and each
   translocation two unlinked `<TRA>` records).

   **Options considered:** A, one record per placement with the accession as ID in each; B, only the `BestAvailable`
   placement, the others in the DGVa JSON; C, the best placement keeps the accession as its ID and each other
   placement gets an ID of its own, linked by an INFO key. **C was chosen**: nothing is hidden, and for breakpoints it
   is what breakend notation needs anyway (one ID per breakend, linked by `MATEID`). VCF 4.4 also names `EVENT` for
   "variant position ambiguity in segmentally duplicated regions", which is case b.

   **`VAR_PLAC_BREAKPOINT_ORDER` doesn't identify case a** (checked 2026-11-20 with `duplicate-placements`, which now
   prints a summary per study). Its values are `From` and `To` ("Links breakpoints using From and to"), not numbers.
   - It is absent from every placement of cases b and c (`nstd116`, `nstd102`, `nstd10`, `nstd40`).
   - In case a it is on the placements of only 3 of the 6 studies: the inter- and intrachromosomal translocation calls
     of `nstd106` (2,890) and `nstd133` (782), and `estd210` (2). It is absent from all of `nstd115`, `nstd99` and
     `nstd107` (66,105 variants), from the `translocation` regions (except one in `nstd133`), and from `nstd133`'s
     inversions and CNVs.
   - Where present, each pair is one `From` and one `To`.
   - `nstd133`'s 245 complex chromosomal rearrangements have up to 114 single-base placements, each `From` or `To`
     (for example `nsv2781790`: 3 `From` and 3 `To` on chromosomes 1 and 3). `From`/`To` alone doesn't say which
     `From` joins which `To`.

   **The strand (`VARIANT_PLACEMENT_STRAND`) gives the join orientation, given which end is `From`** (checked
   2026-11-20). dbVar publishes VCFs with breakend notation for the same `nstd` studies. For every one of their 5,054
   translocation records in `nstd106`, `nstd133` and `nstd107` (GRCh37), the bracket form follows from the strands of
   the two DGVa placements, with no exception. dbVar writes the record at `From`, with `p` the `To` position:

   | `From` strand | `To` strand | breakend at `From` | meaning |
   |---|---|---|---|
   | + | + | `t[p[` | the sequence left of `From` continues right of `To` |
   | + | − | `t]p]` | left of `From`, then the reverse complement of the sequence left of `To` |
   | − | + | `[p[t` | the reverse complement of the sequence right of `To`, then right of `From` |
   | − | − | `]p]t` | the sequence left of `To` continues right of `From` |

   That is, the pair reads as a path along `From` on its strand into `To` on its strand. The mate record at `To` is
   the reverse of the same junction. Limits:
   - **The rule needs to know which end is `From`.** Swapping the two ends changes the result unless both strands
     are equal. Neither the placement ID nor the rank tells it. `nstd106` always has `From` at the lower genomic
     position, and `nstd133` doesn't (390 lower, 391 higher). For `nstd107`, which has no breakpoint order in DGVa,
     dbVar's record is always at the lower position: dbVar seems to take the lower end as `From` when it isn't given.
     The same applies to `nstd99` and `nstd115` (no breakpoint order either): their orientation would rest on that
     convention.
   - **`not reported` and missing strands.** dbVar writes them as `+`. The intrachromosomal translocations of
     `nstd106` (1,211) and `nstd107` (461) have `not reported` on both ends, and those of `nstd99` (475) have no strand
     at all; a few dozen interchromosomal ones have one end missing. Their orientation is unknown, and breakend
     notation has no way to say so: these would stay `<TRA>`, or be written with an assumed orientation and flagged.
   - The complex chromosomal rearrangements of `nstd133` have a strand on every placement, but which `From` joins
     which `To` is still not given.

   So case a has to be recognised from the type (translocations, complex rearrangements) together with single-base
   placements, or for inversions and CNVs from two single-base placements on one chromosome. The breakpoint order is
   worth writing where it exists.

   **Decision (2026-11-20), implemented in `reader.plan_multi_placements()` from `sql/multi_placement_variants.sql`:**
   - IDs: the `BestAvailable` placement keeps the accession (the lowest placement ID if none is best); the others get
     `<accession>_2`, `_3`... in placement ID order. DGVa accessions never contain `_`.
   - `EVENT=<accession>` on every record of the variant; no `DBVARID`. `PLACEMENT_RANK` on each of them.
   - **Translocation given as two single-base placements with both strands:** a breakend pair (`MATEID` both ways,
     `EVENTTYPE=TRA`), orientation from the strands. With no breakpoint order, `From` is the lower position (natural
     chromosome order, then position), as dbVar does.
   - **Translocation with a strand `not reported` or missing:** two `<TRA>` records linked by `EVENT`,
     `EVENTTYPE=TRA`.
   - **Complex chromosomal rearrangement given as single-base placements, each with `From`/`To` and a strand:** one
     single breakend per breakpoint (`From`: `+` → `T.`, `-` → `.T`; `To`: `+` → `.T`, `-` → `T.`), `EVENTTYPE=COMPLEX`,
     declared in each VCF that uses it, as DRAGEN declares its own:
     `##EVENTTYPE=<ID=COMPLEX,Description="Complex chromosomal rearrangement with several breakpoints, ...">`.
     Without the order or a strand: `<COMPLEX_CHROMOSOMAL_REARRANGEMENT>` records linked by `EVENT`.
   - **Any other type given as two single-base placements on one chromosome** (the `nstd133` inversions and CNVs): one
     record from the first breakpoint to the second; the second placement goes on the side list as
     `merged_breakpoint`.
   - **Everything else** (alternative locations, pieces): one record per placement, linked by `EVENT`.
   - Placements with identical coordinates (`estd208`) count once.

   vcf-validator 0.10.2 accepts the breakend pairs, the single breakends, `EVENT`/`EVENTTYPE`/`MATEID` with
   `Number=A`, and the `##EVENTTYPE` line. It warns on a single breakend joined on its right (`.T`), whose ALT doesn't
   start with REF by design.
5. **Key names (decided 2026-11-20).** The names are in the key table above and in `etc/vcf_keys.yaml`. The rules
   that produced them (developer, 2026-11-20):
   - **convertGVFtoVCF's names don't matter.** Only one of the two converters will be used.
   - **Reuse before inventing:** a VCF 4.4 reserved key first (`EVENT`, `EVENTTYPE`, `MATEID`, `AC`,
     `AN`, `SOMATIC`), then a key EVA already writes: the release's `VC` holds the SO accession of the variant class
     (`VC=SO:0000159`), which replaces `TYPESOID`. The release's other keys (`SID`, `RS_VALIDATED`/`SS_VALIDATED`,
     `REMAPPED`, `ALMM`, `ASMM`, `LOE`) describe clustered variants and have no DGVa equivalent.
   - **Variant or analysis:** a property of the whole analysis goes in the metadata (the default placement method, the
     assertion method when uniform, the validation experiment), and can also stay on the variant.
   - **Long, descriptive names** rather than short, ambiguous ones. Exception: `CLNSIG`, the ClinVar/dbVar name, kept
     so that tools recognise it.
   - **No key for a value encoded elsewhere:** `PLACEMENTS` and `SAMPLESET` are dropped; strand and breakpoint order
     are used to build breakends rather than written.

   **`PLACEMENT_ALT_STATUS` is kept** (checked 2026-11-20): it doesn't repeat the rank. Only `nstd19` has it (mouse
   insertions, 10,048 placements: 5,024 submitted on `GCF_000001635.16`, copied to the remapped ones). Every placement
   with a status is `BestAvailable` and the only placement of its variant on the assembly, whether `Known` (4,420) or
   `Ambiguous` (604). `Ambiguous` doesn't follow the chromosome, type (insertion or mobile element insertion), method,
   or the width of the outer window either. The XSD defines it only as the status of the placement "Recipient", and
   its meaning isn't documented further, so the key's description stays literal.

   **Validation (decided 2026-11-20): `VALIDATION_RESULT` only, not the reserved `VALIDATED` flag.** Oracle's optimizer
   statistics on `VARIANT_CALL_VALIDATION` (a frequency histogram over all 2,992,979 rows, analysed 2026-03-11; read
   without scanning the table) give `Pass` 2,523,234, `Fail` 371,029 and `Inconclusive` 98,716: a flag could only say
   `Pass`. The rows cover 2,138,880 calls, so a call can have several results.

   **Experiments can't be assumed one per analysis** (checked 2026-11-20 on the 761 experiments of `DGVA_EXPERIMENT`,
   with `DGVA_METHOD`, `DGVA_ANALYSIS`, `DGVA_DETECTION` and `DGVA_PLATFORM`). Several experiments of one study can be
   the same kind of experiment, so they were compared on their description, not their IDs:

   | Experiment type | Studies with several | All alike: method type + analysis type | + detection | + platform |
   |---|---|---|---|---|
   | Validation | 28 of 59 | 3 | 3 | 1 |
   | Discovery | 73 of 206 | 26 | 24 | 14 |
   | Genotyping | 2 of 2 | 2 | 1 | 1 |

   So 25 of the validated studies use genuinely different validation methods (for example `nstd37`: aCGH,
   karyotyping, SNP array, MLPA, FISH and PCR), and 47 studies different discovery methods. The analysis type alone
   hides this ("Manual observation" covers PCR, FISH and karyotyping); the method type shows it. Each call has exactly
   one discovery experiment (`VARIANT_CALL_EXPERIMENT` has one row per call: 37,536,419), but can have several
   validation results (2,992,979 results for 2,138,880 calls). Whether one call was validated by several different
   methods needs a scan of `VARIANT_CALL_VALIDATION` (outside busy hours); it doesn't change the decision below.

   **Decision (2026-11-20): records reference merged experiment descriptions.**
   - The experiments of a study are merged when their whole description is identical: experiment type, method type,
     analysis type, detection method and platform(s). Each merged description gets an ID, assigned in a fixed order
     (by the lowest submitter experiment ID it covers), and is described once in the metadata, with the submitter
     experiment IDs it merges.
   - **Validation:** `VALIDATION_RESULT=<description ID>:<result>`, one value per result (for example
     `VALIDATION_RESULT=V1:Pass,V2:Fail`).
   - **Discovery (decided 2026-11-20): one EVA calls analysis per discovery description** (the whole description,
     platform included), per submitted assembly. The description lives in the analysis metadata, so the VCF needs no
     discovery key. Genotyping experiments are treated the same way.

   **Calls analyses if there is one per discovery description** (counted 2026-11-20). 586 of the 761 experiments have
   calls (584 discovery, 2 genotyping; checked with one index probe per experiment, no scan), in 206 studies. Today
   there is one calls analysis per study and submitted assembly:

   | Merged on | Calls analyses (per study, before assemblies) | Studies split | Largest |
   |---|---|---|---|
   | (today) | 206 | 0 | 1 |
   | method type + analysis type | 267 | 38 | 7 |
   | + detection method | 316 | 41 | 24 |
   | **+ platform (the whole description: decided)** | 363 | 51 | 30 (`nstd37`), then 27 (`estd59`), 14 (`estd188`) |

   155 studies keep a single calls analysis. The large splits come from platforms: `nstd37`'s 33 experiments have 2
   method types, 2 analysis types and 4 detection methods, but 28 platforms (clinical laboratories' arrays). Not
   counted: the 14 studies submitted on several assemblies, where a description's calls can be on several assemblies
   (`nstd102`: 7 descriptions on 3 assemblies give 12 calls analyses, some with a single call). The exact count over
   all studies is `dgva-extraction checks --step calls_per_experiment` (`calls_analyses.tsv`), run outside busy hours. The regions analyses don't change: regions have no experiment.

6. **The SO terms in use.** `symbolic_alleles.yaml` is keyed by the Exchange types of the XSD.
   `dgva-extraction check-types` compares it with `SEQUENCE_ONTOLOGY_CV` once the `SELECT`-only account exists.
