# VCF representation of the DGVa data (sprint 3, ticket 10)

**Status: proposed, for review by a second developer.** The writer implements these rules, so they can be tested and
checked with vcf-validator, but nothing here is agreed until the review. Changing a rule after the studies are
submitted means resubmitting them.

The keys themselves are in [`dgva_extraction/etc/vcf_keys.yaml`](../dgva_extraction/etc/vcf_keys.yaml), and the
symbolic allele for each variant type in [`dgva_extraction/etc/symbolic_alleles.yaml`](../dgva_extraction/etc/symbolic_alleles.yaml).
Those two files are the single source: the writer builds the header from them, checks every value against them, and
refuses any value they don't list. Each key has its `Number`, `Type`, `Description`, database source, whether VCF
4.4 reserves it, its status and the sprint whose reader fills it.

The proposals were checked with vcf-validator 0.10.2 on VCFs written by the writer (`tests/test_vcf_writer.py`,
`tests/test_cli.py`). They have **not** been run through eva-sub-cli validation yet (ticket 10, step 4).

## Keys

| Key | Field | Number | Type | Reserved (4.4) | Source | Status |
|---|---|---|---|---|---|---|
| `END` | INFO | 1 | Integer | yes | stop (or outer/inner stop); `POS` for insertions | agreed |
| `SVLEN` | INFO | A | Integer | yes | stop − start + 1; `.` for insertions | agreed |
| `IMPRECISE` | INFO | 0 | Flag | yes | any CI bound not 0 | agreed |
| `CIPOS`, `CIEND` | INFO | 2 | Integer | yes | `CIPOS*`/`CIEND*`, else `CI_START`/`CI_STOP`, else inner/outer | agreed |
| `SVCLAIM` | INFO | A | String | yes | type (copy number gain/loss: `D`); otherwise **open**, see below | proposed |
| `TYPESOID` | INFO | 1 | String | no | SO ID of the type (name from the legacy VCF exporter) | proposed |
| `SUBMITTED_TYPE` | INFO | 1 | String | no | submitter type, where it differs beyond case and separators | proposed (sprint 4) |
| `PARENT` | INFO | . | String | no | `LINKED_VARIANT`: call → region and region → region | proposed |
| `ASSERTM` | INFO | 1 | String | no | region assertion method (name from the legacy VCF exporter) | proposed |
| `PLACEMETHOD` | INFO | 1 | String | no | placement method, when not `Submitted genomic` | proposed |
| `AF` | INFO | A | Float | yes | allele frequency (incl. the `estd219` load table); `.` on every line of a VCF with no samples | **agreed (2026-11-19)** |
| `CN` | INFO | A | Float | yes | copy number (incl. the `estd214`/`estd219` load tables) | proposed (sprint 4) |
| `RUC`, `RUS` | INFO | A, . | Float, String | yes | repeat counts and motifs (`nstd128`) | proposed (sprint 4) |
| `HGVS` | INFO | . | String | no | `VARIANT_PLACEMENT_HGVS` | proposed (sprint 4) |
| `PLACRANK`, `PLACEMENTS` | INFO | 1 | String | no | `VAR_PLAC_PLACEMENT_RANK`, `VAR_PLAC_PLACS_PER_ASSEM` | proposed (sprint 4) |
| `CYTOBAND`, `MUTORDER`, `ALTSTATUS` | INFO | 1, ., 1 | String, Integer, String | no | gap rows 105, 52/114, 93 | proposed (sprint 4) |
| `SAMPLESET`, `VALIDMETH`, `VALIDRES` | INFO | . | String | no | from the legacy VCF exporter; gap rows 88, 90, 91 | proposed (sprint 4) |
| `SVSEQ` | INFO | 1 | String | no | a sequence that is not an allele (gap row 89) | proposed (sprint 4) |
| `GT` | FORMAT | 1 | String | yes | zygosity of the call's sample and the species' ploidy (rule 4) | **agreed (2026-11-19)**, read in sprint 4 |
| `CN` | FORMAT | 1 | Float | yes | copy number of the call's sample | proposed (sprint 4) |

Header lines and INFO keys are written in the order of `vcf_keys.yaml`, and only the keys that a VCF uses are in its
header.

## Record rules

1. **One VCF record per DGVa accession and placement**, with the accession (`esv`/`nsv`/`essv`/`nssv`) in the ID
   column. Records at the same position are never merged.
2. **Regions and calls are both written, in separate VCFs for separate EVA analyses** (agreed 2026-11-19). The regions
   are DGVa's calls made on top of the submitted calls, not the raw submitted data, so they get their own analysis.
   A call points at its region(s) with `PARENT`, and so does a region nested in another region (`estd1`, `estd20`,
   `estd59`). A call with no region (W1: `estd236`, `estd237`) has no `PARENT`. A withdrawn parent is left out of
   `PARENT`. See "File layout" below.
3. **A value with no target fails the study.** The writer refuses a record whose type has no symbolic allele, or that
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
5. **Placement method.** `PLACEMETHOD` is written on each record whose method is not `Submitted genomic`
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

Per study and submitted assembly:

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
    sampleset is written in its samples' files only, and `SAMPLESET` in INFO keeps the link.
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
| Sequences that are not an allele (row 89) | INFO `SVSEQ` |
| `PLACEMENTS` always, or only when not `Single` (row 116) | **Always**, when the database has a value. Leaving out a default hides the difference between "Single" and "no value". |
| `SUBMITTED_TYPE` (rows 86, 108) | Only where it differs from the Exchange type beyond case and `_`/`-` |
| Placement method | Per line for non-default methods; header line for the null-method studies (rule 5) |

## Decisions still open

These came out of checking the writer's output with vcf-validator, or have no answer in the gap table. Each one
currently makes the writer refuse the record, so a study it affects fails until it is decided.

1. **VCF version: 4.4 (proposed).** 4.4 reserves the SV keys the gap table uses (`RUC`, `RUS`, `SVCLAIM`, `CNV:TR`). It
   also brings rules that vcf-validator enforces:
   - `SVLEN` on every symbolic allele, including `<INS>`, where it is written `.`;
   - `SVCLAIM` on every `<DEL>` and `<DUP>`;
   - `RUS` or `RUL` on every `<CNV:TR>`;
   - the repeat keys only on `<CNV:TR>`.

   Writing 4.3 would avoid the `SVCLAIM` question but lose the reserved repeat keys.
2. **`SVCLAIM` for `<DEL>`/`<DUP>`.** A copy number gain or loss is by definition an abundance claim (`D`). For
   `deletion`, `duplication`, `tandem duplication` and the mobile element deletions, the type doesn't say.
   - Proposal: derive it from the analysis type of the call's experiment (`VARIANT_CALL_EXPERIMENT` →
     `EXPERIMENT_ANALYSIS` → `DGVA_ANALYSIS.ANALYSIS_TYPE`).
     - `Probe signal intensity`, `Read depth`, `SNP genotyping analysis` and `MCD analysis` give `D`.
     - `Paired-end mapping`, `Split read mapping`, `Split read and paired-end mapping`, `Sequence alignment` and the
       assembly-based analyses give `J`.
     - `Read depth and paired-end mapping` gives `DJ`.
     - `Curated`, `Merging`, `Multiple`, `Not provided` and `Other` remain open, as do regions, which have no
       experiment of their own. Their claim could come from their calls.
   - Until this is decided, a `deletion` or `duplication` record is refused (`missing_svclaim`).
3. **Types with no VCF symbolic allele.**
   - `sequence alteration`: `estd1` alone has 25,193 calls.
   - `complex substitution`, `indel`/`delins` without a sequence, `copy-neutral loss of heterozygosity`.
   - The translocations, which need breakend (`BND`) notation with the mate's position and orientation. That includes
     the `nstd83` translocation regions, which the reader positions from their supporting calls.

   They are refused (`no_symbolic_allele`). Options: a non-standard symbolic allele declared in the header, `<CNV>` with
   the type in `TYPESOID`, or leaving them out of the VCF and listing them in the DGVa JSON.
4. **Duplicate IDs.** A variant with two submitted placements on one assembly (72,857 `Submitted genomic` placements
   with rank `NotBest`) gives two records with the same ID. vcf-validator 0.10.2 accepts this, but VCF 4.4 says no
   identifier should be present in more than one data record. Options: accept it (one record per placement, as
   proposed), or write only the `BestAvailable` placement and list the others in the DGVa JSON.
5. **Are the new key names right?** `TYPESOID`, `ASSERTM`, `PLACEMETHOD`, `PARENT`, `SAMPLESET`, `VALIDMETH` and
   `VALIDRES` reuse the legacy DGVa VCF exporter's names where it had one. The rest (`HGVS`, `PLACRANK`, `PLACEMENTS`,
   `CYTOBAND`, `MUTORDER`, `SUBMITTED_TYPE`, `ALTSTATUS`, `SVSEQ`) are new.
6. **The SO terms in use.** `symbolic_alleles.yaml` is keyed by the Exchange types of the XSD.
   `dgva-extraction check-types` compares it with `SEQUENCE_ONTOLOGY_CV` once the `SELECT`-only account exists.
