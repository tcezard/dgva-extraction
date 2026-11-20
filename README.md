# dgva-extraction

Extracts the structural variants of the Database of Genomic Variants archive (DGVa) from its Oracle database into
VCF, for the migration of DGVa into the European Variation Archive (EVA). Part of EVA-2269.

DGVa's GVF files lack fields on millions of submitted variants, so the migration reads the database instead (option B
of the DGVa migration plan). This repository holds that reader and a VCF writer written from scratch. The GVF-based
converter, [convertGVFtoVCF](https://github.com/EBIvariation/convertGVFtoVCF), is a reference only: nothing is
imported from it.

## What it does

- **Reads** one study on one submitted assembly, streaming its submitted placements in database order. It follows the
  rules measured in the gap analysis:
  - submitted and remapped placements are told apart, and only submitted ones are written;
  - withdrawn variants are left out;
  - the W items and the confidence-interval sources are handled as in the reference document.
- **Writes** one VCF record per DGVa accession and placement, with the accession in the ID column. The output is
  deterministic: the same records give the same bytes.
- **Fails the study** rather than dropping data: a type with no SO term or no symbolic allele, a value with no VCF
  key, or rows that disagree all stop the study.
- **Lists** what is in no VCF (withdrawn, unplaced, remapped-only and no-position variants; placements with no start or
  chromosome; duplicate placements) in a side list per study.

The VCF representation (keys, record rules, open decisions) is in
[`docs/vcf-representation.md`](docs/vcf-representation.md). The keys are in
[`dgva_extraction/etc/vcf_keys.yaml`](dgva_extraction/etc/vcf_keys.yaml).

## Database safety

dgva-extraction **runs only with a `SELECT`-only account** on the `DGVA` schema.

- At start it checks the session's privileges, and stops if it finds any write privilege.
- It never commits, runs no DDL or PL/SQL, and rolls back before closing.
- SQL is kept in [`dgva_extraction/sql/`](dgva_extraction/sql) and takes bind variables only.
- Extracting a large study (`estd214`, `estd219`) and running `list-studies` over all studies read tens of millions of
  rows. Run them outside busy hours.

## Install

```bash
python -m venv .venv && . .venv/bin/activate
pip install -e '.[test]'
pytest
```

The writer's tests run [vcf-validator](https://github.com/EBIvariation/vcf-validator) on the VCFs they write when
`vcf_validator` is on the `PATH`. CI installs it.

## Configuration

Credentials are read from a YAML file **outside the repository** (never commit one):

```yaml
dgva:
  host: <database host>
  port: 1521
  service_name: <service>
  user: <SELECT-only account>
  password: <password>   # or leave it out and set DGVA_PASSWORD
  schema: DGVA
  arraysize: 10000
```

## Usage

```bash
# Check that the configured account cannot write
dgva-extraction --config ~/dgva/config.yaml check-privileges

# Studies and assemblies, with their counts per placement class (what to convert)
dgva-extraction --config ~/dgva/config.yaml list-studies --output studies.tsv

# Variant types in use and their symbolic allele
dgva-extraction --config ~/dgva/config.yaml check-types

# One study on one submitted assembly: a regions VCF and a calls VCF (separate EVA analyses), statistics
# (with run time and peak memory) and side list
dgva-extraction --config ~/dgva/config.yaml extract --study estd1 --assembly GCA_000001405.1 \
    --fasta GCA_000001405.1.fa --assembly-report GCA_000001405.1_assembly_report.txt --output-dir output/

# Variants of a study that are in no VCF
dgva-extraction --config ~/dgva/config.yaml side-list --study estd1 --output estd1.side_list.tsv
```

## Licence

Apache 2.0. See [LICENSE](LICENSE).
