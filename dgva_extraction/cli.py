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

from dgva_extraction import __version__
from dgva_extraction.config import load_database_config
from dgva_extraction.database import DgvaDatabase
from dgva_extraction.keys import SymbolicAlleles
from dgva_extraction.reader import StudyReader, study_side_list, write_side_list
from dgva_extraction.reference import FastaReference, Genome
from dgva_extraction.vcf_writer import StudyVcfWriter, VcfLineBuilder

logger = logging.getLogger('dgva_extraction')


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
    return DgvaDatabase.connect(load_database_config(args.config))


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
    parser.add_argument('--config', help='YAML configuration with a "dgva" section, kept outside the repository')
    commands = parser.add_subparsers(dest='command', required=True)

    commands.add_parser('check-privileges', help='Connect and check that the session cannot write')

    command = commands.add_parser('list-studies', help='Studies and assemblies with their counts per class')
    command.add_argument('--study', default='ALL', help='One study accession, or ALL (default; slow)')
    command.add_argument('--output', required=True)

    commands.add_parser('check-types', help='Variant types in use, and their symbolic allele')

    command = commands.add_parser('side-list', help='Variants of a study that are in no VCF')
    command.add_argument('--study', required=True)
    command.add_argument('--output', required=True)

    command = commands.add_parser('extract', help='Write the VCF of one study on one submitted assembly')
    command.add_argument('--study', required=True)
    command.add_argument('--assembly', required=True, help='Assembly accession, as in VARIANT_PLACEMENT')
    command.add_argument('--fasta', required=True, help='Uncompressed FASTA of the assembly (indexed if needed)')
    command.add_argument('--assembly-report', help='NCBI assembly report: contig order and name aliases')
    command.add_argument('--output-dir', required=True)
    return parser


COMMANDS = {'check-privileges': check_privileges, 'list-studies': list_studies, 'check-types': check_types,
            'side-list': side_list, 'extract': extract}


def main(argv=None):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(name)s: %(message)s')
    args = build_parser().parse_args(argv)
    if not args.config:
        build_parser().error('--config is required')
    return COMMANDS[args.command](args)


if __name__ == '__main__':
    sys.exit(main())
