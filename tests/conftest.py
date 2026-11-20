import os
import shutil
import subprocess

import pytest

from dgva_extraction.reference import FastaReference, Genome

# Two small chromosomes. chr1 has an IUPAC code (R) at position 61 and lower-case bases at 121-130.
CHR1 = ('ACGTACGTAC' * 6 + 'R' + 'TTGCATGCA' + 'GGCCAATTGG' * 5 + 'acgtacgtac' + 'ACGTTGCA' * 5 + 'TGCA' * 5)
CHR2 = 'CCCCAAAAGGGGTTTT' * 20
ASSEMBLY_REPORT = """# Assembly name:  TEST
# Sequence-Name\tSequence-Role\tAssigned-Molecule\tAssigned-Molecule-Location/Type\tGenBank-Accn\tRelationship\tRefSeq-Accn\tAssembly-Unit\tSequence-Length\tUCSC-style-name
2\tassembled-molecule\t2\tChromosome\tCM000002.1\t=\tNC_000002.1\tPrimary Assembly\t320\tchr2
1\tassembled-molecule\t1\tChromosome\tCM000001.1\t=\tNC_000001.1\tPrimary Assembly\t190\tchr1
"""


def write_fasta(path, sequences, width=60):
    with open(path, 'w') as handle:
        for name, sequence in sequences:
            handle.write(f'>{name} test sequence\n')
            for index in range(0, len(sequence), width):
                handle.write(sequence[index:index + width] + '\n')


@pytest.fixture
def reference_files(tmp_path):
    fasta = tmp_path / 'reference.fa'
    # FASTA order (1, 2) differs from the assembly report order (2, 1) on purpose.
    write_fasta(fasta, [('1', CHR1), ('2', CHR2)])
    report = tmp_path / 'assembly_report.txt'
    report.write_text(ASSEMBLY_REPORT)
    return str(fasta), str(report)


@pytest.fixture
def genome(reference_files):
    fasta_path, report_path = reference_files
    fasta = FastaReference(fasta_path)
    yield Genome(fasta, report_path, 'GCA_TEST.1')
    fasta.close()


def vcf_validator_errors(vcf_path):
    """Run vcf_validator on a VCF: None if it is not installed, else the list of error lines (empty when valid)."""
    validator = shutil.which('vcf_validator')
    if validator is None:
        return None
    output_dir = os.path.dirname(vcf_path)
    result = subprocess.run([validator, '-i', vcf_path, '-r', 'text', '-o', output_dir],
                            capture_output=True, text=True, check=False)
    errors = []
    for name in os.listdir(output_dir):
        if name.startswith(os.path.basename(vcf_path) + '.errors') and name.endswith('.txt'):
            with open(os.path.join(output_dir, name)) as handle:
                errors.extend(line.strip() for line in handle if line.startswith('Line '))
    if result.returncode != 0 and not errors:
        errors.append(result.stdout + result.stderr)
    return errors
