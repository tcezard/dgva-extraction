"""dgva-extraction: DGVa structural variants from the DGVa database to VCF, for the migration into EVA."""
try:
    from dgva_extraction._version import version as __version__
except ImportError:  # not installed
    __version__ = '0.0.0'
