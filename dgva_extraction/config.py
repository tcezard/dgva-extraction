"""Configuration, through the shared `cfg` singleton of ebi-eva-common-pyutils, as in eva-submission.

The YAML file is kept outside the repository. It is found, in order, at the path given on the command line, at the
path in the DGVA_EXTRACTION_CONFIG environment variable, or at ~/.dgva_extraction_config.yml. Its DGVA section uses
the same keys as convertGVFtoVCF's config, so one file can serve both (never commit one holding real values)::

    DGVA:
      host: <database host>
      port: 1521
      service_name: <service>
      user: <SELECT-only account>
      password: <password>        # or leave it out and set DGVA_PASSWORD in the environment
      schema: DGVA                 # the schema that holds the DGVa tables
      arraysize: 10000             # rows fetched per round trip
"""
import os
from dataclasses import dataclass

from ebi_eva_common_pyutils.config import cfg

CONFIG_ENVIRONMENT_VARIABLE = 'DGVA_EXTRACTION_CONFIG'
DEFAULT_CONFIG_PATH = '~/.dgva_extraction_config.yml'
CONFIG_SECTION = 'DGVA'
DEFAULT_ARRAYSIZE = 10000
DEFAULT_SCHEMA = 'DGVA'
PASSWORD_ENVIRONMENT_VARIABLE = 'DGVA_PASSWORD'


class ConfigurationError(ValueError):
    pass


@dataclass(frozen=True)
class DatabaseConfig:
    host: str
    port: int
    service_name: str
    user: str
    password: str
    schema: str = DEFAULT_SCHEMA
    arraysize: int = DEFAULT_ARRAYSIZE

    def __repr__(self):
        # Keep the password out of logs and tracebacks.
        return (f'DatabaseConfig(host={self.host!r}, port={self.port!r}, service_name={self.service_name!r}, '
                f'user={self.user!r}, password=***, schema={self.schema!r}, arraysize={self.arraysize!r})')


def load_config(*args):
    """Load the configuration into `cfg` from the first path that exists: the ones given, then
    $DGVA_EXTRACTION_CONFIG, then ~/.dgva_extraction_config.yml."""
    cfg.load_config_file(
        *args,
        os.getenv(CONFIG_ENVIRONMENT_VARIABLE),
        os.path.expanduser(DEFAULT_CONFIG_PATH),
    )


def database_config():
    """The database connection details from the DGVA section of the loaded configuration."""
    section = cfg.get(CONFIG_SECTION)
    if not isinstance(section, dict):
        raise ConfigurationError(f'{cfg.config_file} has no {CONFIG_SECTION} section')
    missing = [key for key in ('host', 'port', 'service_name', 'user') if not section.get(key)]
    password = section.get('password') or os.environ.get(PASSWORD_ENVIRONMENT_VARIABLE)
    if not password:
        missing.append(f'password (or the {PASSWORD_ENVIRONMENT_VARIABLE} environment variable)')
    if missing:
        raise ConfigurationError(f'{cfg.config_file}: missing {", ".join(missing)} in the {CONFIG_SECTION} section')
    return DatabaseConfig(
        host=section['host'],
        port=int(section['port']),
        service_name=section['service_name'],
        user=section['user'],
        password=password,
        schema=section.get('schema', DEFAULT_SCHEMA),
        arraysize=int(section.get('arraysize', DEFAULT_ARRAYSIZE)),
    )
