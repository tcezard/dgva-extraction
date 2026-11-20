"""Configuration, read from a YAML file kept outside the repository.

Example (never commit a file holding real values)::

    dgva:
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

import yaml

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


def load_database_config(config_path):
    """Read the `dgva` section of the configuration file."""
    with open(config_path) as handle:
        config = yaml.safe_load(handle) or {}
    section = config.get('dgva')
    if not isinstance(section, dict):
        raise ConfigurationError(f'{config_path} has no "dgva" section')
    missing = [key for key in ('host', 'port', 'service_name', 'user') if not section.get(key)]
    password = section.get('password') or os.environ.get(PASSWORD_ENVIRONMENT_VARIABLE)
    if not password:
        missing.append(f'password (or the {PASSWORD_ENVIRONMENT_VARIABLE} environment variable)')
    if missing:
        raise ConfigurationError(f'{config_path}: missing {", ".join(missing)} in the "dgva" section')
    return DatabaseConfig(
        host=section['host'],
        port=int(section['port']),
        service_name=section['service_name'],
        user=section['user'],
        password=password,
        schema=section.get('schema', DEFAULT_SCHEMA),
        arraysize=int(section.get('arraysize', DEFAULT_ARRAYSIZE)),
    )
