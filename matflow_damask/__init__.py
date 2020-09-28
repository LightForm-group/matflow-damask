'`matflow_damask.__init__.py`'

from functools import partial

from matflow_damask._version import __version__

from matflow.extensions import (
    input_mapper,
    output_mapper,
    cli_format_mapper,
    sources_mapper,
    software_versions,
    register_output_file,
    func_mapper,
)

SOFTWARE = 'DAMASK'

input_mapper = partial(input_mapper, software=SOFTWARE)
output_mapper = partial(output_mapper, software=SOFTWARE)
cli_format_mapper = partial(cli_format_mapper, software=SOFTWARE)
software_versions = partial(software_versions, software=SOFTWARE)
register_output_file = partial(register_output_file, software=SOFTWARE)
sources_mapper = partial(sources_mapper, software=SOFTWARE)
func_mapper = partial(func_mapper, software=SOFTWARE)

from matflow_damask import main
