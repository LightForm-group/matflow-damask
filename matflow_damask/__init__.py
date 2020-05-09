'`matflow_damask.__init__.py`'

from functools import partial

from matflow_damask._version import __version__

from matflow import (
    input_mapper,
    output_mapper,
    cli_format_mapper,
    func_mapper,
    register_output_file,
)

input_mapper = partial(input_mapper, software='damask')
output_mapper = partial(output_mapper, software='damask')
cli_format_mapper = partial(cli_format_mapper, software='damask')
register_output_file = partial(register_output_file, software='damask')

from matflow_damask import main
