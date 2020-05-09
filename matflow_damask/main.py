'`matflow_damask.main.py`'

from textwrap import dedent
from pathlib import Path

import numpy as np
from damask_parse import (
    read_geom,
    read_table,
    write_load_case,
    write_geom,
    write_material_config,
)
from damask_parse.utils import get_header

from matflow_damask import (
    input_mapper,
    output_mapper,
    cli_format_mapper,
    func_mapper,
    register_output_file,
)

@output_mapper('microstructure_seeds', 'generate_microstructure_seeds', 'random')
def read_seeds_from_random(path):
    'Parse the file from the `seeds_fromRandom` DAMASK command.'

    header_lns = get_header(path)
    num_header_lns = len(header_lns)

    grid_size = None
    random_seed = None
    for ln in header_lns:
        if ln.startswith('grid'):
            grid_size = [int(j) for j in [i for i in ln.split()][1:][1::2]]
        if ln.startswith('randomSeed'):
            random_seed = int(ln.split()[1])

    data = np.loadtxt(path, skiprows=(num_header_lns + 1), ndmin=2)
    position = data[:, 0:3]
    eulers = data[:, 3:6]

    out = {
        'position': position,
        'eulers': eulers,
        'grid_size': grid_size,
        'random_seed': random_seed,
    }

    return out

@output_mapper('volume_element', 'generate_volume_element', 'random_voronoi')
@output_mapper('volume_element', 'generate_volume_element', 'random_voronoi_from_orientations')
def read_damask_geom(path):
    return read_geom(path)


@input_mapper('load.load', 'simulate_volume_element_loading', 'CP_FFT')
def write_damask_load_case(path, load_case):
    write_load_case(path, load_case)


@input_mapper('geom.geom', 'simulate_volume_element_loading', 'CP_FFT')
@input_mapper('geom.geom', 'visualise_volume_element', 'VTK')
def write_damask_geom(path, volume_element):
    write_geom(volume_element, path)


@input_mapper('material.config', 'simulate_volume_element_loading', 'CP_FFT')
def write_damask_material(path, material_properties, volume_element):
    write_material_config(material_properties, Path(path).parent, volume_element)


@output_mapper('volume_element_response', 'simulate_volume_element_loading', 'CP_FFT')
def read_damask_table(path):
    return read_table(path)


@input_mapper('orientation.seeds', 'generate_volume_element', 'random_voronoi')
def write_microstructure_seeds(path, microstructure_seeds):

    grid_size = microstructure_seeds['grid_size']
    position = microstructure_seeds['position']
    eulers = microstructure_seeds['eulers']

    data = np.hstack([position, eulers, np.arange(1, len(position) + 1)[:, None]])

    header = f"""
        3 header
        grid a {grid_size[0]} b {grid_size[1]} c {grid_size[2]}
        microstructures {len(data)}
        1_pos 2_pos 3_pos 1_euler 2_euler 3_euler microstructure
    """

    fmt = ['%20.16f'] * 6 + ['%10g']
    header = dedent(header).strip()
    np.savetxt(path, data, header=header, comments='', fmt=fmt)

@input_mapper(
    'orientation.seeds',
    'generate_volume_element',
    'random_voronoi_from_orientations'
)
def write_microstructure_new_orientations(path, microstructure_seeds, orientations):

    grid_size = microstructure_seeds['grid_size']
    position = microstructure_seeds['position']
    eulers = orientations['euler_angles']

    data = np.hstack([position, eulers, np.arange(1, len(position) + 1)[:, None]])

    header = f"""
        3 header
        grid a {grid_size[0]} b {grid_size[1]} c {grid_size[2]}
        microstructures {len(data)}
        1_pos 2_pos 3_pos 1_euler 2_euler 3_euler microstructure
    """

    fmt = ['%20.16f'] * 6 + ['%10g']
    header = dedent(header).strip()
    np.savetxt(path, data, header=header, comments='', fmt=fmt)

@cli_format_mapper('size', 'generate_volume_element', 'random_voronoi_from_orientations')
def format_rve_size(size):
    return ' '.join(['{}'.format(i) for i in size])


register_output_file('VTR_file', 'geom.vtr', 'visualise_volume_element', 'VTK')
