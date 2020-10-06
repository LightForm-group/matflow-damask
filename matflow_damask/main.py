'`matflow_damask.main.py`'

import json
from textwrap import dedent
from pathlib import Path

import hickle
import numpy as np
import pkg_resources
from damask_parse import (
    read_HDF5_file,
    write_load_case,
    write_geom,
    write_material,
    write_numerics,
    geom_to_volume_element,
)
from damask_parse.utils import (
    get_header_lines,
    parse_damask_spectral_version_info,
    add_volume_element_missing_texture,
    volume_element_from_2D_microstructure,
)
from damask_parse import __version__ as damask_parse_version
from matflow.scripting import get_wrapper_script

from matflow_damask import (
    input_mapper,
    output_mapper,
    cli_format_mapper,
    sources_mapper,
    func_mapper,
    register_output_file,
    software_versions,
)


@input_mapper(
    input_file='orientation_coordinate_system.json',
    task='generate_volume_element',
    method='random_voronoi'
)
def write_orientation_coordinate_system_from_seeds(path, microstructure_seeds):
    with Path(path).open('w') as handle:
        json.dump(
            microstructure_seeds['orientations']['orientation_coordinate_system'],
            handle
        )


@input_mapper(
    input_file='orientation_coordinate_system.json',
    task='generate_volume_element',
    method='random_voronoi_from_orientations'
)
def write_orientation_coordinate_system_from_orientations(path, orientations):
    with Path(path).open('w') as handle:
        json.dump(orientations['orientation_coordinate_system'], handle)


def read_orientation_coordinate_system(path):
    with Path(path).open('r') as handle:
        return json.load(handle)


@output_mapper('microstructure_seeds', 'generate_microstructure_seeds', 'random')
def read_seeds_from_random(seeds_path, orientation_coordinate_system, phase_label):
    'Parse the file from the `seeds_fromRandom` DAMASK command.'

    header_lns = get_header_lines(seeds_path)
    num_header_lns = len(header_lns)

    grid_size = None
    random_seed = None
    for ln in header_lns:
        if ln.startswith('grid'):
            grid_size = [int(j) for j in [i for i in ln.split()][1:][1::2]]
        if ln.startswith('randomSeed'):
            try:
                random_seed = int(ln.split()[1])
            except ValueError:
                # Random seed set to "None" in seeds_fromRandom 2.0.3-1097-ga7fca4df
                pass

    data = np.loadtxt(seeds_path, skiprows=(num_header_lns + 1), ndmin=2)
    position = data[:, 0:3]
    eulers = data[:, 3:6]

    microstructure_seeds = {
        'position': position,
        'orientations': {
            'type': 'euler',
            'euler_angles': eulers,
            'orientation_coordinate_system': orientation_coordinate_system,
        },
        'grid_size': grid_size,
        'random_seed': random_seed,
        'phase_label': phase_label,
    }

    return microstructure_seeds


@output_mapper('volume_element', 'generate_volume_element', 'random_voronoi')
@output_mapper('volume_element', 'generate_volume_element', 'random_voronoi_from_orientations')
def read_damask_geom(geom_path, ori_coord_system_path, phase_label_path, homog_label_path,
                     model_coordinate_system, buffer_phase_label):

    # TODO: coord systems:
    # volume_element['model_coordinate_system'] = model_coordinate_system
    # ori_coord_sys = read_orientation_coordinate_system(ori_coord_system_path)
    # volume_element['orientation_coordinate_system'] = ori_coord_sys

    with Path(phase_label_path).open('r') as handle:
        phase_labels = [handle.read().strip()]

    with Path(homog_label_path).open('r') as handle:
        homog_label = handle.read().strip()

    if buffer_phase_label:
        phase_labels.append(buffer_phase_label)

    volume_element = geom_to_volume_element(
        geom_path,
        phase_labels=phase_labels,
        homog_label=homog_label,
    )

    return volume_element


@input_mapper('load.load', 'simulate_volume_element_loading', 'CP_FFT')
def write_damask_load_case(path, load_case):
    write_load_case(path, load_case)


@input_mapper('geom.geom', 'simulate_volume_element_loading', 'CP_FFT')
@input_mapper('geom.geom', 'visualise_volume_element', 'VTK')
def write_damask_geom(path, volume_element):
    write_geom(volume_element, path)


@input_mapper('material.config', 'simulate_volume_element_loading', 'CP_FFT')
def write_damask_material(path, homogenization_schemes, volume_element,
                          single_crystal_parameters, phases, texture_alignment_method):

    # TODO: sort out texture alignment

    # Merge single-crystal properties into phases:
    for phase_label in phases.keys():
        SC_params_name = phases[phase_label].pop('single_crystal_parameters', None)
        if SC_params_name:
            phases[phase_label].update(**single_crystal_parameters[SC_params_name])

    write_material(
        homog_schemes=homogenization_schemes,
        phases=phases,
        volume_element=volume_element,
        dir_path=Path(path).parent,
    )


@input_mapper('numerics.config', 'simulate_volume_element_loading', 'CP_FFT')
def write_damask_numerics(path, numerics):
    if numerics:
        write_numerics(Path(path).parent, numerics)


@output_mapper('volume_element_response', 'simulate_volume_element_loading', 'CP_FFT')
def read_damask_hdf5_file(hdf5_path, incremental_data, operations=None):
    return read_HDF5_file(hdf5_path, incremental_data, operations=operations)


@input_mapper('phase_label.txt', 'generate_volume_element', 'random_voronoi')
@input_mapper('phase_label.txt', 'generate_volume_element', 'random_voronoi_from_orientations')
def write_phase_label(path, microstructure_seeds):
    with Path(path).open('w') as handle:
        handle.write(f'{microstructure_seeds["phase_label"]}\n')


@input_mapper('homog_label.txt', 'generate_volume_element', 'random_voronoi')
@input_mapper('homog_label.txt', 'generate_volume_element', 'random_voronoi_from_orientations')
def write_homog_label(path, homog_label):
    with Path(path).open('w') as handle:
        handle.write(f'{homog_label}\n')


@input_mapper('orientation.seeds', 'generate_volume_element', 'random_voronoi')
def write_microstructure_seeds(path, microstructure_seeds):

    grid_size = microstructure_seeds['grid_size']
    position = microstructure_seeds['position']
    eulers = microstructure_seeds['orientations']['euler_angles']

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


@func_mapper(task='generate_volume_element', method='extrusion')
def volume_element_from_microstructure_image(microstructure_image, depth, image_axes,
                                             phase_label, homog_label):
    out = {
        'volume_element': volume_element_from_2D_microstructure(
            microstructure_image,
            phase_label,
            homog_label,
            depth,
            image_axes,
        )
    }
    return out


@func_mapper(task='visualise_volume_element', method='VTK')
def visualise_volume_element(volume_element):
    from damask import Geom
    geom_obj = Geom.from_file('geom.geom')
    geom_obj.to_vtr('geom.vtr')


@func_mapper(task='generate_volume_element', method='random_voronoi_2')
def generate_volume_element(microstructure_seeds, size, homog_label, scale_morphology,
                            buffer_phase_size, buffer_phase_label):
    out = generate_volume_element_random_voronoi(
        microstructure_seeds,
        size,
        homog_label,
        scale_morphology,
        buffer_phase_size,
        buffer_phase_label,
        orientations=None,
    )
    return out


@func_mapper(task='generate_volume_element', method='random_voronoi_from_orientations_2')
def generate_volume_element(microstructure_seeds, size, homog_label, scale_morphology,
                            buffer_phase_size, buffer_phase_label, orientations):
    out = generate_volume_element_random_voronoi(
        microstructure_seeds,
        size,
        homog_label,
        scale_morphology,
        buffer_phase_size,
        buffer_phase_label,
        orientations=orientations,
    )
    return out


def generate_volume_element_random_voronoi(microstructure_seeds, size, homog_label,
                                           scale_morphology, buffer_phase_size,
                                           buffer_phase_label, orientations=None):
    from damask import Geom

    geom_obj = Geom.from_Voronoi_tessellation(
        grid=np.array(microstructure_seeds['grid_size']),
        size=np.array(size),
        seeds=np.array(microstructure_seeds['position']),
    )

    if scale_morphology is not None:
        # scale morphology: keep the same "elements per volume", but scale morphology

        scale_morphology = np.array(scale_morphology)
        original_size = geom_obj.get_size()
        original_grid = geom_obj.get_grid()

        new_size = original_size * scale_morphology
        new_grid = original_grid * scale_morphology

        geom_scaled = geom_obj.scale(new_grid)
        geom_scaled.set_size(new_size)

        geom_obj = geom_scaled

    phase_labels = [microstructure_seeds['phase_label']]
    if buffer_phase_size is not None:

        original_grid = geom_obj.get_grid()
        original_size = geom_obj.get_size()

        new_grid = original_grid + np.array(buffer_phase_size)
        new_size = original_size * (new_grid / original_grid)

        geom_canvased = geom_obj.canvas(grid=new_grid)
        geom_canvased.set_size(new_size)

        geom_obj = geom_canvased
        phase_labels.append(buffer_phase_label)

    geom_obj.to_file('geom.geom')

    if orientations is None:
        oris = microstructure_seeds['orientations']['euler_angles']
    else:
        oris = orientations['euler_angles']

    volume_element = geom_to_volume_element(
        'geom.geom',
        phase_labels=phase_labels,
        homog_label=homog_label,
        orientations=oris,
    )

    return {'volume_element': volume_element}


@software_versions()
def get_versions(executable='DAMASK_spectral'):
    'Get versions of pertinent software associated with this extension.'

    out = {
        'DAMASK_spectral': parse_damask_spectral_version_info(executable=executable),
        'damask (Python)': {'version': pkg_resources.get_distribution('damask').version},
        'damask-parse (Python)': {'version': damask_parse_version},
    }
    return out


register_output_file('VTR_file', 'geom.vtr', 'visualise_volume_element', 'VTK')
