'`matflow_damask.main.py`'

import copy
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
    volume_element_from_2D_microstructure,
    add_volume_element_buffer_zones,
    validate_orientations,
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
from matflow_damask.utils import get_by_path, set_by_path


def read_orientation_coordinate_system(path):
    with Path(path).open('r') as handle:
        return json.load(handle)


@func_mapper(task='generate_microstructure_seeds', method='random')
def seeds_from_random(size, num_grains, phase_label, grid_size=None,
                      orientation_coordinate_system=None):
    from damask import seeds
    from damask import Rotation

    size = np.array(size)
    grid_size = grid_size and np.array(grid_size)

    position = seeds.from_random(size, num_grains, cells=grid_size)
    rotation = Rotation.from_random(shape=(num_grains,))

    out = {
        'microstructure_seeds': {
            'position': position,
            'orientations': {
                'type': 'quat',
                'quaternions': rotation.quaternion,
                'orientation_coordinate_system': orientation_coordinate_system,
                'unit_cell_alignment': {
                    'x': 'a',
                    'z': 'c',
                },
                'P': -1,
            },
            'size': size,
            'random_seed': None,
            'phase_label': phase_label,
        }
    }
    return out


@func_mapper(task='sample_texture', method='from_random')
def orientations_from_random(num_orientations,
                             orientation_coordinate_system=None):
    from damask import Rotation

    rotation = Rotation.from_random(shape=(num_orientations,))

    out = {
        'orientations': {
            'type': 'quat',
            'quaternions': rotation.quaternion,
            'orientation_coordinate_system': orientation_coordinate_system,
            'unit_cell_alignment': {
                'x': 'a',
                'z': 'c',
            },
            'P': -1,
        },
    }
    return out


@input_mapper('load.yaml', 'simulate_volume_element_loading', 'CP_FFT')
@input_mapper('load.yaml', 'simulate_orientrations_loading', 'Taylor')
def write_damask_load_case(path, load_case):
    path = Path(path)
    write_load_case(path.parent, load_case, name=path.name)


@input_mapper('geom.vtr', 'simulate_volume_element_loading', 'CP_FFT')
def write_damask_geom(path, volume_element):
    path = Path(path)
    write_geom(path.parent, volume_element, name=path.name)


@input_mapper('material.yaml', 'simulate_volume_element_loading', 'CP_FFT')
def write_damask_material(path, homogenization_schemes, volume_element,
                          single_crystal_parameters,
                          single_crystal_parameter_perturbation, phases,
                          texture_alignment_method):

    # TODO: sort out texture alignment

    # Apply a perturbation to a specific single-crystal parameter:
    if (
        single_crystal_parameters and
        single_crystal_parameter_perturbation and
        single_crystal_parameter_perturbation['perturbation']
    ):
        single_crystal_parameters = copy.deepcopy(single_crystal_parameters)
        scale_factor = (1 + single_crystal_parameter_perturbation['perturbation'])
        address = single_crystal_parameter_perturbation.get('address')
        set_by_path(
            single_crystal_parameters,
            address,
            get_by_path(single_crystal_parameters, address) * scale_factor,
        )
    phases = copy.deepcopy(phases)

    # Merge single-crystal properties into phases:
    for phase_label in phases.keys():

        SC_params_name = phases[phase_label]['mechanical']['plastic'].pop(
            'single_crystal_parameters',
            None
        )
        if SC_params_name:
            SC_params = single_crystal_parameters[SC_params_name]
            phases[phase_label]['mechanical']['plastic'].update(**SC_params)

    path = Path(path)
    write_material(
        homog_schemes=homogenization_schemes,
        phases=phases,
        volume_element=volume_element,
        dir_path=path.parent,
        name=path.name,
    )


@input_mapper('material.yaml', 'simulate_orientrations_loading', 'Taylor')
def write_damask_taylor_material(path, orientations, phases):
    # Convert orientations to quats and check size
    orientations = validate_orientations(orientations)
    num_oris = orientations.get('quaternions').shape[0]
    if num_oris % 8 != 0:
        msg = ('Number of orientations must be a multiple of 8.')
        raise ValueError(msg)
    num_oris_per_mat = num_oris // 8

    # Check only 1 phase and get name
    if len(phases) != 1:
        msg = ('Only one phase should be specified.')
        raise ValueError(msg)
    pahse_name = list(phases.keys())[0]

    volume_element = {
        'element_material_idx': np.arange(8).reshape([2, 2, 2]),
        'grid_size': np.array([2, 2, 2]),
        'size': [1.0, 1.0, 1.0],

        'constituent_material_idx': np.arange(8).repeat(num_oris_per_mat),
        'constituent_material_fraction': np.full(num_oris, 1. / num_oris_per_mat),
        'constituent_phase_label': np.full(num_oris, pahse_name),
        'constituent_orientation_idx': np.arange(num_oris),
        'material_homog': np.full(8, 'Taylor'),
        'orientations': orientations,
    }
    homogenization_schemes = {'Taylor': {
        'mechanical': {'type': 'isostrain'},
        'N_constituents': num_oris_per_mat,
    }}

    path = Path(path)
    write_geom(path.parent, volume_element, name='geom.vtr')
    write_material(
        homog_schemes=homogenization_schemes,
        phases=phases,
        volume_element=volume_element,
        dir_path=path.parent,
        name=path.name,
    )


@input_mapper('numerics.yaml', 'simulate_volume_element_loading', 'CP_FFT')
@input_mapper('numerics.yaml', 'simulate_orientrations_loading', 'Taylor')
def write_damask_numerics(path, numerics):
    if numerics:
        path = Path(path)
        write_numerics(path.parent, numerics, name=path.name)


@output_mapper('volume_element_response', 'simulate_volume_element_loading', 'CP_FFT')
@output_mapper('orientrations_response', 'simulate_orientrations_loading', 'Taylor')
def read_damask_hdf5_file(hdf5_path, incremental_data=None, field_data=None,
                          grain_data=None, operations=None):
    return read_HDF5_file(hdf5_path, incremental_data=incremental_data,
                          field_data=field_data, grain_data=grain_data,
                          operations=operations)


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


@func_mapper(task='modify_volume_element', method='add_buffer_zones')
def modify_volume_element_add_buffer_zones(volume_element, buffer_sizes,
                                           phase_ids, phase_labels, homog_label, order):
    out = {
        'volume_element': add_volume_element_buffer_zones(
            volume_element, buffer_sizes, phase_ids, phase_labels, homog_label, order
        )
    }
    return out


@func_mapper(task='visualise_volume_element', method='VTK')
def visualise_volume_element(volume_element):
    path = Path('geom.vtr')
    write_geom(path.parent, volume_element, name=path.name)


@func_mapper(task='generate_volume_element', method='random_voronoi')
def generate_volume_element_random_voronoi_2(microstructure_seeds, grid_size, homog_label,
                                             scale_morphology, scale_update_size,
                                             buffer_phase_size, buffer_phase_label):
    out = generate_volume_element_random_voronoi(
        microstructure_seeds,
        grid_size,
        homog_label,
        scale_morphology,
        scale_update_size,
        buffer_phase_size,
        buffer_phase_label,
        orientations=None,
    )
    return out


@func_mapper(task='generate_volume_element', method='random_voronoi_from_orientations')
def generate_volume_element_random_voronoi_orientations_2(microstructure_seeds, grid_size,
                                                          homog_label, scale_morphology,
                                                          scale_update_size,
                                                          buffer_phase_size,
                                                          buffer_phase_label,
                                                          orientations):
    out = generate_volume_element_random_voronoi(
        microstructure_seeds,
        grid_size,
        homog_label,
        scale_morphology,
        scale_update_size,
        buffer_phase_size,
        buffer_phase_label,
        orientations=orientations,
    )
    return out


def generate_volume_element_random_voronoi(microstructure_seeds, grid_size, homog_label,
                                           scale_morphology, scale_update_size,
                                           buffer_phase_size, buffer_phase_label,
                                           orientations=None):
    from damask import Grid

    grid_obj = Grid.from_Voronoi_tessellation(
        cells=np.array(grid_size),
        size=np.array(microstructure_seeds['size']),
        seeds=np.array(microstructure_seeds['position'])
    )

    if scale_morphology is not None:
        scale_morphology = np.array(scale_morphology)

        original_cells = grid_obj.cells
        new_cells = original_cells * scale_morphology
        grid_scaled = grid_obj.scale(new_cells)

        if scale_update_size:
            original_size = grid_obj.size
            new_size = original_size * scale_morphology
            grid_scaled.size = new_size

        grid_obj = grid_scaled

    phase_labels = [microstructure_seeds['phase_label']]
    if buffer_phase_size is not None:
        original_cells = grid_obj.cells
        original_size = grid_obj.size

        new_cells = original_cells + np.array(buffer_phase_size)
        new_size = original_size * (new_cells / original_cells)

        grid_canvased = grid_obj.canvas(cells=new_cells)
        grid_canvased.size = new_size

        grid_obj = grid_canvased
        phase_labels.append(buffer_phase_label)

    grid_obj.save_ASCII('geom.geom')

    volume_element = geom_to_volume_element(
        'geom.geom',
        phase_labels=phase_labels,
        homog_label=homog_label,
        orientations=(orientations or microstructure_seeds['orientations']),
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
