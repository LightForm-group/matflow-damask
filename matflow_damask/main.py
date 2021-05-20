'`matflow_damask.main.py`'

import os
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
    validate_orientations,
    volume_element_from_2D_microstructure,
    add_volume_element_buffer_zones,
)
from damask_parse import __version__ as damask_parse_version
from matflow.scripting import get_wrapper_script
from matflow.utils import working_directory
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


@input_mapper(
    input_file='orientation_coordinate_system.json',
    task='generate_volume_element',
    method='random_voronoi_OLD'
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
    method='random_voronoi_from_orientations_OLD'
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

    size = None
    random_seed = None
    for ln in header_lns:
        if ln.startswith('size'):
            size = [float(j) for j in [i for i in ln.split()][1:][1::2]]
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
            'euler_degrees': True,
            'orientation_coordinate_system': orientation_coordinate_system,
            # DAMASK uses x//a alignment for hexagonal system (unlike the more
            # common y//b):
            'unit_cell_alignment': {
                'x': 'a',
                'z': 'c',
            }
        },
        'size': size,
        'random_seed': random_seed,
        'phase_label': phase_label,
    }

    return microstructure_seeds


@func_mapper(task='generate_microstructure_seeds', method='random_NEW')
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


@output_mapper('volume_element', 'generate_volume_element', 'random_voronoi_OLD')
@output_mapper('volume_element', 'generate_volume_element', 'random_voronoi_from_orientations_OLD')
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

        SC_params_name = phases[phase_label]['plasticity'].pop(
            'single_crystal_parameters',
            None
        )
        if SC_params_name:
            SC_params = single_crystal_parameters[SC_params_name]
            phases[phase_label]['plasticity'].update(**SC_params)
    write_material(
        homog_schemes=homogenization_schemes,
        phases=phases,
        volume_element=volume_element,
        dir_path=Path(path).parent,
    )


@input_mapper('numerics.yaml', 'simulate_volume_element_loading', 'CP_FFT')
def write_damask_numerics(path, numerics):
    if numerics:
        write_numerics(Path(path).parent, numerics)


@output_mapper('volume_element_response', 'simulate_volume_element_loading', 'CP_FFT')
def read_damask_hdf5_file(hdf5_path, incremental_data, operations=None, visualise=None):

    out = read_HDF5_file(hdf5_path, incremental_data, operations=operations)

    if visualise is not None:

        if visualise is True:
            visualise = {}

        os.mkdir('viz')
        with working_directory('viz'):

            from damask import Result

            result = Result(hdf5_path)

            incs = visualise.pop('increments', None)
            if incs:
                if not isinstance(incs, list):
                    incs = [incs]
                incs_normed = []
                for i in incs:
                    if i >= 0:
                        i_normed = i
                    else:
                        i_normed = len(result.increments) + i
                    incs_normed.append(i_normed)
                result.pick('increments', incs_normed)
            result.to_vtk(**visualise)

    return out


@input_mapper('phase_label.txt', 'generate_volume_element', 'random_voronoi_OLD')
@input_mapper('phase_label.txt', 'generate_volume_element', 'random_voronoi_from_orientations_OLD')
def write_phase_label(path, microstructure_seeds):
    with Path(path).open('w') as handle:
        handle.write(f'{microstructure_seeds["phase_label"]}\n')


@input_mapper('homog_label.txt', 'generate_volume_element', 'random_voronoi_OLD')
@input_mapper('homog_label.txt', 'generate_volume_element', 'random_voronoi_from_orientations_OLD')
def write_homog_label(path, homog_label):
    with Path(path).open('w') as handle:
        handle.write(f'{homog_label}\n')


@input_mapper('orientation.seeds', 'generate_volume_element', 'random_voronoi_OLD')
def write_microstructure_seeds(path, microstructure_seeds):

    grid_size = microstructure_seeds['grid_size']
    position = microstructure_seeds['position']

    eulers = microstructure_seeds['orientations']['euler_angles']
    if not microstructure_seeds['orientations']['euler_degrees']:
        eulers = np.rad2deg(eulers)

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
    'random_voronoi_from_orientations_OLD'
)
def write_microstructure_new_orientations(path, microstructure_seeds, orientations):

    grid_size = microstructure_seeds['grid_size']
    position = microstructure_seeds['position']

    if (
        'unit_cell_alignment' not in orientations or
        'x' not in orientations['unit_cell_alignment'] or
        orientations['unit_cell_alignment']['x'] != 'a'
    ):
        msg = ('Orientations to be written to a DAMASK seeds file must have '
               'DAMASK-compatible `unit_cell_alignment`: x parallel to a.')
        # TODO: allow conversion here.
        raise NotImplementedError(msg)

    eulers = orientations['euler_angles']
    if not orientations['euler_degrees']:
        eulers = np.rad2deg(eulers)

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


@cli_format_mapper('size', 'generate_volume_element', 'random_voronoi_from_orientations_OLD')
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

@func_mapper(task='modify_volume_element', method='add_buffer_zones')
def modify_volume_element_add_buffer_zones(volume_element, buffer_sizes,
                                           phase_ids, phase_labels, homog_label, order):
    out = {
        'volume_element': add_volume_element_buffer_zones(
            volume_element, buffer_sizes, phase_ids, phase_labels, homog_label, order
        )
    }
    return out

@func_mapper(task='modify_volume_element', method='new_orientations')
def modify_volume_element_new_orientations(volume_element, volume_element_response):

    n_grains = volume_element['orientations']['quaternions'].shape[0]
    n_fragments = volume_element_response['orientations']['data']['quaternions'][-1].shape[0]

    old_oris = volume_element_response['orientations']['data']['quaternions'][-1]
    random_index = np.random.randint(n_fragments, size=n_grains) ; print("randomindex array size: ", random_index.shape)
    volume_element['orientations']['quaternions'] = old_oris[random_index, :]

    print("old orientations array shape: ", old_oris.shape, "\nnew orientations array shape: ", old_oris[random_index, :].shape) # DEBUG

    # return volume_element with new oris...
    out = {'volume_element': volume_element} ; print("volume_element['orientations']['quaternions']:\n", volume_element['orientations']['quaternions'])
    return out

@func_mapper(task='modify_volume_element', method='geometry')
def modify_volume_element_geometry(volume_element, volume_element_response):

    print("\nold_geomsize:", volume_element['size']) # DEBUG
    volume_element['size'] = np.matmul(volume_element_response['def_grad']['data'][-1], volume_element['size'])
    print("\nnew_geomsize:", volume_element['size']) # DEBUG
    
    # return volume_element with new size ...
    out = { 'volume_element': volume_element }
    return out
    

@func_mapper(task='visualise_volume_element', method='VTK')
def visualise_volume_element(volume_element):
    try:
        from damask import Geom

        geom_obj = Geom.from_file('geom.geom')
        geom_obj.to_vtr('geom.vtr')
    except ImportError:
        from damask import Grid

        grid_obj = Grid.load_ASCII('geom.geom')
        grid_obj.save('geom.vtr')


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
    try:
        from damask import Geom

        geom_obj = Geom.from_Voronoi_tessellation(
            grid=np.array(grid_size),
            size=np.array(microstructure_seeds['size']),
            seeds=np.array(microstructure_seeds['position']),
        )

        if scale_morphology is not None:
            scale_morphology = np.array(scale_morphology)

            original_grid = geom_obj.get_grid()
            new_grid = original_grid * scale_morphology
            geom_scaled = geom_obj.scale(new_grid)

            if scale_update_size:
                original_size = geom_obj.get_size()
                new_size = original_size * scale_morphology
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

        # specifying pack ensures consistent behaviour:
        geom_obj.to_file('geom.geom', pack=False)

    except ImportError:
        from damask import Grid as Geom

        grid_obj = Geom.from_Voronoi_tessellation(
            cells=np.array(grid_size),
            size=np.array(microstructure_seeds['size']),
            seeds=np.array(microstructure_seeds['position']),
            material=np.arange(1, microstructure_seeds['position'].shape[0]+1)
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
