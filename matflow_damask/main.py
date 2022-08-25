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
    read_geom,
    read_material,
    read_HDF5_file,
    write_load_case,
    write_geom,
    write_material,
    write_numerics,
)
from damask_parse.utils import (
    get_header_lines,
    parse_damask_spectral_version_info,
    validate_orientations,
    volume_element_from_2D_microstructure,
    add_volume_element_buffer_zones,
    validate_orientations,
    validate_volume_element,
    spread_orientations,
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
from matflow_damask.utils import apply_single_crystal_parameter_perturbations, get_by_path, set_by_path


def read_orientation_coordinate_system(path):
    with Path(path).open('r') as handle:
        return json.load(handle)


@func_mapper(task='generate_microstructure_seeds', method='random')
def seeds_from_random(
    size,
    num_grains,
    phase_label=None,
    phase_labels=None,
    phase_fractions=None,
    grid_size=None,
    RNG_seed=None,
    orientation_coordinate_system=None,
    orientations_use_max_precision=False,
):
    from damask import seeds
    from damask import Rotation

    if phase_label is not None and phase_labels is not None:
        raise ValueError(f"Specify exactly one of `phase_label` and `phase_labels`.")

    if (
        (phase_label is not None and not isinstance(phase_label, str)) or 
        (phase_labels is not None and not isinstance(phase_labels, list))
    ):
        raise ValueError(
            f"Specify `phase_label` as a string, or `phase_labels` as a list of strings."
        )

    if phase_labels is None:
        phase_labels = [phase_label]
        phase_labels_idx = np.zeros(num_grains)

    num_phase_labels = len(phase_labels)
    if phase_fractions is None:        
        phase_fractions = [1/num_phase_labels for _ in phase_labels]
    else:
        if len(phase_fractions) != len(phase_labels):
            raise ValueError(
                f"Length of `phase_fractions` ({len(phase_fractions)}) must be equal to "
                f"length of `phase_labels` ({len(phase_labels)})."
            )
        if sum(phase_fractions) != 1.0:
            raise ValueError(
                f"Sum of `phase_fractions` ({sum(phase_fractions)}) must sum to one."
            )
    
    # Assign phases to labels, randomly:
    rng = np.random.default_rng(seed=RNG_seed)
    phase_labels_idx = rng.choice(a=num_phase_labels, size=num_grains, p=phase_fractions)

    size = np.array(size)
    grid_size = grid_size and np.array(grid_size)

    position = seeds.from_random(size, num_grains, cells=grid_size, rng_seed=RNG_seed)
    rotation = Rotation.from_random(shape=(num_grains,))

    oris = {
        'type': 'quat',
        'quaternions': rotation.quaternion,
        'quat_component_ordering': 'scalar-vector',
        'orientation_coordinate_system': orientation_coordinate_system,
        'unit_cell_alignment': {
            'x': 'a',
            'z': 'c',
        },
        'P': -1,
        'use_max_precision': orientations_use_max_precision,
    }
    oris = validate_orientations(oris)

    out = {
        'microstructure_seeds': {
            'position': position,
            'orientations': oris,
            'size': size,
            'random_seed': None,
            'phase_labels': phase_labels,
            'phase_labels_idx': phase_labels_idx,
            'phase_fractions': phase_fractions,
        }
    }
    return out


@func_mapper(task='sample_texture', method='from_random')
def orientations_from_random(num_orientations,
                             orientation_coordinate_system=None,
                             orientations_use_max_precision=False):
    from damask import Rotation

    rotation = Rotation.from_random(shape=(num_orientations,))

    oris = {
        'type': 'quat',
        'quaternions': rotation.quaternion,
        'quat_component_ordering': 'scalar-vector',
        'orientation_coordinate_system': orientation_coordinate_system,
        'unit_cell_alignment': {
            'x': 'a',
            'z': 'c',
        },
        'P': -1,
        'use_max_precision': orientations_use_max_precision,
    }
    out = {'orientations': validate_orientations(oris)}
    return out


@input_mapper('load.yaml', 'simulate_volume_element_loading', 'CP_FFT')
@input_mapper('load.yaml', 'simulate_orientations_loading', 'Taylor')
def write_damask_load_case(path, load_case, solver, initial_conditions):
    path = Path(path)
    write_load_case(path.parent, load_case, solver, initial_conditions, name=path.name)


@input_mapper('geom.vtr', 'simulate_volume_element_loading', 'CP_FFT')
def write_damask_geom(path, volume_element):
    path = Path(path)
    write_geom(path.parent, volume_element, name=path.name)


@input_mapper('material.yaml', 'simulate_volume_element_loading', 'CP_FFT')
def write_damask_material(path, homogenization_schemes, volume_element,
                          single_crystal_parameters,
                          single_crystal_parameter_perturbation, phases,
                          texture_alignment_method,
                          orientations_use_max_precision):

    path = Path(path)

    # Apply perturbations to single-crystal parameters:
    if single_crystal_parameters and single_crystal_parameter_perturbation:
        single_crystal_parameters = apply_single_crystal_parameter_perturbations(
            single_crystal_parameters,
            single_crystal_parameter_perturbation,
        )
        # HACK: write out single crystal parameters as a JSON file, so they can be parsed
        # as an output of the task (e.g. as "perturbed_single_crystal_parameters"), if we
        # want:
        with path.parent.joinpath('single_crystal_parameters.json').open('wt') as fh:
            json.dump(single_crystal_parameters, fh)

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

    if orientations_use_max_precision is not None:
        volume_element['orientations'].update({
            'use_max_precision': orientations_use_max_precision
        })
    write_material(
        homog_schemes=homogenization_schemes,
        phases=phases,
        volume_element=volume_element,
        dir_path=path.parent,
        name=path.name,
    )


@input_mapper('material.yaml', 'simulate_orientations_loading', 'Taylor')
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
    phase_name = list(phases.keys())[0]

    volume_element = {
        'element_material_idx': np.arange(8).reshape([2, 2, 2]),
        'grid_size': np.array([2, 2, 2]),
        'size': [1.0, 1.0, 1.0],

        'constituent_material_idx': np.arange(8).repeat(num_oris_per_mat),
        'constituent_material_fraction': np.full(num_oris, 1. / num_oris_per_mat),
        'constituent_phase_label': np.full(num_oris, phase_name),
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
@input_mapper('numerics.yaml', 'simulate_orientations_loading', 'Taylor')
def write_damask_numerics(path, numerics):
    if numerics:
        path = Path(path)
        write_numerics(path.parent, numerics, name=path.name)


@output_mapper('perturbed_single_crystal_parameters', 'simulate_volume_element_loading', 'CP_FFT')
def read_single_crystal_parameters_JSON(json_path):
    json_path = Path(json_path)
    if json_path.is_file():
        with Path(json_path).open("rt") as fh:
            return json.load(fh)
    else:
        return None

@output_mapper('volume_element_response', 'simulate_volume_element_loading', 'CP_FFT')
def read_damask_hdf5_file(hdf5_path, incremental_data=None, volume_data=None,
                          phase_data=None, field_data=None, grain_data=None,
                          operations=None, visualise=None):

    out = read_HDF5_file(hdf5_path, incremental_data=incremental_data,
                         volume_data=volume_data, phase_data=phase_data,
                         field_data=field_data, grain_data=grain_data,
                         operations=operations)

    if visualise is not None:

        from damask import Result
        from damask_parse.utils import parse_inc_specs_using_result_obj

        if visualise is True:
            visualise = [{}]
        elif isinstance(visualise, dict):
            visualise = [visualise]

        os.mkdir('viz')
        with working_directory('viz'):

            result = Result(hdf5_path)

            for viz_dict_idx, viz_dict in enumerate(visualise, 1):

                if len(visualise) > 1:
                    viz_dir = str(viz_dict_idx)
                    os.mkdir(viz_dir)
                else:
                    viz_dir = '.'
                with working_directory(viz_dir):

                    try:
                        # all incs if not specified:
                        incs_spec = viz_dict.get('increments', None)
                        parsed_incs = parse_inc_specs_using_result_obj(incs_spec, result)
                        result = result.view('increments', parsed_incs)

                        # all phases if not specified:
                        phases = viz_dict.get('phases', True)
                        result = result.view('phases', phases)

                        # all homogs if not specified:
                        homogs = viz_dict.get('homogenizations', True)
                        result = result.view('homogenizations', homogs)

                        # all outputs if not specified:
                        outputs = viz_dict.get('fields', '*')

                        outputs = visualise_static_outpurts(outputs, result, out)

                        result.save_VTK(output=outputs)

                    except Exception as err:
                        print(f'Could not save VTK files for visualise item: {viz_dict}. '
                              f'Exception was: {err}')
                        continue

    return out


def visualise_static_outpurts(outputs, result, out):
    """Create separate VTK file for grain and phase maps."""

    static_outputs = ['grain', 'phase']

    if isinstance(outputs, list):
        static_outputs = list(set(outputs).intersection(static_outputs))
        if len(static_outputs) > 0:
            v = result.geometry0

            for static_output in static_outputs:
                outputs.remove(static_output)
                dat_array = out['field_data'][static_output]['data']
                v.add(dat_array.flatten(order='F'), label=static_output)

            v.save('static_outputs')

    return outputs


@output_mapper('orientations_response', 'simulate_orientations_loading', 'Taylor')
def read_damask_hdf5_file_2(hdf5_path, incremental_data=None, volume_data=None,
                            phase_data=None, operations=None):
    return read_HDF5_file(hdf5_path, incremental_data=incremental_data,
                          volume_data=volume_data, phase_data=phase_data,
                          operations=operations)


@func_mapper(task='generate_volume_element', method='extrusion')
def volume_element_from_microstructure_image(microstructure_image, depth,
                                             image_axes, phase_label,
                                             phase_label_mapping, homog_label):
    out = {
        'volume_element': volume_element_from_2D_microstructure(
            microstructure_image, homog_label, phase_label,
            phase_label_mapping, depth, image_axes,
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


@func_mapper(task='modify_volume_element', method='spread_orientations')
def modify_volume_element_spread_orientations(volume_element, phases, stddev_degrees):
    VE = spread_orientations(volume_element, phases, sigmas=stddev_degrees)
    return {"volume_element": VE}


@func_mapper(task='modify_volume_element', method='new_orientations')
def modify_volume_element_new_orientations(volume_element, volume_element_response):

    n_grains = volume_element['orientations']['quaternions'].shape[0]
    n_fragments = volume_element_response['incremental_data']['orientations']['data']['quaternions'].shape[1]

    old_oris = volume_element_response['incremental_data']['orientations']['data']['quaternions'][-1]
    random_index = np.random.randint(n_fragments, size=n_grains)
    print("randomindex array size: ", random_index.shape)
    volume_element['orientations']['quaternions'] = old_oris[random_index, :]

    print("old orientations array shape: ", old_oris.shape,
          "\nnew orientations array shape: ", old_oris[random_index, :].shape)  # DEBUG

    out = {'volume_element': volume_element}
    return out


@func_mapper(task='modify_volume_element', method='geometry')
def modify_volume_element_geometry(volume_element, volume_element_response):

    print("\nold_geomsize:", volume_element['size'])  # DEBUG
    volume_element['size'] = np.matmul(
        volume_element_response['incremental_data']['def_grad']['data'][-1], volume_element['size'])
    print("\nnew_geomsize:", volume_element['size'])  # DEBUG

    out = {'volume_element': volume_element}
    return out


@func_mapper(task='modify_volume_element', method='grid_size')
def modify_volume_element_geometry(volume_element, new_grid_size):    
    
    from scipy.ndimage import zoom

    new_grid_size = np.array(new_grid_size)
    zoom_factor = np.array(new_grid_size) / volume_element['grid_size']
    new_elem_mat_idx = zoom(volume_element['element_material_idx'], zoom_factor, order=0)
    
    volume_element['grid_size'] = new_grid_size
    volume_element['element_material_idx'] = new_elem_mat_idx

    volume_element = validate_volume_element(volume_element)    
    out = {'volume_element': volume_element}
    return out    


@func_mapper(task='visualise_volume_element', method='VTK')
def visualise_volume_element(volume_element):
    path = Path('geom.vtr')
    write_geom(path.parent, volume_element, name=path.name)


@func_mapper(task='generate_volume_element', method='random_voronoi')
def generate_volume_element_random_voronoi_2(microstructure_seeds, grid_size, homog_label,
                                             scale_morphology, scale_update_size,
                                             buffer_phase_size, buffer_phase_label,
                                             orientations_use_max_precision, periodic):
    out = generate_volume_element_random_voronoi(
        microstructure_seeds,
        grid_size,
        homog_label,
        scale_morphology,
        scale_update_size,
        buffer_phase_size,
        buffer_phase_label,
        orientations_use_max_precision,
        orientations=None,
        periodic=periodic,
    )
    return out


@func_mapper(task='generate_volume_element', method='random_voronoi_from_orientations')
def generate_volume_element_random_voronoi_orientations_2(microstructure_seeds, grid_size,
                                                          homog_label, scale_morphology,
                                                          scale_update_size,
                                                          buffer_phase_size,
                                                          buffer_phase_label,
                                                          orientations,
                                                          orientations_use_max_precision, 
                                                          periodic):
    out = generate_volume_element_random_voronoi(
        microstructure_seeds,
        grid_size,
        homog_label,
        scale_morphology,
        scale_update_size,
        buffer_phase_size,
        buffer_phase_label,
        orientations_use_max_precision,
        orientations=orientations,
        periodic=periodic,
    )
    return out


def generate_volume_element_random_voronoi(
    microstructure_seeds,
    grid_size,
    homog_label,
    scale_morphology,
    scale_update_size,
    buffer_phase_size,
    buffer_phase_label,
    orientations_use_max_precision,
    orientations=None,
    orientations_idx=None,
    periodic=True,
):
    from damask import Grid

    grid_obj = Grid.from_Voronoi_tessellation(
        cells=np.array(grid_size),
        size=np.array(microstructure_seeds['size']),
        seeds=np.array(microstructure_seeds['position']),
        periodic=periodic,
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

    const_phase_lab = np.array(microstructure_seeds['phase_labels'])[
        np.array(microstructure_seeds['phase_labels_idx'])
    ]

    num_grains = len(microstructure_seeds['position'])
    if orientations is not None:
        oris = orientations
    else:
        oris = copy.deepcopy(microstructure_seeds['orientations'])
    oris.update({'use_max_precision': orientations_use_max_precision})

    if orientations_idx is not None:
        ori_idx = orientations_idx
    else:
        ori_idx = np.arange(num_grains)

    if buffer_phase_size is not None:
        original_cells = grid_obj.cells
        original_size = grid_obj.size

        new_cells = original_cells + np.array(buffer_phase_size)
        new_size = original_size * (new_cells / original_cells)

        grid_canvased = grid_obj.canvas(cells=new_cells)
        grid_canvased.size = new_size

        grid_obj = grid_canvased

        if buffer_phase_label is None:
            raise ValueError(
                f"Must specify `buffer_phase_label` if specifying `buffer_phase_size`."
            )
        const_phase_lab = np.append(const_phase_lab, [buffer_phase_label])
        oris['quaternions'] = np.vstack([oris['quaternions'], np.array([[1, 0, 0, 0]])])
        ori_idx = np.append(ori_idx, [oris['quaternions'].shape[0] - 1])
        num_grains += 1

    volume_element = {
        'size': grid_obj.size.astype(float).tolist(),
        'grid_size': grid_obj.cells.tolist(),
        'orientations': oris,
        'element_material_idx': grid_obj.material,
        'constituent_material_idx': np.arange(num_grains),
        'constituent_material_fraction': np.ones(num_grains),
        'constituent_phase_label': const_phase_lab,
        'constituent_orientation_idx': ori_idx,
        'material_homog': np.full(num_grains, homog_label),
    }
    volume_element = validate_volume_element(volume_element)
    return {'volume_element': volume_element}

@func_mapper(
    task='generate_volume_element',
    method='random_voronoi_from_dual_phase_orientations'
)
def generate_volume_element_from_random_voronoi_dual_phase_orientations(
    microstructure_seeds,
    grid_size,
    homog_label,
    scale_morphology,
    scale_update_size,
    buffer_phase_size,
    buffer_phase_label,
    orientations_use_max_precision,
    orientations_phase_1,
    orientations_phase_2,
    RNG_seed,
    periodic,
):

    if len(microstructure_seeds['phase_labels']) != 2:
        raise ValueError(
            f"There are not two phase labels to correspond to two orientation sets. "
            f"`phase_labels` is: {microstructure_seeds['phase_labels']}.")

    ori_1 = validate_orientations(orientations_phase_1)
    ori_2 = validate_orientations(orientations_phase_2)
    oris = copy.deepcopy(ori_1)
    
    phase_labels_idx = microstructure_seeds['phase_labels_idx']
    phase_labels = microstructure_seeds['phase_labels']
    num_grains = len(phase_labels_idx)
    _, counts = np.unique(phase_labels_idx, return_counts=True)
    
    num_ori_1 = ori_1['quaternions'].shape[0]
    num_ori_2 = ori_2['quaternions'].shape[0]
    sampled_oris_1 = ori_1['quaternions']
    sampled_oris_2 = ori_2['quaternions']

    rng = np.random.default_rng(seed=RNG_seed)

    # If there are more orientations than phase label assignments, choose a random subset:
    if num_ori_1 != counts[0]:
        try:
            ori_1_idx = rng.choice(a=num_ori_1, size=counts[0], replace=False)
        except ValueError as err:
            raise ValueError(
                f"Probably an insufficient number of `orientations_phase_1` "
                f"({num_ori_1} given for phase {phase_labels[0]!r}, whereas {counts[0]} "
                f"needed). Caught ValueError is: {err}"
            )
        sampled_oris_1 = sampled_oris_1[ori_1_idx]
    if num_ori_2 != counts[1]:
        try:
            ori_2_idx = rng.choice(a=num_ori_2, size=counts[1], replace=False)
        except ValueError as err:
            raise ValueError(
                f"Probably an insufficient number of `orientations_phase_2` "
                f"({num_ori_2} given for phase {phase_labels[1]!r}, whereas {counts[1]} "
                f"needed). Caught ValueError is: {err}"
            )
        sampled_oris_2 = sampled_oris_2[ori_2_idx]

    ori_idx = np.ones(num_grains) * np.nan
    for idx, i in enumerate(counts):
        ori_idx[phase_labels_idx == idx] = np.arange(i) + np.sum(counts[:idx])

    if np.any(np.isnan(ori_idx)):
        raise RuntimeError("Not all phases have an orientation assigned!")
    ori_idx = ori_idx.astype(int)
    
    oris['quaternions'] = np.vstack([sampled_oris_1, sampled_oris_2])
    
    outputs = generate_volume_element_random_voronoi(
        microstructure_seeds,
        grid_size,
        homog_label,
        scale_morphology,
        scale_update_size,
        buffer_phase_size,
        buffer_phase_label,
        orientations_use_max_precision,
        orientations=oris,
        orientations_idx=ori_idx,
        periodic=periodic,
    )
    return outputs


@func_mapper(task='generate_volume_element', method='from_damask_input_files')
def generate_volume_element_from_damask_input_files(geom_path, material_path, orientations):
    geom_dat = read_geom(geom_path)
    material_data = read_material(material_path)
    volume_element = {
        'element_material_idx': geom_dat['element_material_idx'],
        'grid_size': geom_dat['grid_size'],
        'size': geom_dat['size'],
        **material_data['volume_element'],
    }

    if orientations is not None:
        orientations = validate_orientations(orientations)
        num_supplied_ori = orientations['quaternions'].shape[0]
        num_material_ori = volume_element['orientations']['quaternions'].shape[0]
        if num_supplied_ori != num_material_ori:
            raise ValueError(
                f'Number of orientations supplied {num_supplied_ori} is different to '
                f'number in the material file {num_material_ori}.'
            )

        volume_element['orientations'] = orientations

    volume_element = validate_volume_element(volume_element)
    out = {'volume_element': volume_element}
    return out


@func_mapper(task='generate_volume_element', method='dual_phase_ti_alpha_colony')
def generate_RVE_dual_phase_ti_alpha_colony(grid_size, alpha_particle_axes_ratio,
                                            alpha_particle_centres, alpha_orientation,
                                            beta_orientation):

    from damask import seeds, Grid

    my_seeds = seeds.from_random([1, 1, 1], 1)
    my_grid = Grid.from_Voronoi_tessellation(
        np.array(grid_size),
        [1, 1, 1],
        my_seeds,
        periodic=True
    )

    # Ellipsoid axes ratio: (2 : 1 : 0.4)
    plate_dims = np.array([0.35] * 3)

    plate_dims_norm_factor = 2 / np.max(alpha_particle_axes_ratio)
    plate_dims[0] *= alpha_particle_axes_ratio[0] * plate_dims_norm_factor
    plate_dims[1] *= alpha_particle_axes_ratio[1] * plate_dims_norm_factor
    plate_dims[2] *= alpha_particle_axes_ratio[2] * plate_dims_norm_factor

    # Fixed for now until/if we investigate effect of proximity of the alpha laths:
    centres = np.array([
        [0.3, 0.10, 0.3],
        [0.45, 0.65, 0.6],
        [0.85, 0.30, 0.9],
    ])

    for centre in centres:
        my_grid = my_grid.add_primitive(
            dimension=plate_dims,
            center=centre,
            exponent=1,
        )

    oris = {
        'type': 'quat',
        'quaternions': np.array([
            beta_orientation,
            alpha_orientation,
        ]),
        'quat_component_ordering': 'scalar-vector',
        'unit_cell_alignment': {'x': 'a'},
        'P': 1,
    }

    volume_element = {
        'constituent_material_idx': np.arange(1 + len(centres)),
        'constituent_phase_label': ['Ti-beta', 'Ti-alpha', 'Ti-alpha', 'Ti-alpha'],
        'constituent_orientation_idx': [0] + [1] * len(centres),
        'material_homog': ['SX'] * (1 + len(centres)),
        'element_material_idx': my_grid.material,
        'grid_size': my_grid.cells,
        'orientations': oris,
    }
    volume_element = validate_volume_element(volume_element)

    return {'volume_element': volume_element}


@func_mapper(task='generate_volume_element', method='single_voxel_grains')
def generate_volume_element_single_voxel_grains(grid_size, size, homog_label,
                                                scale_morphology, scale_update_size,
                                                phase_label, buffer_phase_size,
                                                buffer_phase_label, orientations,
                                                orientations_use_max_precision):

    from damask import Grid

    grid_size = np.array(grid_size)
    size = np.array(size)

    grid_obj = Grid(
        material=np.arange(np.product(grid_size)).reshape(grid_size, order='F'),
        size=size,
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

    phase_labels = [phase_label]
    if buffer_phase_size is not None:
        original_cells = grid_obj.cells
        original_size = grid_obj.size

        new_cells = original_cells + np.array(buffer_phase_size)
        new_size = original_size * (new_cells / original_cells)

        grid_canvased = grid_obj.canvas(cells=new_cells)
        grid_canvased.size = new_size

        grid_obj = grid_canvased
        phase_labels.append(buffer_phase_label)

    orientations.update({'use_max_precision': orientations_use_max_precision})
    volume_element = {
        'orientations': orientations,
        'element_material_idx': grid_obj.material,
        'grid_size': grid_obj.cells.tolist(),
        'size': grid_obj.size.astype(float).tolist(),
        'phase_labels': phase_labels,
        'homog_label': homog_label,
    }
    volume_element = validate_volume_element(volume_element)
    return {'volume_element': volume_element}

@func_mapper(task='optimise_single_crystal_parameters', method='bayesian')
def optimise_single_crystal_parameters_bayesian(
    perturbed_volume_element_responses,
    perturbed_single_crystal_parameters,
    single_crystal_parameter_perturbations,
    experimental_tensile_test,
    prior_distribution,
):    
    # select new parameters at random for now:
    single_crystal_parameters_new = perturbed_single_crystal_parameters[0]
    posterior_distribution = None

    return {
        'single_crystal_parameters': single_crystal_parameters_new,
        'posterior_distribution': posterior_distribution,
    }

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
