# Change Log

## [0.1.24] - 2022.02.15

### Changed

- Support `solver` and `initial_conditions` keys in the load YAML file.

## [0.1.23] - 2021.12.20

### Added

- Support passing `orientations` to task `generate_volume_element/from_damask_input_files`. Orientation in the passed `material.yml` file will be overwritten by those passed directly in `orientations`.
- Add task method: `generate_volume_element/single_voxel_grains` for generating an RVE that associates each grain with exactly one voxel.

## [0.1.22] - 2021.11.09

### Added

- Add saving of VE phases to a VTK file.

## [0.1.21] - 2021.09.24

### Changed

- Use `quat_component_ordering` key in orientations dict.

## [0.1.20] - 2021.08.14

### Changed

- Add parameter `RNG_seed` to task `generate_microstructure_seeds/random`.

## [0.1.19] - 2021.08.12

### Added

- Add task `generate_volume_element/dual_phase_ti_alpha_colony`.

## [0.1.18] - 2021.08.06

### Added

- Add implementation of Taylor simulation task.
- Add implementation of task `sample_texture` method `from_random`.
- Add implementation of task `generate_volume_element` method `from_damask_input_files`

### Changed

- Use new output map options (`volume_data`, `phase_data`, `field_data` and `grain_data`) for reading results from a simulation.
- Remove methods that made use of the DAMASK command line scripts.
- Add `orientations_use_max_precision` parameter to task/methods: `generate_volume_element/random_voronoi`, `generate_volume_element/random_voronoi_from_orientations`, `sample_texture/from_random`, `generate_microstructure_seeds/random` and `simulate_volume_element_loading/CP_FFT`.

## [0.1.17] - 2021.05.11

### Added

- Add visualise option to `simulate_volume_element_loading` output mapper to generate VTK file using DAMASK.

## [0.1.16] - 2021.04.10

### Added

- Add implementation of task `generate_microstructure_seeds` with method `random_NEW` to use DAMASK library functions instead of command line script which has been removed in DAMASK v3a2.

## [0.1.15] - 2020.01.19

### Added
- Add implementation of task `modify_volume_element` by method `add_buffer_zones`.

### Fixed 
 
- Include an `euler_degrees` bool dict item in `orientations` to signify if the dict item `euler_angles` is represented in degrees or radians.

## [0.1.14] - 2020.01.11

### Fixed

- Use `pack=False` in DAMASK's `geom_obj.to_file` to ensure consistent format of geometry file in function `main.generate_volume_element_random_voronoi`.

## [0.1.13] - 2020.01.10

### Changed

- Input and output mapper functions for writing/reading DAMASK microstructure seeds files have been updated to respect the `unit_cell_alignment` key of `orientations`.
- Changed some `generate_volume_element` task method names:
  - `random_voronoi` -> `random_voronoi_OLD`
  - `random_voronoi_from_orientations` -> `random_voronoi_from_orientations_OLD`
  - `random_voronoi_2` -> `random_voronoi`
  - `random_voronoi_from_orientations_2` -> `random_voronoi_from_orientations`

### Fixed

- Fix function name used twice for function mappers for `generate_volume_element` using method `random_voronoi_2` and `random_voronoi_from_orientations_2`.

## [0.1.12] - 2020.12.16

### Changed

- Add optional `single_crystal_parameter_perturbation` parameter to `simulate_volume_element_loading` task (used for parameter fitting).
- Change material/numerics.config to material/numerics.yaml.

## [0.1.11] - 2020.10.06

### Changed

- Standardise `orientations` key in `microstructure_seeds`.

### Added

- Add functions for `generate_volume_element` (`random_voronoi_2` and `random_voronoi_from_orientations_2`) using DAMASK Python package functions.

## [0.1.10] - 2020.09.29

### Changed

- Use `write_numerics` instead of `write_numerics_config`.

## [0.1.9] - 2020.09.29

### Changed

- Update for changes to damask-parse that support DAMASK v3.

## [0.1.8] - 2020.08.25

### Fixed

- Fix bug if numerics not specified in `write_damask_numerics` input mapper.

## [0.1.7] - 2020.08.25

### Added

- Add input mapper for writing the `numerics.config` file.

## [0.1.6] - 2020.08.22

### Changed

- Moved function `main.read_damask_hdf5_file` to `damask-parse` as `readers.read_HDF5_file`.
- Moved function `utils.get_HDF5_incremental_quantity` to `damask-parse` as `utils.get_HDF5_incremental_quantity` (and fixed use of correct Numpy function for "sum_along_axes" transform.)

## [0.1.5] - 2020.08.18

### Added

- Add support for specifying orientation and model coordinate systems with respect to e.g. RD/TD/ND.
- Add phase label attribute to `microstructure_seeds` and add `phase_labels` and `grain_phase_label_idx` to `volume_element`.

### Changed

- Function `write_damask_material` updated to reflect upstream changes.

## [0.1.4] - 2020.07.28

### Changed

- Add `image_axes` argument to function mapper function that generates an RVE via 2D microstructure extrusion: `volume_element_from_microstructure_image`.

## [0.1.3] - 2020.06.26

### Added

- Add function mapper for generating a volume element via extrusion of a 2D microstructure image.

## [0.1.2] - 2020.06.09

### Fixed

- Close file properly in `get_HDF5_incremental_quantity`

## [0.1.1] - 2020.05.12

### Added

- Added output mapper for DAMASK HDF5 files: `read_damask_hdf5_file`.

## [0.1.0] - 2020.05.09

- Initial release.
