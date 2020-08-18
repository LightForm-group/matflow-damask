# Change Log

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
