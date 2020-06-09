import h5py
import numpy as np


def get_HDF5_incremental_quantity(hdf5_path, dat_path, transforms=None):
    'Accessing HDF5 file directly.'

    with h5py.File(str(hdf5_path), 'r') as f:

        incs = [i for i in f.keys() if 'inc' in i]
        data = np.array([f[i][dat_path][()] for i in incs])

        if transforms:
            for i in transforms:
                if 'mean_along_axes' in i:
                    data = np.mean(data, i['mean_along_axes'])
                if 'sum_along_axes' in i:
                    data = np.mean(data, i['sum_along_axes'])

        return data
