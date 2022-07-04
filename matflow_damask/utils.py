import copy

import numpy as np


def get_by_path(root, path):
    """Get a nested dict or list item according to its "key path"
    
    Parameters
    ----------
    root : dict or list
        Can be arbitrarily nested.
    path : list of str
        The address of the item to get within the `root` structure.
        
    Returns
    -------
    sub_data : any
        
    """
    
    sub_data = root
    for key in path:
        sub_data = sub_data[key]
    
    return sub_data

def set_by_path(root, path, value):
    """Set a nested dict or list item according to its "key path"
    
    Parmaeters
    ----------
    root : dict or list
        Can be arbitrarily nested.
    path : list of str
        The address of the item to set within the `root` structure.
    value : any
        The value to set.
        
    """
    
    sub_data = root
    for key in path[:-1]:
        sub_data = sub_data[key]
    sub_data[path[-1]] = value


def apply_single_crystal_parameter_perturbations(single_crystal_params, pert):
    
    if 'address' in pert and 'addresses' in pert:
        raise ValueError(f"Specify either `address` or `addresses`.")
    if 'perturbation' in pert and 'perturbations' in pert:
        raise ValueError(f"Specify either `perturbation` or `perturbations`.")
    if 'std_deviation' in pert and 'std_deviations' in pert:
        raise ValueError(f"Specify either `std_deviation` or `std_deviations`.")

    is_scale = 'perturbation' in pert or 'perturbations' in pert
    is_std_dev = 'std_deviation' in pert or 'std_deviations' in pert
    if (is_scale and is_std_dev) or (not is_scale and not is_std_dev):
        raise ValueError("Specify perturbation as either `perturbations` or `std_deviations`")

    if is_scale:
        if ('perturbation' in pert and 'addresses' in pert) or ('perturbations' in pert and 'address' in pert):
            raise ValueError(f"Specify `perturbation` and `address` or `perturbations` and `addresses`.")
        if 'perturbations' in pert and len(pert['perturbations']) != len(pert['addresses']):
            raise ValueError(f"The lengths of `perturbations` and `addresses` must be equal")

    elif is_std_dev:
        if ('std_deviation' in pert and 'addresses' in pert) or ('std_deviations' in pert and 'address' in pert):
            raise ValueError(f"Specify `std_deviation` and `address` or `std_deviations` and `addresses`.")
        if 'std_deviations' in pert and len(pert['std_deviations']) != len(pert['addresses']):
            raise ValueError(f"The lengths of `perturbations` and `addresses` must be equal")    

    if 'perturbation' in pert:
        perturbations = [pert['perturbation']]
        addresses = [pert.get('address')]
    elif 'std_deviation' in pert:
        std_devs = [pert['std_deviations']]
        addresses = [pert.get('address')]
    else:
        perturbations = pert.get('perturbations')
        std_devs = pert.get('std_deviations')
        addresses = pert.get('addresses')
    
    params = copy.deepcopy(single_crystal_params)
    if is_scale:
        for idx, (pert_i, add_i) in enumerate(zip(perturbations, addresses)):
            if add_i is not None:
                scale_factor_i = (1 + pert_i)
                new_val = get_by_path(params, add_i) * scale_factor_i
                for eq_con in pert.get('equal_constraints', []):
                    eq_con = sorted(eq_con)
                    if idx > eq_con[0] and idx in eq_con[1:]:
                        new_val = get_by_path(params, addresses[eq_con[0]])
                set_by_path(params, add_i, new_val)
    elif is_std_dev:
        for idx, (std_dev_i, add_i) in enumerate(zip(std_devs, addresses)):
            if add_i is not None:
                # draw a perturbation from a zero-mean normal distribution
                rng = np.random.default_rng()
                pert_i = rng.normal(loc=0.0, scale=std_dev_i)
                scale_factor_i = (1 + pert_i)
                new_val = get_by_path(params, add_i) * scale_factor_i
                for eq_con in pert.get('equal_constraints', []):
                    eq_con = sorted(eq_con)
                    if idx > eq_con[0] and idx in eq_con[1:]:
                        new_val = get_by_path(params, addresses[eq_con[0]])
                set_by_path(params, add_i, new_val)


    return params
