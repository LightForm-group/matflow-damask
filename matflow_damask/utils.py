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
