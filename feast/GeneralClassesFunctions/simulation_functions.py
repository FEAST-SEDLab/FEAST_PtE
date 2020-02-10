"""
All functions that are required in the simulations and independent of input data and detection method are stored here.
"""

import numpy as np
import os
import pickle


def save_results(dir_out, results):
    """
    Save results to a file
    Inputs:
        dir_out             Name of output file to save
        results             A results object
    """
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)
    n_realization = len(os.listdir(dir_out))
    file_out = dir_out + '/realization' + str(n_realization) + '.p'
    pickle.dump(results, open(file_out, 'wb'))


def set_kwargs_attrs(obj_in, kwargs, only_existing=True):
    """
    Function for overwriting parameters with key word arguments
    Inputs:
        obj_in          Object with parameters to be updated
        kwargs          Dict containing new parameter values
        only_existing   Determines with new parameters can be created with kwargs
    """
    for key in kwargs.keys():
        # If only_existing is true, only set attributes that already exist
        if only_existing:
            if not hasattr(obj_in, key):
                raise ValueError("Tried to set invalid attribute. Class: ", type(obj_in), 'attempted attribute:', key)
        setattr(obj_in, key, kwargs[key])
