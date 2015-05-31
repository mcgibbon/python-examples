# -*- coding: utf-8 -*-
"""
load_netcdf_full.py: an example on loading data from a netCDF file. Some of
the operations described in this file may be rarely used. If you would like
a shorter example, see load_netcdf_simple.py

This code is released into the public domain. For full details, see
UNLICENSE at the root of this project. In essence, you can copy and use this
code for whatever purpose, with or without attribution.

If you have any questions about the code or could use any sort of help using
it, feel free to e-mail me at mcgibbon (at) uw {dot} edu.
"""
# we will be using the netCDF4 library
import netCDF4 as nc4

# This 'if' statement should always be put at the start of your
# 'script' section after any function or class definitions, or global
# variables.
if __name__ == '__main__':
    # netCDF4 Dataset objects are by default opened in read-only mode
    # put the filename you would like to load as the first argument
    dataset = nc4.Dataset('../../sample_data/ARM_sounding_example.nc')
    # data is stored in the dataset.variables dictionary
    # its keys are variable names, and values are Variable objects
    p = dataset.variables['pres']
    T = dataset.variables['tdry']
    # attributes of netcdf variables are stored as attributes of Variable
    # objects
    print(p.long_name)
    print(p.units)
    print(T.long_name)
    print(T.units)
    # the same is true of Dataset objects
    print(dataset.history)
    # Math operations do not work on variables, and need numpy
    # arrays instead. To get the array only, use [:] or any kind of slicing
    # operation. This can be done within your mathematical operation.
    try:
        p*p
    except TypeError:
        print('Raised TypeError for Variable object.')
    try:
        p[:]*p[:]
    except TypeError:
        print('Raised TypeError for numpy array.')
    # Note that since these are now numpy arrays, we *cannot* access
    # attributes of those variables
    print(hasattr(p[:], 'long_name'))
    print(hasattr(p, 'long_name'))
    # Datasets will be closed at the end of your script, but can also be
    # closed manually.
    dataset.close()
