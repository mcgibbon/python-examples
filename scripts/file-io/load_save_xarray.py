# -*- coding: utf-8 -*-
"""
load_netcdf_simple.py: an example on loading data from a netCDF file. Only
the most commonly-performed operations are covered in this file. For a more
extensive set of operations, see load_netcdf_full.py.

This code is released into the public domain. For full details, see
UNLICENSE at the root of this project. In essence, you can copy and use this
code for whatever purpose, with or without attribution.

If you have any questions about the code or could use any sort of help using
it, feel free to e-mail me at mcgibbon (at) uw {dot} edu.
"""
# we will be using the xarray library
import xarray as xr
from glob import glob  # used to get file paths

# This 'if' statement should always be put at the start of your
# 'script' section after any function or class definitions, or global
# variables.
if __name__ == '__main__':
    dataset = xr.open_dataset('../../sample_data/ARM_sounding_example.nc')
    # data is stored in the dataset.data_vars dictionary
    # its keys are variable names, and values are Variable objects
    print('\n----- Variable object:')
    print(dataset.data_vars['pres'])

    print('\n----- Array object:')
    print(dataset.data_vars['pres'].values)
    
    # You can open multiple files into a single dataset object
    filenames = glob('../../sample_data/sgpmet*.cdf')
    print('Loading filenames: {}'.format(filenames))
    dataset = xr.open_mfdataset(filenames)

    # We can change the units of Pressure and save the dataset out to a new file
    print(dataset.data_vars.keys())
    print(dataset.data_vars['atmos_pressure'])  # who uses kPa anyways???

    # Wrong:
    # dataset.data_vars['atmos_pressure'] *= 1000.
    # dataset['atmos_pressure'].units = 'Pa'
    # Right:
    dataset['atmos_pressure'].values *= 1000.
    dataset['atmos_pressure']['units'] = 'Pa'

    # This dataset also failed to identify 'time_offset' as the time axis. Let's fix it.
    print('\n----- Old time variable:')
    print(dataset['time'])

    dataset['time'] = dataset['time_offset']
    # Wrong:
    # dataset.drop(['time_offset', 'base_time')
    # Right:
    dataset = dataset.drop(['time_offset', 'base_time'])

    print('\n----- New time variable:')
    print(dataset['time'])

    # We can resample the data to be hourly instead of 1-minute:
    dataset = dataset.resample('1H', dim='time', how='mean')
    print('\n----- Resampled hourly dataset:')
    print(dataset)

    dataset.to_netcdf('output.nc')  # saves to current directory
