# -*- coding: utf-8 -*-
"""
load_UW_text.py: an example for loading UW formatted text data using the
genfromtxt command in numpy. For data formatting see the file loaded in
the script below.

This code is released into the public domain. For full details, see
UNLICENSE at the root of this project. In essence, you can copy and use this
code for whatever purpose, with or without attribution.

If you have any questions about the code or could use any sort of help using
it, feel free to e-mail me at mcgibbon (at) uw {dot} edu.
"""
import numpy as np
# "datetime" is an object in the "datetime" module"
from datetime import datetime

if __name__ == '__main__':
    # name of an example data file
    filename = '../../UW_apu30_2015_0328.txt'
    # first argument is name of file
    # skip_header says how many lines to skip at the start of the file,
    #     default is 0.
    # usecols takes in a list of integers that say which columns to load.
    #     numpy can't load text data, so we skip the first two columns.
    data_array = np.genfromtxt(filename, skip_header=1,
                               usecols=range(2, 17))
    # we could now use the data directly from this array, but there are a few
    # ways to make it easier to read your code

    # if you know the names of the columns, you can make a tuple in your code
    column_names = ['Year', 'MM', 'DD', 'Jday', 'hh', 'mm', 'Conc', 'LWC',
                    'Rain', 'dBZ', 'Dm', 'DMax', 'Sigma_M', 'Wx_Code4677']
    # OR
    # if you don't know the column names, you can read them from the text file
    # "with" will open the file at the start of the block, and close it when
    #     you exit the block
    with open(filename, 'r') as f:
        # This reads the first line of the file, including the newline
        # character at the end of the line
        header_line = f.readline()
    # we want to remove the newline character at the end of the line
    # note that we could also have put the [:-1] at the end of the f.readline()
    #     command
    header_line = header_line[:-1]
    # this splits the string by whitespace into a list
    column_names = header_line.split()
    # remove the first two headers, since we did not load them
    column_names = column_names[2:]

    # Note that column_names (for this file) is the same whether we
    # use dictionary comprehension to make a data dictionary
    data_dict = {column_names[i]: data_array[:, i]
                 for i in range(len(column_names))}
    # Now we have a dictionary and can access values with expressions like
    # data_dict['Conc'][0], data_dict['Sigma_M'][10:20], etc.

    # In Python you would usually rather use what are called "datetime" objects
    # than floats to represent time. To generate datetimes:
    # initialize a list in which to put our datetimes
    datetime_list = []
    for i in range(len(data_dict['Year'])):
        # make the new datetime
        new_datetime = datetime(data_dict['Year'][i], data_dict['MM'][i],
                                data_dict['DD'][i], data_dict['hh'][i],
                                data_dict['mm'][i])
        # add the datetime to our list
        datetime_list.append(new_datetime)
    # If you like, you can make an array of datetimes
    datetime_array = np.array(datetime_list)
    # You could add it to your data dictionary, if you like
    data_dict['datetime'] = datetime_array
