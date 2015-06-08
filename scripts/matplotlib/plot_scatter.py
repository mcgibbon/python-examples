# -*- coding: utf-8 -*-
"""
plot_scatter.py: an example for making a scatterplot using matplotlib.

This code is released into the public domain. For full details, see
UNLICENSE at the root of this project. In essence, you can copy and use this
code for whatever purpose, with or without attribution.

If you have any questions about the code or could use any sort of help using
it, feel free to e-mail me at mcgibbon (at) uw {dot} edu.
"""
# we will use numpy to load data
import numpy as np
# we will use matplotlib for plotting
import matplotlib.pyplot as plt
# we will use datetime objects to keep track of time
from datetime import datetime


# We will use this function to load some data to plot
# For more details on loading data, see examples in the "file-io" folder.
def load_UW_apu_text_data(filename):
    """
    Takes in a filename corresponding to a UW apu timeseries. Returns a
    dictionary whose keys are the headings for non-string fields in the file,
    and values are 1-dimensional numpy arrays containing the values in those
    fields. The dictionary also contains a key 'datetime' whose values are
    datetime objects corresponding to the time of that index.
    """
    # We skip the first two columns because they are text fields, and numpy
    # cannot load text data to arrays.
    data_array = np.genfromtxt(filename, skip_header=1,
                               usecols=range(2, 17))
    # We know the column names, so we can put them in our code
    column_names = ['Year', 'MM', 'DD', 'Jday', 'hh', 'mm', 'Conc', 'LWC',
                    'Rain', 'dBZ', 'Dm', 'DMax', 'Sigma_M', 'Wx_Code4677']

    # let's make a data dictionary using our column name information
    data_dict = {column_names[i]: data_array[:, i]
                 for i in range(len(column_names))}
    # Now we have a dictionary and can access values with expressions like
    # data_dict['Conc'][0], data_dict['Sigma_M'][10:20], etc.

    # let's generate a datetime field
    datetime_list = []
    for i in range(len(data_dict['Year'])):
        # make the new datetime
        new_datetime = datetime(*[int(data_dict[key][i]) for key in
                                  ('Year', 'MM', 'DD', 'hh', 'mm')])
        # add the datetime to our list
        datetime_list.append(new_datetime)
    # we will pack it into an array and add it to our data dictionary
    datetime_array = np.array(datetime_list)
    data_dict['datetime'] = datetime_array
    # return the data dictionary
    return data_dict

# This 'if' statement should always be put at the start of your
# 'script' section after any function or class definitions, or global
# variables.
if __name__ == '__main__':
    # First we need to load data to plot.
    data_dict = load_UW_apu_text_data(
        '../../sample_data/UW_apu30_2015_0328.txt')

    # One way to plot is with the Matlab-like interface, where you use
    # the same commands as you would in Matlab but from the matplotlib.pyplot
    # library
    # initialize a new figure
    plt.figure()
    # plot LWC against Rain
    plt.plot(data_dict['LWC'], data_dict['Rain'], 'x')
    # Set the title
    plt.title('LWC against Rain')
    # label the x and y axes
    plt.xlabel('LWC')
    plt.ylabel('Rain')
    # set x and y plottinglimits
    plt.xlim(0, 4000)
    plt.ylim(0, 2)
    # display the figure
    plt.show()

    # A second way to plot is with the object-oriented subplot() interface.
    # If you pass in no arguments, matplotlib assumes you want just 1 plot.
    fig, ax = plt.subplots()
    # plot LWC against Rain
    ax.plot(data_dict['LWC'], data_dict['Rain'], 'x')
    # Set the title
    ax.set_title('LWC against Rain')
    # label the x and y axes
    ax.set_xlabel('LWC')
    ax.set_ylabel('Rain')
    # set x and y plottinglimits
    ax.set_xlim(0, 4000)
    ax.set_ylim(0, 2)
    plt.show()
