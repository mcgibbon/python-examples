# -*- coding: utf-8 -*-
"""
importing.py: an example script on how to import packages in Python.

This code is released into the public domain. For full details, see
UNLICENSE at the root of this project. In essence, you can copy and use this
code for whatever purpose, with or without attribution.

If you have any questions about the code or could use any sort of help using
it, feel free to e-mail me at mcgibbon (at) uw {dot} edu.
"""
# To do things in Python, you will often want to load functions and other
# objects from modules. Some commonly used modules in geoscience are numpy,
# netCDF4, and matplotlib.

# There are a couple different conventions for importing modules
# One is import [packagename]
import numpy
print(numpy.exp(1))  # prints e

# To make code more concise, you might want to shorten the package name with
# import [packagename] as [alias]
import numpy as np
print(np.exp(1))  # also prints e

# The normal convention for matplotlib.pyplot is plt:
import matplotlib.pyplot as plt

# To avoid referencing the package a function comes from when you use a
# function (or class), you can import the name itself with 
# from [packagename] import [object1], [object2], ...
from numpy import exp
# or
from numpy import log, zeros
print(log(2))  # prints the natural log of 2

# You can also use aliases for functions:
from numpy import log as logarithm
print(logarithm(2))  # also prints the natural log of 2
