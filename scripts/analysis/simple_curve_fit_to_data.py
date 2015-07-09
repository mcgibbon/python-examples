""" simple_curve_fit_to_data.py
    A very simple example of scipy curve_fit to data from a file.

    Last modified 4 January by David Bailey <dbailey@physics.utoronto.ca>
        to be compatible with Python 2.7 and Python 3+
    Originally based in 2010 on example at
        http://www.scipy.org/scipy_Example_List.

    Requires data file: gauss.dat

    Note: The code is overcommented since it is is a pedagogical example for
        students new to python and scipy.
    Note: Scipy has a non-standard definition of cov, which means that
        parameter uncertainties here are incorrect if fit is poor.
        See more detailed discussion in extended_curve_fit_to_data.py

"""
from __future__ import print_function
import pylab
import numpy
from scipy.optimize import curve_fit

# Function to fit: 'x' is the independent variable(s), 'p' the parameter vector
#   Note:   A one line lambda function definition  can be used for very simple
#           functions, but using "def" always works.
#   Note:   "*p" unpacks p into its elements; needed for curvefit


def func(x, *p):
    return p[0] + p[1] * numpy.exp(-1 * (x - p[2])**2 / (2 * p[3]**2))

# Load data to fit
x, y_data, y_sigma = numpy.loadtxt('../../sample_data/gauss.txt', unpack=True)

# Fit function to data

# For information of curve_fit.py, see
#   http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
# This fits the function "func" to the data points (x, y_data) with y
#   uncertainties "y_sigma", and initial parameter values p0.
p, cov = curve_fit(func, x, y_data, p0=(10., 40., 1171., 1,), sigma=y_sigma)

# Output results

print("Covariance Matrix : \n", cov, "\n")
print("Estimated parameters: ", p)
try:
    print("Estimated uncertainties: ", numpy.sqrt(cov.diagonal()))
# If cov has not been calculated because of a bad fit, the above print
#   statement will cause a python AttributeError which is caught by
#   this try-except. This could be checked with an if statement, but
#   Python leans to asking forgiveness afterwards, not permission before.
except AttributeError:
    print("Not calculated; fit is bad.")

# Plot data as red circles, and fitted function as (default) line
pylab.plot(x, y_data, 'ro', x, func(x, *p))
# Add error bars on data as red crosses.
pylab.errorbar(x, y_data, yerr=y_sigma, fmt='r+')
pylab.show()

# End simple_curve_fit_to_data.py
