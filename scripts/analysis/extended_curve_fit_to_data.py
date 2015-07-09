"""
extended_curve_fit_to_data.py:
    Fit x,y data to a model using the Levenberg-Marquardt method.
    Includes use of parameter classes, possibility of ning, some
    error checking, and Kolmogorov-Smirnov test.

    The example provided is a fit of a Gaussian function
    to data in file gauss.txt.

    To fit your own data, you need to change:
    (1) def func(x,*p) to return the function you are trying to fit,
    (2) the name of the data file read in by numpy.loadtxt,
    (3) the initial p0 values in the scipy.optimize.curve_fit call.

    Copyright (c) 2011-2014 University of Toronto
    Last Modification:   24 September 2014 by David Bailey
        Several rebinning bugs fixed.
    Original Version:    5 October   2011 by David Bailey
    Contact: David Bailey <dbailey@physics.utoronto.ca>
                            (http://www.physics.utoronto.ca/~dbailey)
    License: Released under the MIT License; the full terms are this license
                are appended to the end of this file, and are also available
                at http://www.opensource.org/licenses/mit-license.php.

    With inspiration from:
        http://www.scipy.org/Cookbook/FittingData and
        http://www.phy.uct.ac.za/courses/python/examples/fitresonance.py

    Thanks to Wes Watters <wwatters@wellesley.edu> (16 September 2013) for
        pointing out problem with scipy cov matrix. (See WARNING below.)

    To Do: Contour Probability Plots,
            automatic freezing of irrelvant variables.

"""
from __future__ import print_function

import matplotlib               # http://matplotlib.sourceforge.net
# Explicitly choose normal TkAgg Backend. This can be deleted if TkAgg or
#      another suitable backend is chosen in the user matplotlibrc file. See
# http://matplotlib.sourceforge.net/faq/installing_faq.html#what-is-a-backend.
matplotlib.use('TKAgg')
import numpy                    # http://numpy.scipy.org/
import scipy                    # http://numpy.scipy.org/

import scipy.stats
import scipy.optimize
import scipy.special
from matplotlib import pyplot

# Define expanded parameter class and fit definition,
#   using ideas from www.scipy.org/Cookbook/FittingData example and
#              www.phy.uct.ac.za/courses/python/examples/fitresonance.py


class Parameter(object):

    """ Parameter for fitting. """

    def __init__(self, name, value):
        self.guess = value
        self.name = name
        self.value = value

    def set(self, value):
        self.value = value

    def __call__(self):
        return self.value

    def __len__(self):          # What is returned by len(name)
        try:
            return len(self.value)
        except:
            return 1

    def __str__(self):          # String returned by str or print
        return "('" + self.name + "'," + str(self.value) + ")"

    def __repr__(self):         # A string is returned by repr
        return "('" + self.name + "'," + str(self.value) + ")"


def fit(function, parameters, x, y, y_uncertainty=None, **kwargs):

    def f(params):
        [p.set(params[i])for i, p in enumerate(parameters)]
        # Minimize difference between data and function,
        #   normalized by data uncertainties
        return (y - function(x)) / y_uncertainty

    # if no y uncertainties, assign equal weights (1) to all points
    if y_uncertainty is None:
        y_uncertainty = numpy.ones(len(y), "float")

    # Return Levenberg-Marquardt nonlinear least squares fit
    #   Using leastsq instead of curve_fit gives more potential control
    #           and access to the fit.
    # For more informatin on leastsq, see
    #   docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.leastsq.html
    # Keyword arguments ("**kwargs") are just passed through to leastsq
    p_fit, cov, infodict, mesg, ier = scipy.optimize.leastsq(
        f, [p.value for p in parameters], full_output=True, **kwargs)
    # WARNING : Scipy seems to use non-standard poorly documented notation for cov,
    #   which misleads many people. See "cov_x" on
    #   http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.leastsq.html#scipy.optimize.leastsq
    #   (which underlies curve_fit) and also discussion at
    #   http://stackoverflow.com/questions/14854339/in-scipy-how-and-why-does-curve-fit-calculate-the-covariance-of-the-parameter-es.
    #   I think this agrees with @cdeil at http://nbviewer.ipython.org/5014170/.
    #   THANKS to Wes Watters <wwatters@wellesley.edu> for pointing this out to me (16 September 2013)
    #
    # Convert Scipy cov matrix to standard covariance matrix.
    chisq = sum(infodict["fvec"]**2).sum()
    dof = len(x) - len(p_fit)
    cov_standard = cov * dof / chisq
    # Load fitted values into p
    for i in range(len(p_fit)):
        parameters[i].set(p_fit[i])
    return parameters, cov_standard, infodict, mesg, ier

# define function to be fitted
#   Note: A couple of examples are given here, but normally you will
#           write your own function describing the model you are fitting to.


def gaussian(x):
    # A gaussian with
    #   central value:                  x0
    #   standard deviation:             sigma_x
    #   peak height above background:   y0
    #   constant background:            B
    return B() + y_peak() * numpy.exp(-(x - x0())**2 / (2 * x_width()**2))


def lorentzian(x):
    # A lorentzian peak with:
    #   Constant Background          : B
    #   Peak height above background : y_peak
    #   Central value                : x0
    #   Full Width at Half Maximum   : width
    return B() + (y_peak() / numpy.pi) / (1.0 + ((x - x0()) / x_width())**2)

# Choose function
func = gaussian

# Define the fit parameters, guessing some starting values
B = Parameter("B", 10)
y_peak = Parameter("y_peak", 40.0)
# This converges to a negative y_peak if x0=1170, but fits perfectly if x0=1171.
#   One can't assume that a fit will produce sensible results.
x0 = Parameter("x0", 1171.0)
# The width is defined differently depending on choice of function
if func.__name__ == "gaussian":
    x_width = Parameter("sigma_x", 1.0)
elif func.__name__ == "lorentzian":
    x_width = Parameter("fwhm_x", 2.0)
else:
    print("\nInvalid width parameter name.\n")


p = [B, y_peak, x0, x_width]

# Read Data
#   unpack=True transposes lines in columns, e.g. if each line has an x, y, z
#               triplet, then x, y, z = numpy.loadtxt will load the data into
#               x[i], y[i], z[i] instead of xyz[i]
x_data, y_data, dy_data = numpy.loadtxt(
    '../../sample_data/gauss.txt', unpack=True, skiprows=0)

'''Problem Note: If you have created your own data file, which looks fine
    but produces "ValueError: too many values to unpack" when read with
    above statement, then it may have invisible characters
    or non-standard linebreaks, e.g. Mac CR instead
    of Windows style CRLF. This may be fixed by opening the data file in
    a text editor and resaving it with standard line breaks.
'''

# Rebin data if desired


def rebinned(rebin, x, y, dy):
    # Use Bresenham's algorithm to rebin data into (almost) equal slices.
    #  (see http://en.wikipedia.org/wiki/Bresenham%27s_line_algorithm)
    import numpy
    # Note type of input data, so same type can be returned
    data_type = type(x)
    # convert input data to numpy array (does nothing if already true)
    x, y, dy = numpy.array(x), numpy.array(y), numpy.array(dy)
    # initialize rebinned data and error
    x_rebin, y_rebin, dy_rebin = [], [], []
    # initialize data index for beginning of current rebin
    bin_start = 0
    # initialize rebinned index+1
    rebin_index = 1
    while bin_start < len(y):
        bin_end = min(int(round(rebin * rebin_index)), len(y))
        x_rebin.append(sum(x[bin_start:bin_end]) / (bin_end - bin_start))
        y_rebin.append(sum(y[bin_start:bin_end] / (bin_end - bin_start)))
        # Uncertainties combined as square-root of sum of squares
        dy_rebin.append(
            (sum(dy[bin_start:bin_end]**2))**0.5 / (bin_end - bin_start))
        bin_start = bin_end
        rebin_index += 1
    if data_type == numpy.ndarray:
        return numpy.array(x_rebin), numpy.array(y_rebin), numpy.array(dy_rebin)
    else:
        return x_rebin, y_rebin, dy_rebin

# Number of bins to be combined into each new bin
rebin = 1
# rebin data
x_data, y_data, dy_data = rebinned(rebin, x_data, y_data, dy_data)

# Record initial function guess
#   Create many finely spaced points for function plot
#           linspace(start,stop,num)
#               returns num evenly spaced samples over interval [start, stop]
x_func = numpy.linspace(min(x_data), max(x_data), 150)
initial_plot = func(x_func)

# Do fit using Levenberg-Marquardt
p, cov, infodict, mesg, ier = fit(func, p, x_data, y_data, dy_data)
#   Check for unsuccessful fit
#   	(see http://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.leastsq.html)
if ier > 4 or ier < 1:
    print("!!!! UNSUCCESSFULL FIT !!!!")
    print("Error Code (ier) = ", ier)
    print(mesg)
    print("!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
    cov = numpy.zeros((len(p), len(p)))

# calculate final chi square and degrees of freedom
chisq = sum(infodict["fvec"]**2)
dof = len(x_data) - len(p)

# calculate fit values and residuals
for i in enumerate(x_data):
    y_fit = func(x_data)
    y_data_res = y_data - y_fit

# Kolmogorov-Smirnov test is binning independent test of whether two
#       distributions are statistically different.
#     (e.g. see http://www.physics.csbsju.edu/stats/KS-test.html)
y_fit_sum = sum(y_fit)          # sum fit values
y_data_sum = sum(y_data)        # sum data values
# create test statistic array for all data points
D_ks = [0] * len(y_data_res)
y_fit_ks = 0                    # zero fit integral
y_data_ks = 0                   # zero data integral
for i in range(len(x_data)):
    y_fit_ks += y_fit[i]      # fit integral up to point i
    y_data_ks += y_data[i]    # data integral up to point i
    # normalized difference
    D_ks[i] = y_data_ks / y_data_sum - y_fit_ks / y_fit_sum
# maximum integrated difference between data and fit
Dn = max(max(D_ks), -min(D_ks))

# Plot

# create figure with light gray background
fig = pyplot.figure(facecolor="0.98")
# 3 rows, 1 column, subplot 1
#   3 rows are declared, but there are only 2 plots; this leaves room for text
#       in the empty 3rd  row
fit = fig.add_subplot(311)
fit.set_xticklabels(())  # remove tick labels from upper plot (for clean look)
# plot data
fit.errorbar(x_data, y_data, yerr=dy_data, fmt='b+', label="Data")
#   fmt='b+' creates blue crosses for data
#       (see http://matplotlib.sourceforge.net/api/pyplot_api.html)
# plot fit
#   draw starting guess as dashed green line ('r-')
fit.plot(x_func, initial_plot, 'g-', label="Start", linestyle="--")
#   draw final fit as solid red line ('r-')
fit.plot(x_func, func(x_func), 'r-', label="Fit")
pyplot.ylabel("Counts")
# Set fontsize for y tickmarks
pyplot.setp(fit.get_yticklabels(), fontsize="small")
fit.legend(loc='upper left',
           prop=matplotlib.font_manager.FontProperties(size="small"))

# separate plot to show residuals
residuals = fig.add_subplot(312)  # 3 rows, 1 column, subplot 2
residuals.errorbar(x_data, y_data_res, yerr=dy_data, fmt='b+',
                   label="Residuals")
# draw horizontal line at 0 on vertical axis
residuals.axhline(y=0, color="r")
pyplot.xlabel("Energy [KeV]")
pyplot.setp(residuals.get_xticklabels(), fontsize="small")
pyplot.setp(residuals.get_yticklabels(), fontsize="small")
residuals.legend(loc='upper left',
                 prop=matplotlib.font_manager.FontProperties(size="small"))

# These data look better if 'plain', not scientific, notation is used,
#   and if the tick labels are not offset by a constant (as is done by default).
#   Note: This only works for matplotlib version 1.0 and newer, so it is
#           enclosed in a "try" to avoid errors.
try:
    pyplot.ticklabel_format(style='plain', useOffset=False, axis='x')
except:
    print(
        "Your plots will look nicer if you upgrade to the newest matplotlib.")

# Output statistics to terminal and plot (using pyplot.figtext)
#   Only selected statistics are printed to plot
#  Note: If the fit is poor, i.e. chisq/dof is large, the uncertainties
#   are scaled up. If the fit is too good, i.e. chisq/dof << 1, it suggests
#   that the uncertainties have been overestimated, but the uncertainties
#   are not scaled down.
print("******** RESULTS FROM FIT ******** (by extended_curve_fit_to_data.py)")
print("Fit Function: ", func.__name__)
print("\nDegrees of Freedom (dof): {0:8.4g}".format(dof))
print("\nFitted parameters :")
try:
    for i in range(len(p)):
        print("  {0:10s} {1:12f} +/- {2:10f}     (Initial Guess: {3:6g})".\
              format(p[i].name, p[i].value,
                     numpy.sqrt(cov[i, i]) * max(1, numpy.sqrt(chisq / dof)), p[i].guess))

    print("\nCorrelation matrix")
    print(16 * (" "), end=' ')
    for i in range(len(p)):
        print("{0:10s}".format(p[i].name), end=' ')

    for i, row in enumerate(cov):
        print("\n    {0:10s}".format(p[i].name), end=' ')
        for j in range(len(p)):
            print("{0:>10f}".format(cov[i, j] / numpy.sqrt(cov[i, i] * cov[j, j])), end=' ')

    print("\n\nChi-Squared Statistics : ")
    print("  Chi-squared (chisq)          : {0:8.4g}".format(chisq))
    print("  Reduced chisq (chisq/dof)    : {0:8.4g}".format(chisq / dof))
    # The "cdf" is the probability that if the model is correct, and the
    #	parameters have the values output by the fit, and that the uncertainties,
    #	are normally distributed (i.e. Gaussian), that the experiment
    #	would result in a chi-squared greater than the value
    #	observed in this fit.
    cdf = scipy.special.chdtrc(dof, chisq)
    print("  Cumulative Probability (cdf) : {0:8.4g}%".format(
        100 * cdf))
    if cdf < 0.05:
        print("\nNOTE: This does not appear to be a great fit, so the")
        print("      parameter uncertainties may be underestimated.")
    elif cdf > 0.95:
        print("\nNOTE: This fit seems better than expected, so the")
        print("      data uncertainties may have been overestimated.")

    print("\nKolmogorov-Smirnov Statistics : ")
    print("  Integrals: Data, Fit : ", y_data_sum, y_fit_sum)
    print("  +/- maximum integrated differences : ", max(D_ks), min(D_ks))
    print("  Dsum, Dn = ", y_data_sum * Dn, Dn)
    # Kolmogorov-Smirnov one sample test
    #   This test is only valid if the uncertainties are statistical, i.e. root(n),
    #	on each data point.
    #   Also, since the parameters are determined from the data, the KS test
    #	statistics are even then not exactly valid.
    #   It is not clear if this scipy cdf calculation can be trusted in all
    #   cases, but it is suitable for flagging discrepancies.
    print("  Cumulative Probability (cdf) ",\
        1.0 - float(scipy.stats.ksone.cdf(Dn, y_data_sum)))

    # print selected information in empty 3rd plot row
    pyplot.figtext(0.05, 0.25, "Converged with ChiSq = " + str(chisq) + ", DOF = " +
                   str(dof) + ", CDF = " + str(100 * scipy.special.chdtrc(dof, chisq)) + "%")
    pyplot.figtext(0.05, 0.21, "Kolmogorov-Smirnov (Approximate): Sum = " +
                   '{0:.5g}'.format(float(y_fit_sum)) + ", Dn = " +
                   '{0:.5g}'.format(float(Dn)) + ", CDF = " +
                   '{0:.5g}'.format(
                       100. - 100. * float(scipy.stats.ksone.cdf(Dn, y_data_sum)))
                   + "%")
    for i, pmin in enumerate(p):
        pyplot.figtext(0.08, 0.16 - i * 0.03, (p[i].name).ljust(18) + " = " +
                       str(p[i].value).ljust(1) + " +/- " +
                       str(numpy.sqrt(cov[i, i])
                           * max(1, numpy.sqrt(chisq / dof))),
                       fontdict=None, family='monospace')
# If fit is bad, the above block will produce a python TypeError which is
#	caught by this try-except.
#   (Python leans to asking forgiveness afterwards, not permission before.)
# Note: Including family="Monospace" in the above figtext call will
#       produce nicer looking output, but can cause problems with
#       some older python installations.
except TypeError:
    print("**** BAD FIT ****")
    print("Parameters were: {0}".format(p))
    print("Uncertainties not calculated.")
    print("Try a different initial guess for the fit parameters.")
    pyplot.figtext(0.05, 0.25, "BAD FIT.  Guess again.")

pyplot.show()  # Show completed plot

# End extended_curve_fit_to_data.py

"""
Full text of MIT License:

    Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""
