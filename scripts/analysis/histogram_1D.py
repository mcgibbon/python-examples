"""
histogram_1D.py
    
    A simple example to histogram data from a file and output it.
    
    To histogram your own data, you need to change:
    (1) name of file containing the data ("hist_data.txt" in this example)
    (2) select the column containing the data to be histogrammed.
    (3) any desired optional parameters or labels.

    The program produces a histogram, which is also saved as a text file
    (e.g. for input to a fitting program) and as a png graphic.

    Copyright (c) 2013 University of Toronto
    Original Version:   2 March  2013 by David Bailey
    Contact: David Bailey <dbailey@physics.utoronto.ca>
                            (http://www.physics.utoronto.ca/~dbailey)
    License: Released under the MIT License; the full terms are this license
                are appended to the end of this file, and are also available
                at http://www.opensource.org/licenses/mit-license.php.

    Note: The code is over-commented since it is is a pedagogical example for
    	students new to python and scipy.

"""
import numpy
from matplotlib import pyplot

## Load data from file (Change name to your data file)
#   Data file is a simple text file, with columns of numbers.
#   The "columns" do not need to be aligned.
columns = numpy.loadtxt("hist_data.txt", unpack=True)

# Plot and axis titles (optional, but very strongly recommended)
pyplot.title("Example Histogram")
pyplot.xlabel("x")
pyplot.ylabel("Entries")

## Histogram the selected data column
#    Only the column number is required, all other parameters are optional.
#       For all the pyplot.hist parameters and options, see
#           http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.hist      
bin_contents, bin_ends, patches = pyplot.hist(
                columns[0],          # Selected data column
                # All the following parameters are optional
                bins      = 100,      # number of hisogram bins
                # Minimum and maximum of data to be plotted
                range     = (2000., 20000.),
                label     = "Run 1", # label for this data
                normed    = False,   # normalize histogram to unity
                log       = True,    # Choose if y log axis
                histtype  = "bar",   # Type of histogram
                facecolor = "green", # bar colour
                edgecolor = "blue",  # bar border colour
                linewidth = 1,       # bar border width
                alpha     = 0.8,     # bar transparency
                )

print sum(bin_contents), " entries in histogram"

# A legend is optional, but is useful is more than one set of data plotted.
pyplot.legend()

# Save plot as a graphic file, if desired.
pyplot.savefig('histogram_demo',dpi=72)

## Write histogram to file ("hist.txt"), if desired, e.g. for fitting.
#   The format is centre-of-bin position, width of bin, number of counts in bin.
bin_width     = (bin_ends[1]-bin_ends[0])
bin_widths    = len(bin_contents)*[bin_width]
bin_positions = bin_ends[0:-1]+bin_width/2.0
numpy.savetxt("hist.txt", zip(bin_positions, bin_widths, bin_contents))

# Display data plot in terminal
pyplot.show()

##### End histogram_1D.py
