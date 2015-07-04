# -*- coding: utf-8 -*-
"""
main_block.py: this script describes what the main block is and why you should
use it.

This code is released into the public domain. For full details, see
UNLICENSE at the root of this project. In essence, you can copy and use this
code for whatever purpose, with or without attribution.

If you have any questions about the code or could use any sort of help using
it, feel free to e-mail me at mcgibbon (at) uw {dot} edu.
"""
# Always import your packages at the top of your script, so when someone reads
# your code they know what you've imported from other files
import numpy as np

# below that, define any global variables you want access to. These should
# be things that don't really change, that you might want to import from
# other files, like scientific constants
Rd = 287.


# below that, define any functions or classes. See other files for more details
# about functions and classes
def foo(bar):
    """Takes in an integer, and returns the integer that is one greater."""
    return bar + 1

# then it's time for the main block
# this is where all of your "script" code should go. You should avoid putting
# scripts outside of main blocks. I can think of at least two good reasons for
# this:
#     it makes it clear where you're done defining things and where your
#     script begins, and what you consider to be "shared" (importable) vs.
#     script-specific code
# and
#     if you import anything from this file in another file, the import will
#     not run anything inside the if __name__ == '__main__' block. This lets
#     you import functionality from your script to use in other scripts,
#     without having to run all the analysis (for instance) in that file.
if __name__ == '__main__':
    bar = 10.
    while bar < 20:
        bar = foo(bar)
        print(bar)
