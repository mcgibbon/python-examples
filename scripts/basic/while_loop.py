# -*- coding: utf-8 -*-
"""
while_loop.py: this script shows examples of while loops.

This code is released into the public domain. For full details, see
UNLICENSE at the root of this project. In essence, you can copy and use this
code for whatever purpose, with or without attribution.

If you have any questions about the code or could use any sort of help using
it, feel free to e-mail me at mcgibbon (at) uw {dot} edu.
"""

if __name__ == '__main__':
    # The basic format is while [condition]
    done = False
    i = 0
    while not done:
        if i > 10:
            done = True
        i += 1  # add 1 to i

    # perhaps cleaner than the above example is this:
    i = 0
    while i <= 10:
        i += 1

    # You can also break out of while loops with the break statement
    i = 0
    while True:
        if i > 10:
            break
        i += 1

    # as with for statements, you can use the else clause on while loops,
    # which will only trigger if the while loop is exited by the condition
    # being false, and not if you break from the loop
    i = 1
    j = 1
    while i < 10:
        if j > 20:
            print('doubling wins!')
            break
        i += 1  # add 1 to i
        j *= 2  # double j
    else:
        print('adding wins!')
