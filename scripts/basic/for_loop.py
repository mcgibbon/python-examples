# -*- coding: utf-8 -*-
"""
for_loop.py: this script shows examples of for loops.

This code is released into the public domain. For full details, see
UNLICENSE at the root of this project. In essence, you can copy and use this
code for whatever purpose, with or without attribution.

If you have any questions about the code or could use any sort of help using
it, feel free to e-mail me at mcgibbon (at) uw {dot} edu.
"""

if __name__ == '__main__':
    # The basic syntax for a for loop is for [name] in [iterator]
    for fruit in ('apples', 'pears', 'bananas'):
        print("I like {}".format(fruit))

    # You can skip to the next item in the iterator with "continue"
    for pet in ('dog', 'cat', 'mouse'):
        if pet == 'cat':
            continue
        print("I don't have a {}".format(pet))

    # You can break out of the for loop entirely with "break"
    for i in (0, 5, 10, 2):  # range(start, stop, step) works only for ints
        if i > 7:
            break
        print("I can count to {}".format(i))

    # You can also include an "else" clause that only runs if the for loop
    # ends *without* a break statement triggering
    for animal in ('lions', 'tigers', 'bears'):
        if animal == 'bears':
            print('Eep! Run away!')
            break
        print("I'm not afraid of {}".format(animal))
    else:
        print("I'm not afraid of anything!")

    # You can use tuple assignment in iterators
    for name, age in (('Alice', 10), ('Bob', 9)):
        print('Meet {}, who is {} years old'.format(name, age))

    # This also applies to the dictionary .items() function
    ages = {'Alice': 10, 'Bob': 9}
    for name, age in ages.items():
        print('Meet {}, who is {} years old'.format(name, age))

    # You could iterate through lists (or particularly arrays) with an index:
    cities = ['Toronto', 'Seattle', 'London']
    for i in range(len(cities)):
        print("I've been to {}".format(cities[i]))
