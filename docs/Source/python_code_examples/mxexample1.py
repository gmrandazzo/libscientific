#!/usr/bin/env python3
import libscientific
from random import random

# Create a random list of list 
a = [[random() for j in range(2)] for i in range(10)]

# Convert the list of list matrix into a libscientific matrix
m = libscientific.matrix.Matrix(a)

# Get the value at row 1, column 1
print("Get value example")
print(m[1, 1])

# Modify the value at row 1, column 1
print("Set value example")
m[1, 1] = -2.
m.debug()


# Convert the matrix again to a list of list
mlst = m.tolist()
for row in mlst:
    print(row)
