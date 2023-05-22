#!/usr/bin/env python3
import libscientific
from random import random

# Create a random list of list 
a = [[[random() for j in range(2)] for i in range(10)] for k in range(3)]

# Convert the list of list of lists into a libscientific tensor
t = libscientific.tensor.Tensor(a)

# Get the value at row 1, column 1
print("Get value example")
print(t[1, 1, 1])

# Modify the value at row 1, column 1
print("Set value example")
t[1, 1, 1] = -2.
t.debug()


# Convert the matrix again to a list of list
tlst = t.tolist()
i = 1
for block in tlst:
    print("Block %d" % (i))
    for row in block:
        print(row)
    i+=1
