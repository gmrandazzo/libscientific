#!/usr/bin/env python3
import scientific
from random import random

# Create a list of values that you whant to convert to a double vector
a = [random() for j in range(5)]

# Transform the list a into a double vector d
d = scientific.vector.DVector(a)

# Just print to video the content of vector d
d.debug()

# If you want to catch the value in position 1
print(d[1])

# If you want to modify the value in position 1
d[1] = -2

#If you want to get back the result as a "list" 
dlst = d.tolist()

for item in dlst:
    print(item)
