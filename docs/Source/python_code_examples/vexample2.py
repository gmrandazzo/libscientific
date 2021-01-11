#!/usr/bin/env python3
import scientific
from random import random

# Create a list of values that you whant to convert to a double vector
a = [random() for j in range(5)]
d = scientific.vector.DVector(a)
# print the output of the double vector d
print("orig vector")
d.debug()


# append the value 0.98765 at the end of d
d.append(0.98765)
print("append 0.98765 at the end")
d.debug()

# extend the vector d with more other values from a list
d.extend([0.4362, 0.34529, 0.99862])
print("extent the vector with 3 more values")
d.debug()
