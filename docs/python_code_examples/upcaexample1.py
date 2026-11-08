#!/usr/bin/env python3
"""Multi-way PCA: 15 objects x 4 features measured in 3 conditions (layers)."""

import math

import libscientific

# a list of 3 layers, each 15 objects x 4 features
# (a numpy array of shape (3, 15, 4) works too)
x = [
    [
        [3.0 * (j + 1) + math.sin(0.9 * (i + 1)) * (1 + 0.5 * j) * (1 + 0.7 * k) for j in range(4)]
        for i in range(15)
    ]
    for k in range(3)
]

model = libscientific.upca.UPCA(npc=1, scaling=1)
model.fit(x)
print("Explained variance (%):", model.get_exp_variance())
print("Loading matrix of component 0 (features x layers):", model.get_loadings()[0])

scores = model.predict(x)
print("Scores of the first 3 objects:", scores[:3])

reconstructed = model.reconstruct_original_tensor()
print("X[0][0][0] =", round(x[0][0][0], 4), "reconstructed:", round(reconstructed[0][0][0], 4))
