#!/usr/bin/env python3
"""Multi-way PLS: predict one response from 2 layers of 4 features."""

import libscientific

# X: 2 layers of 20 objects x 4 features; Y: 1 layer of 20 objects x 1 response
x = [[[((i * (j + 3) + k * 7) % 11) / 3.0 for j in range(4)] for i in range(20)] for k in range(2)]
y = [[[x[0][i][0] + 2.0 * x[0][i][1] + x[1][i][2]] for i in range(20)]]

model = libscientific.upls.UPLS(nlv=3, xscaling=1, yscaling=0)
model.fit(x, y)
print("X variance explained by each latent variable (%):", model.get_exp_variance())

y_pred = model.predict(x)  # list of Y layers
for i in range(3):
    print(f"object {i}: y = {y[0][i][0]:.3f}, predicted = {y_pred[0][i][0]:.3f}")
