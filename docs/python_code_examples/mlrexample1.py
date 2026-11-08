#!/usr/bin/env python3
"""Multiple Linear Regression: fit y = 1 + 2 x1 - 0.5 x2 and predict new objects."""

import libscientific

x = [[float(i), float((i * 7) % 4)] for i in range(6)]
y = [[1.0 + 2.0 * a - 0.5 * b] for a, b in x]

model = libscientific.mlr.MLR()
model.fit(x, y)

# coefficients: row 0 is the intercept, then one row per feature
print("Coefficients (intercept, x1, x2):", [round(row[0], 3) for row in model.get_coefficients()])
print("R2 of the fit:", model.get_r2())

print("Predictions:", model.predict([[10.0, 2.0], [-1.0, 0.0]]))
