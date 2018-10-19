from sklearn.datasets import make_blobs
import numpy as np

import math

n_samples = 100
centers = [(-5, -5), (3,5)]

X, y = make_blobs(n_samples=n_samples, n_features=2, cluster_std=1.0, centers=centers, shuffle=False, random_state=42)

"""
mean = (1, 2)
cov = [[1, -2], [0, 1]]

X = []
for i in range(n_samples):
    X.append(np.random.multivariate_normal(mean, cov))

X = np.array(X)
"""
for row in X:
    print("%.3f %.3f" % (row[0], row[1]))

print("AAAAAA")
covX = np.cov(X, rowvar=0)
print(covX.shape)
invcovX = np.linalg.inv(covX)

u = np.mean(X, axis=0)

for row in X:
    ru = row-u
    r = ru.dot(invcovX)
    dst = math.sqrt(r.dot(ru))
    print("%.3f %.3f %.3f" % (r[0], r[1], dst))
