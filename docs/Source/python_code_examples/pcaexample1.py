#!/usr/bin/env python3

import libscientific
import random

def mx_to_video(m, decimals=5):
    for row in m:
        print("\t".join([str(round(x, decimals)) for x in row]))

random.seed(123456)

# Create a random matrix of 10 objects and 4 features
a = [[random.random() for j in range(4)] for i in range(10)]
print("Original Matrix")
mx_to_video(a)

# Compute 2 Principal components using the UV scaling (unit variance scaling)
model = libscientific.pca.PCA(scaling=1, npc=2)
# Fit the model
model.fit(a)

# Show the scores
print("Showing the PCA scores")
scores = model.get_scores()
mx_to_video(scores, 3)

# Show the loadings
print("Showing the PCA loadings")
loadings = model.get_loadings()
mx_to_video(loadings, 3)

# Show the explained variance
print(model.get_exp_variance())

# Show the loadings
print("Predict/Project new data into the PCA model")
p_scores = model.predict(a)
mx_to_video(p_scores)

# Reconstruct the original PCA matrix from the 2 principal components
print("Reconstruct the original PCA matrix using the PCA Model")
ra = model.reconstruct_original_matrix()
mx_to_video(ra)
