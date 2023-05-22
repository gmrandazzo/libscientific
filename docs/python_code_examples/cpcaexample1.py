#!/usr/bin/env python3

import libscientific
import random

def mx_to_video(m, decimals=5):
    for row in m:
        print("\t".join([str(round(x, decimals)) for x in row]))

def t_to_video(t):
    i = 1
    for m in t:
        print("Block: %d" % (i))
        mx_to_video(m, 3)
        i+=1

random.seed(123456)

# Create a random matrix of 10 objects and 4 features
a = [[[random.random() for j in range(4)] for i in range(10)] for k in range(4)]

print("Original Matrix")
t_to_video(a)

# Compute 2 Principal components using the UV scaling (unit variance scaling)
model = libscientific.cpca.CPCA(scaling=1, npc=2)
# Fit the model
model.fit(a)

# Show the super scores
print("Showing the CPCA super scores")
sscores = model.get_super_scores()
mx_to_video(sscores, 3)

# Show the super weights
print("Showing the CPCA super weights")
sweights = model.get_super_weights()
mx_to_video(sweights, 3)

# Show the block scores
print("Showing the CPCA block scores")
block_scores = model.get_block_scores()
t_to_video(block_scores)

# Show the block loadings
print("Showing the CPCA block loadings")
block_loadings = model.get_block_loadings()
t_to_video(block_loadings)

# Show the total variance explained by the super scores
print("Showing the CPCA total variance explained")
print(model.get_total_exp_variance())

# Predict/Project new data into the model
print("Project/Predict new data into the CPCA model")
p_ss, p_bs = model.predict(a)
print("Showing the predicted super scores")
mx_to_video(p_ss, 3)
print("Showing the predicted block scores")
t_to_video(p_bs)
