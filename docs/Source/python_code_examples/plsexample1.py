#!/usr/bin/env python3

import scientific
import random

def mx_to_video(m, decimals=5):
    for row in m:
        print("\t".join([str(round(x, decimals)) for x in row]))

random.seed(123456)
x = [[random.random() for j in range(4)] for i in range(10)]
y = [[random.random() for j in range(1)] for i in range(10)]
xp = [[random.random() for j in range(4)] for i in range(10)]

print("Original Matrix")
print("X")
mx_to_video(x)
print("Y")
mx_to_video(y)
print("XP")
mx_to_video(xp)
print("Computing PLS ...")
model = scientific.pls.PLS(nlv=2, xscaling=1, yscaling=0)
model.fit(x, y)
print("Showing the PLS T scores")
tscores = model.get_tscores()
mx_to_video(tscores, 3)

print("Showing the PLS U scores")
uscores = model.get_uscores()
mx_to_video(uscores, 3)

print("Showing the PLS P loadings")
ploadings = model.get_ploadings()
mx_to_video(ploadings, 3)

print("Showing the X Variance")
print(model.get_exp_variance())


print("Predict XP")
py, pscores = model.predict(xp)
print("Predicted Y for all LVs")
mx_to_video(py, 3)
print("Predicted Scores")
mx_to_video(pscores, 3)
