#!/usr/bin/env python3
"""Linear Discriminant Analysis: classify objects of 3 classes labelled 1, 2, 3."""

import libscientific

x = [[10.0 * (i // 4) + 0.1 * (i % 4), -5.0 * (i // 4) + 0.2 * ((i * 3) % 4)] for i in range(12)]
y = [1 + i // 4 for i in range(12)]  # labels 1, 2, 3 (0, 1, 2 works too)

model = libscientific.lda.LDA()
model.fit(x, y)
print("Classes:", model.classes, "priors:", model.get_priors())

print("Predicted classes:", model.predict(x))
print("Discriminant scores of the first object:", model.predict_scores(x[:1]))
print("Projection on the discriminant directions:", model.transform(x[:1]))
