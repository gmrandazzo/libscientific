#!/usr/bin/env python3
"""Save and load models: single model files (.lsm) and bundles (.lsb)."""

import os
import tempfile

from libscientific.model_io import ModelBundle, import_auto, import_model, inspect_model_file
from libscientific.pca import PCA
from libscientific.pls import PLS

x = [[1.0, 2.0, 3.0], [2.0, 1.0, 4.0], [3.0, 4.0, 2.0], [4.0, 3.0, 1.0], [5.0, 5.0, 0.5]]
y = [[1.0], [2.0], [3.0], [4.0], [5.0]]
outdir = tempfile.mkdtemp()

pca = PCA(scaling=1, npc=2)
pca.fit(x)
pls = PLS(nlv=2)
pls.fit(x, y)

# 1. single model: the save/load methods of every model
pca_path = os.path.join(outdir, "pca_model.lsm")
pca.save(pca_path)
loaded = PCA()
loaded.load(pca_path)
print("Same scores after loading:", loaded.predict(x) == pca.predict(x))

# ... or let the library create the right object
model = import_model(pca_path)
print("import_model created:", model)

# 2. several models in one bundle file
bundle_path = os.path.join(outdir, "project.lsb")
ModelBundle({"screening_pca": pca, "activity_pls": pls}).save(bundle_path)
print("File content:", inspect_model_file(bundle_path))

bundle = import_auto(bundle_path)  # a ModelBundle for multi-model files
print("Models in the bundle:", list(bundle.keys()))
y_pred, _ = bundle["activity_pls"].predict(x, nlv_=2)
print("Predictions of the loaded PLS:", [round(v[0], 3) for v in y_pred])
