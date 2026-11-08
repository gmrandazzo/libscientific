"""Regression tests for the model fixes and the new model bindings.

Each test documents a behaviour that used to be wrong:
- refitting a model must give the same result as a new model;
- PLS predictions with nlv_ must return every response for multi-Y models;
- PLS beta coefficients must refuse multi-Y models;
- PLS VIP: the average of the squared values is 1;
- PCA/CPCA: decomposition and recomposition give back the data;
- LDA with classes numbered from 1;
- MLR, LDA, UPCA, UPLS bindings, including save/load.
"""

import math
import os

import numpy as np
import pytest

from libscientific.cpca import CPCA
from libscientific.lda import LDA
from libscientific.mlr import MLR
from libscientific.model_io import import_auto, import_model
from libscientific.pca import PCA
from libscientific.pls import PLS
from libscientific.upca import UPCA
from libscientific.upls import UPLS
from libscientific.vectlist import dvector_list_tolist

RNG = np.random.default_rng(42)
A = RNG.normal(size=(15, 4))
B = RNG.normal(size=(15, 4)) + 5.0
YA = RNG.normal(size=(15, 2))
YB = RNG.normal(size=(15, 2)) + 3.0


def test_pca_refit_equals_fresh_model():
    refit = PCA(npc=3)
    refit.fit(A)
    refit.fit(B)
    fresh = PCA(npc=3)
    fresh.fit(B)
    assert np.allclose(refit.get_scores(), fresh.get_scores())


def test_pca_reconstruction_is_exact_with_all_components():
    model = PCA(scaling=1, npc=4)
    model.fit(B)
    assert np.allclose(model.reconstruct_original_matrix(), B, atol=1e-8)
    assert sum(model.get_exp_variance()) == pytest.approx(100.0, abs=1e-6)


def test_pca_reconstruct_from_predicted_scores():
    model = PCA(npc=4)
    model.fit(B)
    new = B[:3] + 0.5
    assert np.allclose(model.reconstruct_original_matrix(scores_input=model.predict(new)), new)


def test_pca_npc_is_capped_to_features():
    model = PCA(npc=10)
    model.fit(A)
    assert model.npc == 4


def test_pls_refit_equals_fresh_model():
    refit = PLS(nlv=2)
    refit.fit(A, YA)
    refit.fit(B, YB)
    fresh = PLS(nlv=2)
    fresh.fit(B, YB)
    assert np.allclose(refit.predict(B)[0], fresh.predict(B)[0])


def test_pls_predict_nlv_multi_y():
    model = PLS(nlv=3)
    model.fit(A, YA)
    y_all, _ = model.predict(A)
    y_2, _ = model.predict(A, nlv_=2)
    assert np.shape(y_2) == (15, 2)
    assert np.allclose(y_2, np.array(y_all)[:, 2:4])
    with pytest.raises(ValueError):
        model.predict(A, nlv_=4)


def test_pls_beta_coefficients_single_y_only():
    model = PLS(nlv=2)
    model.fit(A, YA)
    with pytest.raises(ValueError):
        model.get_beta_coefficients(2)
    model.fit(A, YA[:, [0]])
    assert len(model.get_beta_coefficients(2)) == 4


@pytest.mark.parametrize("nlv", [1, 2, 3])
def test_pls_vip_mean_square_is_one(nlv):
    y = A[:, [0]] * 2.0 - A[:, [2]]
    model = PLS(nlv=nlv)
    model.fit(A, y)
    vip = np.array(model.get_vip())
    assert vip.shape == (4,)
    assert np.mean(vip**2) == pytest.approx(1.0, abs=1e-8)


def test_cpca_refit_and_recomposition():
    lat = RNG.normal(size=(20, 2))
    blocks = [lat @ RNG.normal(size=(2, 3)) + 1.0, lat @ RNG.normal(size=(2, 4)) - 2.0]
    model = CPCA(scaling=0, npc=2)
    model.fit([A[:, :3], A])  # first fit on other data
    model.fit(blocks)
    assert len(model.get_total_exp_variance()) == 2
    assert sum(model.get_total_exp_variance()) == pytest.approx(100.0, abs=1e-4)
    t = np.array(model.get_super_scores())
    for b, (block, loadings) in enumerate(zip(blocks, model.get_block_loadings())):
        avg = np.array(dvector_list_tolist(model.model.contents.colaverage)[b])
        assert np.allclose(t @ np.array(loadings).T + avg, block, atol=1e-6)


def test_mlr_fit_predict_refit_and_save(tmp_path):
    x = [[1.0, 0.0], [2.0, 1.0], [3.0, 4.0], [4.0, 2.0], [5.0, 3.0]]
    y = [[1.0 + 2.0 * a - 0.5 * b] for a, b in x]
    model = MLR()
    model.fit(x, y)
    model.fit(x, y)  # refit
    assert np.allclose(model.get_coefficients(), [[1.0], [2.0], [-0.5]])
    assert model.get_r2() == pytest.approx([1.0])
    assert np.allclose(model.predict([[10.0, 2.0], [0.0, 0.0], [1.0, 1.0]]), [[20.0], [1.0], [2.5]])
    with pytest.raises(ValueError):
        model.predict([[1.0]])
    for ext in ("lsm", "sqlite3"):
        path = str(tmp_path / f"mlr.{ext}")
        model.save(path)
        loaded = MLR()
        loaded.load(path)
        assert loaded.predict([[3.0, 3.0]]) == model.predict([[3.0, 3.0]])


@pytest.mark.parametrize("first_label", [0, 1])
def test_lda_labels_and_refit(first_label, tmp_path):
    x = [
        [1.0, 1.0],
        [1.5, 2.0],
        [2.0, 1.5],
        [6.0, 5.0],
        [6.5, 6.0],
        [7.0, 5.5],
        [1.0, 6.0],
        [1.5, 7.0],
        [2.0, 6.5],
    ]
    y = [first_label + i // 3 for i in range(9)]
    model = LDA()
    model.fit(np.array(x) + 10.0, y)
    model.fit(x, y)  # refit
    assert model.classes == [first_label, first_label + 1, first_label + 2]
    assert model.predict(x) == y
    assert np.shape(model.predict_scores(x)) == (9, 3)
    assert np.shape(model.transform(x)) == (9, 2)
    path = str(tmp_path / "lda.lsm")
    model.save(path)
    assert import_model(path).predict(x) == y


def _trilinear(n=15, f=4, k=3):
    return [
        [
            [
                3.0 * (j + 1) + math.sin(0.9 * (i + 1)) * (1 + 0.5 * j) * (1 + 0.7 * kk)
                for j in range(f)
            ]
            for i in range(n)
        ]
        for kk in range(k)
    ]


@pytest.mark.parametrize("scaling", [0, 1])
def test_upca_recomposition_predict_and_refit(scaling, tmp_path):
    x = _trilinear()
    model = UPCA(npc=1, scaling=scaling)
    model.fit(np.array(x) + 3.0)
    model.fit(x)  # refit
    assert model.get_exp_variance() == pytest.approx([100.0])
    assert np.allclose(model.reconstruct_original_tensor(), x, atol=1e-8)
    assert np.allclose(model.predict(x), model.get_scores(), atol=1e-8)
    path = str(tmp_path / "upca.lsm")
    model.save(path)
    loaded = UPCA()
    loaded.load(path)
    assert np.allclose(loaded.predict(x), model.predict(x))


def test_upls_fit_predict_refit_and_save(tmp_path):
    x = RNG.normal(size=(2, 20, 4))
    y = (x[0] @ [1.0, 2.0, 0.0, 0.0] + x[1] @ [0.0, 0.0, 1.0, -1.0]).reshape(1, 20, 1)
    model = UPLS(nlv=3)
    model.fit(x + 3.0, y)
    model.fit(x, y)  # refit
    fresh = UPLS(nlv=3)
    fresh.fit(x, y)
    pred = np.array(model.predict(x))
    assert pred.shape == (1, 20, 1)
    assert np.allclose(pred, fresh.predict(x))
    r2 = 1 - ((pred - y) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    assert r2 > 0.95
    path = str(tmp_path / "upls.lsm")
    model.save(path)
    assert np.allclose(import_auto(path).predict(x), pred)
    assert os.path.getsize(path) > 0
