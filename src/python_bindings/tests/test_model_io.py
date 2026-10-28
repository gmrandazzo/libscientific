"""
test_model_io.py - Unit tests for Python model I/O bindings (.lsm / .lsb)
"""

import os
import random
import tempfile

import pytest

import libscientific.matrix as mx
from libscientific.model_io import (
    ModelBundle,
    ModelType,
    import_auto,
    import_model,
    import_model_bundle,
    inspect_model_file,
    is_bundle,
    save_model,
    save_model_bundle,
)
from libscientific.pca import PCA
from libscientific.pls import PLS


def matrix_sum(m_input):
    s = 0.0
    for row in m_input:
        s += sum(row)
    return s


@pytest.fixture
def sample_pca_data():
    random.seed(42)
    return [[random.random() * 10.0 + j for j in range(5)] for _ in range(20)]


@pytest.fixture
def sample_pls_data():
    random.seed(99)
    x = [[random.random() * 5.0 + j * 0.5 for j in range(6)] for _ in range(25)]
    # Dependent y values with linear relation plus noise
    y = [[row[0] * 1.5 + row[1] * 0.8 + 0.1 * random.random()] for row in x]
    return x, y


def test_pca_single_model_lsm(sample_pca_data):
    """Test saving and importing a single PCA model using .lsm format."""
    pca = PCA(scaling=1, npc=2)
    pca.fit(sample_pca_data)

    orig_scores = pca.get_scores()
    orig_loadings = pca.get_loadings()
    orig_varexp = pca.get_exp_variance()
    orig_pred = pca.predict(sample_pca_data[:5])

    with tempfile.TemporaryDirectory() as tmpdir:
        filepath = os.path.join(tmpdir, "test_pca.lsm")

        # 1. Save model
        save_model(filepath, pca, model_name="pca_v1")
        assert os.path.exists(filepath)

        # 2. Inspect metadata without full payload
        info = inspect_model_file(filepath)
        assert info["num_models"] == 1
        assert info["is_single"] is True
        assert info["primary_type"] == ModelType.PCA
        assert info["primary_name"] == "pca_v1"
        assert is_bundle(filepath) is False

        # 3. Import back into a new instance
        loaded_pca = import_model(filepath)
        assert isinstance(loaded_pca, PCA)
        assert loaded_pca.npc == 2

        # 4. Verify exact numerical match
        loaded_scores = loaded_pca.get_scores()
        loaded_loadings = loaded_pca.get_loadings()
        loaded_varexp = loaded_pca.get_exp_variance()
        loaded_pred = loaded_pca.predict(sample_pca_data[:5])

        assert abs(matrix_sum(loaded_scores) - matrix_sum(orig_scores)) < 1e-12
        assert abs(matrix_sum(loaded_loadings) - matrix_sum(orig_loadings)) < 1e-12
        assert abs(sum(loaded_varexp) - sum(orig_varexp)) < 1e-10
        assert abs(matrix_sum(loaded_pred) - matrix_sum(orig_pred)) < 1e-12

        # 5. Test PCA class-level save and load methods
        pca_direct_path = os.path.join(tmpdir, "pca_direct.lsm")
        pca.save(pca_direct_path)
        pca2 = PCA(scaling=1, npc=2)
        pca2.load(pca_direct_path)
        assert abs(matrix_sum(pca2.get_scores()) - matrix_sum(orig_scores)) < 1e-12


def test_pls_single_model_lsm(sample_pls_data):
    """Test saving and importing a single PLS model using .lsm format."""
    x, y = sample_pls_data
    pls = PLS(nlv=2, xscaling=1, yscaling=0)
    pls.fit(x, y)

    orig_pred_y, orig_scores = pls.predict(x[:5])

    with tempfile.TemporaryDirectory() as tmpdir:
        filepath = os.path.join(tmpdir, "test_pls.lsm")

        # 1. Save model
        save_model(filepath, pls, model_name="pls_chem")
        assert os.path.exists(filepath)

        # 2. Inspect metadata
        info = inspect_model_file(filepath)
        assert info["num_models"] == 1
        assert info["is_single"] is True
        assert info["primary_type"] == ModelType.PLS
        assert info["primary_name"] == "pls_chem"
        assert is_bundle(filepath) is False

        # 3. Import back into a new instance
        loaded_pls = import_model(filepath)
        assert isinstance(loaded_pls, PLS)
        assert loaded_pls.nlv == 2

        # 4. Verify numerical predictions
        loaded_pred_y, loaded_scores = loaded_pls.predict(x[:5])
        assert abs(matrix_sum(loaded_pred_y) - matrix_sum(orig_pred_y)) < 1e-12
        assert abs(matrix_sum(loaded_scores) - matrix_sum(orig_scores)) < 1e-12

        # 5. Test PLS class-level save and load methods
        pls_direct_path = os.path.join(tmpdir, "pls_direct.lsm")
        pls.save(pls_direct_path)
        pls2 = PLS(nlv=2, xscaling=1, yscaling=0)
        pls2.load(pls_direct_path)
        p_y2, p_s2 = pls2.predict(x[:5])
        assert abs(matrix_sum(p_y2) - matrix_sum(orig_pred_y)) < 1e-12


def test_multi_model_bundle_lsb(sample_pca_data, sample_pls_data):
    """Test packaging multiple models into a single .lsb bundle."""
    x_pls, y_pls = sample_pls_data

    # Train PCA model
    pca = PCA(scaling=1, npc=2)
    pca.fit(sample_pca_data)
    orig_pca_pred = pca.predict(sample_pca_data[:5])

    # Train PLS model
    pls = PLS(nlv=2, xscaling=1, yscaling=0)
    pls.fit(x_pls, y_pls)
    orig_pls_pred_y, _ = pls.predict(x_pls[:5])

    with tempfile.TemporaryDirectory() as tmpdir:
        bundle_path = os.path.join(tmpdir, "project_pipeline.lsb")

        # Create bundle and save
        bundle = ModelBundle(
            {
                "stage1_pca": pca,
                "stage2_pls": pls,
            }
        )
        assert len(bundle) == 2
        assert "stage1_pca" in bundle
        assert "stage2_pls" in bundle

        save_model_bundle(bundle_path, bundle)
        assert os.path.exists(bundle_path)

        # Inspect bundle metadata
        info = inspect_model_file(bundle_path)
        assert info["num_models"] == 2
        assert info["is_single"] is False
        assert is_bundle(bundle_path) is True

        # Import bundle
        loaded_bundle = import_model_bundle(bundle_path)
        assert len(loaded_bundle) == 2
        assert "stage1_pca" in loaded_bundle
        assert "stage2_pls" in loaded_bundle

        # Verify models inside the bundle work correctly
        loaded_pca = loaded_bundle["stage1_pca"]
        loaded_pls = loaded_bundle["stage2_pls"]

        l_pca_pred = loaded_pca.predict(sample_pca_data[:5])
        l_pls_pred_y, _ = loaded_pls.predict(x_pls[:5])

        assert abs(matrix_sum(l_pca_pred) - matrix_sum(orig_pca_pred)) < 1e-12
        assert abs(matrix_sum(l_pls_pred_y) - matrix_sum(orig_pls_pred_y)) < 1e-12


def test_import_auto_detection(sample_pca_data, sample_pls_data):
    """Test automatic inference/detection on both single .lsm and multi-model .lsb files."""
    pca = PCA(scaling=1, npc=2)
    pca.fit(sample_pca_data)

    x_pls, y_pls = sample_pls_data
    pls = PLS(nlv=2, xscaling=1, yscaling=0)
    pls.fit(x_pls, y_pls)

    with tempfile.TemporaryDirectory() as tmpdir:
        single_path = os.path.join(tmpdir, "single.lsm")
        bundle_path = os.path.join(tmpdir, "bundle.lsb")

        save_model(single_path, pca, model_name="my_pca")
        save_model_bundle(bundle_path, {"pca": pca, "pls": pls})

        # Test import_auto on single model
        auto_res1 = import_auto(single_path)
        assert isinstance(auto_res1, PCA)
        assert auto_res1.npc == 2

        # Test import_auto on bundle
        auto_res2 = import_auto(bundle_path)
        assert isinstance(auto_res2, ModelBundle)
        assert len(auto_res2) == 2
        assert "pca" in auto_res2
        assert "pls" in auto_res2


def test_error_handling(sample_pca_data):
    """Test error handling on invalid files, bad magics, and type mismatches."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # 1. Non-existent file
        non_existent = os.path.join(tmpdir, "missing.lsm")
        with pytest.raises(FileNotFoundError):
            inspect_model_file(non_existent)
        with pytest.raises(FileNotFoundError):
            import_model(non_existent)

        # 2. Corrupt file / Invalid magic
        bad_file = os.path.join(tmpdir, "bad.lsm")
        with open(bad_file, "wb") as f:
            f.write(b"CORRUPT_HEADER_BYTES_1234567890123456")

        with pytest.raises(IOError):
            inspect_model_file(bad_file)
        with pytest.raises(IOError):
            import_model(bad_file)

        # 3. Type mismatch
        pca = PCA(scaling=1, npc=2)
        pca.fit(sample_pca_data)
        pca_path = os.path.join(tmpdir, "pca.lsm")
        save_model(pca_path, pca)

        with pytest.raises(ValueError):
            import_model(pca_path, expected_type=ModelType.PLS)
