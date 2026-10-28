from libscientific.loadlibrary import *


def test_loadlibrary():
    a = load_libscientific_library()
    assert a is not None
    assert hasattr(a, "PCA")
    assert hasattr(a, "PLS")
    assert hasattr(a, "CPCA")
    assert hasattr(a, "LDA")
    assert hasattr(a, "MLR")
    assert hasattr(a, "SaveModel")
    assert hasattr(a, "ImportModel")
    assert hasattr(a, "SaveModelBundle")
    assert hasattr(a, "ImportModelBundle")
    assert hasattr(a, "InspectModelFile")
