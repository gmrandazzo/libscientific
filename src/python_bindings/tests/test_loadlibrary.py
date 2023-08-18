
from libscientific.loadlibrary import * 


def test_loadlibrary():
    a = load_libscientific_library()
    assert a != None
    assert hasattr(a, "PCA")
    assert hasattr(a, "PLS")
    assert hasattr(a, "CPCA")
    assert hasattr(a, "LDA")
    assert hasattr(a, "MLR")
