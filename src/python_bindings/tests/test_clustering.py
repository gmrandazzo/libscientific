import random

from utils import matrix_sum, raw_vector_sum

from libscientific.clustering import *
from libscientific.matrix import *
from libscientific.vector import *


def test_clustering():
    random.seed(123456)
    x = [[random.random() for j in range(2)] for i in range(20)]

    # sample from list with MDC
    selobj_ids = most_descriptive_compound(x, 5)
    assert sum(selobj_ids) == 61

    # sample from matrix with MaxDissimilaritySelection
    selobj_ids = max_dissimilarity_selection(x, 5)
    assert sum(selobj_ids) == 41

    lbls = k_means_plus_plus(x, 2)
    assert len(lbls) == 20

    # Test passing Matrix objects directly
    xm = Matrix(x)
    sel_mdc_m = most_descriptive_compound(xm, 5)
    assert sel_mdc_m == selobj_ids or len(sel_mdc_m) == 5

    sel_maxdis_m = max_dissimilarity_selection(xm, 5)
    assert len(sel_maxdis_m) == 5

    lbls_m = k_means_plus_plus(xm, 3)
    assert len(lbls_m) == 20

    # Test fast MDC
    sel_mdc_fast = most_descriptive_compound_fast(xm, 5, metric=0, nthreads=2)
    assert len(sel_mdc_fast) == 5

    # Test KMeans++ centers
    centers = k_means_pp_centers(xm, 3, nthreads=2)
    assert len(centers) == 3
