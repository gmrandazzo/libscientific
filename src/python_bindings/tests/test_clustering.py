from libscientific.clustering import * 
from libscientific.matrix import *
from libscientific.vector import *
from utils import(matrix_sum,
                  raw_vector_sum)
import random

def test_clustering():
    random.seed(123456)
    x = [[random.random() for j in range(2)] for i in range(20)]

    # sample from matrix with MDC
    selobj_ids = most_descriptive_compound(x, 5)
    assert sum(selobj_ids) == 61

    # sample from matrix with MaxDissimilaritySelection
    selobj_ids = max_dissimilarity_selection(x, 5)
    assert sum(selobj_ids) == 41

    lbls = k_means_plus_plus(x, 2)
    #assert sum(lbls) == 11
