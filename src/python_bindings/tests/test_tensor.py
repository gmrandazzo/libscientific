from libscientific.tensor import *
import numpy as np
from utils import(matrix_sum,
                  raw_vector_sum)
import random


def test_tensor():
    a = [[[1,2,3],[1,2,3]],[[4,5,6],[4,5,6]]]
    t = new_tensor(a)
    #print_tensor(t)
    lst = tensor_tolist(t)
    print(lst)
    add_tensor_matrix(t, 5, 6)
    #print_tensor(t)
    b = init_tensor()
    tensor_copy(t, b)
    #print_tensor(b)
    del_tensor(t)
    del_tensor(b)

    t = init_tensor()
    add_tensor_matrix(t, 3, 2)
    #print_tensor(t)
    del_tensor(t)

    t = Tensor(a)
    #print(type(t))
    #print(t.order())
    #print(t.nrow(0))
    #print(t.ncol(0))
    #print(t.tolist())
    a[0][0][0] = -9876
    t.fromlist(a)
    #t.debug()

    a = np.array(a)
    a[0][1][1] = 765
    t.fromnumpy(a)
    #t.debug()
