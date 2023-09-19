from libscientific.cpca import * 
import random
from utils import(matrix_sum,
                  tensor_sum)

def test_cpca():       
    random.seed(123456)
    # create a tensor of 4 layers with 10 objects and 4 columns
    a = [[[random.random() for j in range(4)] for i in range(10)] for k in range(4)]
    t = tns.Tensor(a)
    # compute cpca from list of list of list a
    cpca1 = CPCA(1, 2)
    cpca1.fit(a)
    # compute cpca from libscientific t
    cpca2 = CPCA(1, 2)
    cpca2.fit(t)
    #compare results
    sscores1 = cpca1.get_super_scores()
    sscores2 = cpca2.get_super_scores()
    check_value_1 = -6.661338147750939e-16
    assert abs(matrix_sum(sscores1)- check_value_1) <= 1e-14
    assert abs(matrix_sum(sscores2)- check_value_1) <= 1e-14
    
    sweights1 = cpca1.get_super_weights()
    sweights2 = cpca1.get_super_weights()
    check_value_2 = 3.8134031061053206
    assert abs(matrix_sum(sweights1)- check_value_2) <= 1e-14
    assert abs(matrix_sum(sweights2)- check_value_2) <= 1e-14

    block_scores1 = cpca1.get_block_scores()
    block_scores2 = cpca1.get_block_scores()
    check_value_3 = -3.941291737419306e-15
    assert abs(tensor_sum(block_scores1)- check_value_3) <= 1e-14
    assert abs(tensor_sum(block_scores2)- check_value_3) <= 1e-14

    block_loadings1 = cpca1.get_block_loadings()
    block_loadings2 = cpca1.get_block_loadings()
    check_value_4 = -0.9901586826305672
    assert abs(tensor_sum(block_loadings1)- check_value_4) <= 1e-14
    assert abs(tensor_sum(block_loadings2)- check_value_4) <= 1e-14

    exp_var1 = cpca1.get_total_exp_variance()
    exp_var2 = cpca2.get_total_exp_variance()
    check_value_5 = 48.07849638312624
    assert abs(sum(exp_var1)- check_value_5) <= 1e-14
    assert abs(sum(exp_var2)- check_value_5) <= 1e-14
    
    pred_sscores1, pred_bsscores1 = cpca1.predict(a)
    pred_sscores2, pred_bsscores2 = cpca2.predict(a)
    check_value_6 = -4.440892098500626e-16
    check_value_7 = -1.0547118733938987e-15
    assert abs(matrix_sum(pred_sscores1)- check_value_6) <= 1e-14
    assert abs(tensor_sum(pred_bsscores1)- check_value_7) <= 1e-14
    assert abs(matrix_sum(pred_sscores2)- check_value_6) <= 1e-14
    assert abs(tensor_sum(pred_bsscores2)- check_value_7) <= 1e-14

    cpca1.save('cpca.sqlite3')
    cpca3 = CPCA()
    cpca3.load("cpca.sqlite3")
    sscores3 = cpca1.get_super_scores()
    assert abs(matrix_sum(sscores3)- check_value_1) <= 1e-14
    del cpca1
    del cpca2
    del cpca3
