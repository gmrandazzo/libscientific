from libscientific.pca import * 
import random

def matrix_sum(m_input):
    s = 0
    for row in m_input:
        s += sum(row)
    return s

def test_pca():
    random.seed(123456)
    a = [[random.random() for j in range(4)] for i in range(10)]
    m = mx.Matrix(a)
    # computing pca from a and pca from m
    model1 = PCA(1, 2)
    model1.fit(a)
    model2 = PCA(1, 2)
    model2.fit(m)

    scores1 = model1.get_scores()
    scores2 = model2.get_scores()
    check_value_1 = 2.3314683517128287e-15
    assert abs(matrix_sum(scores1)- check_value_1) <= 1e-14
    assert abs(matrix_sum(scores2)- check_value_1) <= 1e-14
    
    loadings1 = model1.get_loadings()
    loadings2 = model1.get_loadings()
    check_value_2 = 0.8473526872707682
    assert abs(matrix_sum(loadings1)- check_value_2) <= 1e-14
    assert abs(matrix_sum(loadings2)- check_value_2) <= 1e-14

    exp_var_1 = model1.get_exp_variance()
    exp_var_2 = model2.get_exp_variance()
    check_value_3 = 68.22757533439007
    assert abs(sum(exp_var_1)- check_value_3) <= 1e-14
    assert abs(sum(exp_var_2)- check_value_3) <= 1e-14

    print("Reconstruct the original PCA matrix using the PCA Model")
    ra1 = model1.reconstruct_original_matrix()
    ra2 = model2.reconstruct_original_matrix()
    check_value_4 = 18.196996048027845
    assert abs(matrix_sum(ra1)- check_value_4) <= 1e-14
    assert abs(matrix_sum(ra2)- check_value_4) <= 1e-14

    pred_scores1 = model1.predict(a)
    pred_scores2 = model2.predict(m)
    assert abs(matrix_sum(pred_scores1)- check_value_1) <= 1e-14
    assert abs(matrix_sum(pred_scores2)- check_value_1) <= 1e-14
    del model1
    del model2
