from libscientific.pls import * 
import random
from utils import (vector_sum,
                   raw_vector_sum,
                   matrix_sum)

def test_pls():
    random.seed(123456)
    x = [[random.random() for j in range(4)] for i in range(10)]
    y = [[random.random() for j in range(1)] for i in range(10)]
    xp = [[random.random() for j in range(4)] for i in range(10)]

    # fit directly from x and y
    model = PLS(nlv=2, xscaling=1, yscaling=0)
    model.fit(x, y)
    check_value_1 = 1.2140365662161021
    assert abs(raw_vector_sum(model.model.contents.xcolscaling)-check_value_1) <=1e-14

    # fit from libscientific.matrix.Matrix x and y
    xm = mx.Matrix(x)
    ym = mx.Matrix(y)
    model2 = PLS(nlv=2, xscaling=1, yscaling=0)
    model2.fit(xm, ym)
    tscores = model.get_tscores()
    check_value_2 = -4.440892098500626e-15
    assert abs(matrix_sum(tscores)-check_value_2) <= 1e-14

    uscores = model.get_uscores()
    check_value_3 = 1.8318679906315083e-15
    assert abs(matrix_sum(uscores)-check_value_3) <= 1e-14

    ploadings = model.get_ploadings()
    check_value_4 = 0.8962458837779853
    assert abs(matrix_sum(ploadings)-check_value_4) <= 1e-14

    exp_var = model.get_exp_variance()
    check_value_5 = 52.49759758877033
    assert abs(sum(exp_var)-check_value_5) <= 1e-14

    b = model.get_beta_coefficients(2)
    check_value_6 = 0.16159894366126876
    assert abs(sum(b)-check_value_6) <= 1e-14
    del b

    py, pscores = model.predict(xp)
    check_value_7 = 6.963935779767737
    check_value_8 = -2.55854976478211
    assert abs(matrix_sum(py)-check_value_7) <= 1e-14
    assert abs(matrix_sum(pscores)-check_value_8) <= 1e-14

    # Predict using beta coefficients
    pby = pls_beta_inference(xp,
                            x_averages=model.get_x_averages(),
                            x_scaling=model.get_x_column_scaling(),
                            y_averages=model.get_y_averages(),
                            y_scaling=model.get_y_column_scaling(),
                            beta_coeff=model.get_beta_coefficients(2))    
    for i, _ in enumerate(pby):
        assert  abs(pby[i][-1]-py[i][-1]) <=1e-14

    model.save('pls.sqlite3')
    model3 = PLS()
    model3.load("pls.sqlite3")
    tscores3 = model3.get_tscores()
    assert abs(matrix_sum(tscores3)-check_value_2) <= 1e-14
    del model
    del model2
    del model3
