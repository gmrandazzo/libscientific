/* Multiple Linear Regression: fit y = 1 + 2 x1 - 0.5 x2 and predict new objects. */
#include <stdio.h>
#include <scientific.h>

int main(void)
{
    matrix *x, *y, *xnew, *ypred;
    MLRMODEL *model;
    size_t i;

    NewMatrix(&x, 6, 2);
    NewMatrix(&y, 6, 1);
    for(i = 0; i < 6; i++){
        x->data[i][0] = (double)i;
        x->data[i][1] = (double)((i * 7) % 4);
        y->data[i][0] = 1.0 + 2.0 * x->data[i][0] - 0.5 * x->data[i][1];
    }

    NewMLRModel(&model);
    MLR(x, y, model, NULL);

    /* model->b: (n_features + 1) x n_y; row 0 is the intercept */
    printf("intercept %.3f  b1 %.3f  b2 %.3f  R2 %.4f\n",
           model->b->data[0][0], model->b->data[1][0], model->b->data[2][0],
           model->r2y_model->data[0]);

    NewMatrix(&xnew, 2, 2);
    xnew->data[0][0] = 10.0; xnew->data[0][1] = 2.0;
    xnew->data[1][0] = -1.0; xnew->data[1][1] = 0.0;
    initMatrix(&ypred);
    MLRPredictY(xnew, NULL, model, ypred, NULL, NULL, NULL);
    for(i = 0; i < ypred->row; i++)
        printf("predicted y[%zu] = %.3f\n", i, ypred->data[i][0]);

    DelMatrix(&ypred);
    DelMatrix(&xnew);
    DelMLRModel(&model);
    DelMatrix(&y);
    DelMatrix(&x);
    return 0;
}
