/* Multi-way PLS: predict one response from 2 layers of 4 features. */
#include <stdio.h>
#include <scientific.h>

int main(void)
{
    tensor *X, *Y, *ypred;
    matrix *tscores;
    UPLSMODEL *model;
    size_t i, j, k;

    NewTensor(&X, 2);
    NewTensor(&Y, 1);             /* a single Y matrix is a tensor with one layer */
    NewTensorMatrix(Y, 0, 20, 1);
    for(k = 0; k < 2; k++)
        NewTensorMatrix(X, k, 20, 4);
    for(i = 0; i < 20; i++){
        for(k = 0; k < 2; k++)
            for(j = 0; j < 4; j++)
                X->m[k]->data[i][j] = (double)((i * (j + 3) + k * 7) % 11) / 3.0;
        Y->m[0]->data[i][0] = X->m[0]->data[i][0] + 2.0 * X->m[0]->data[i][1] + X->m[1]->data[i][2];
    }

    NewUPLSModel(&model);
    UPLS(X, Y, 3, 1, 0, model, NULL);   /* 3 latent variables, autoscale X, center Y */

    initMatrix(&tscores);
    initTensor(&ypred);
    UPLSScorePredictor(X, model, 3, tscores);
    UPLSYPredictor(tscores, model, 3, ypred);
    for(i = 0; i < 3; i++)
        printf("object %zu: y %.3f, predicted %.3f\n", i, Y->m[0]->data[i][0], ypred->m[0]->data[i][0]);

    DelTensor(&ypred);
    DelMatrix(&tscores);
    DelUPLSModel(&model);
    DelTensor(&Y);
    DelTensor(&X);
    return 0;
}
