/* Multi-way PCA: 15 objects x 4 features measured in 3 conditions (layers). */
#include <stdio.h>
#include <math.h>
#include <scientific.h>

int main(void)
{
    tensor *X, *Xr;
    matrix *pscores;
    UPCAMODEL *model;
    size_t i, j, k;

    NewTensor(&X, 3);
    for(k = 0; k < 3; k++){
        NewTensorMatrix(X, k, 15, 4);   /* every layer: 15 objects x 4 features */
        for(i = 0; i < 15; i++)
            for(j = 0; j < 4; j++)
                X->m[k]->data[i][j] = 3.0 * (double)(j + 1)
                    + sin(0.9 * (double)(i + 1)) * (1.0 + 0.5 * (double)j) * (1.0 + 0.7 * (double)k);
    }

    NewUPCAModel(&model);
    UPCA(X, 1, 1, model, NULL);   /* 1 component, autoscaling */
    printf("explained variance: %.2f%%\n", model->varexp->data[0]);
    printf("loading matrix of component 0: %zu features x %zu layers\n",
           model->loadings->m[0]->row, model->loadings->m[0]->col);

    initMatrix(&pscores);
    UPCAScorePredictor(X, model, 1, pscores);
    printf("score of object 0: %.4f (model %.4f)\n", pscores->data[0][0], model->scores->data[0][0]);

    initTensor(&Xr);
    UPCAIndVarPredictor(model->scores, model->loadings, model->colaverage, model->colscaling, 1, Xr);
    printf("X[0][0][0] = %.4f, reconstructed %.4f\n", X->m[0]->data[0][0], Xr->m[0]->data[0][0]);

    DelTensor(&Xr);
    DelMatrix(&pscores);
    DelUPCAModel(&model);
    DelTensor(&X);
    return 0;
}
