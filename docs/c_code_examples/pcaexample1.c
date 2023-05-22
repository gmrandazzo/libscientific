#include <stdio.h>
#include <scientific.h>

int main(void)
{
    matrix *m; // Definition of the input matrix 
    PCAMODEL *model; // Definition of the PCA model
    int i, j;
    int nobj = 20;
    int nvars = 8;
    NewMatrix(&m, nobj, nvars);

    // Fill with random values the matrix m
    srand(nobj);
    for(size_t i = 0; i < nobj; i++){
        for(size_t j = 0; j < nvars; j++){
            m->data[i][j] = randDouble(0,20);
        }
    }


    NewPCAModel(&model); // Allocation of the PCA model
    PCA(m, 1, 5, model, NULL); // Calculation of the PCA on matrix m using unit variance scaling (1) and the extraction of 5 principal components 

    PrintPCA(model); // Print to video the PCA results

    /* Of course you can print in a separate way the different results contained in the model variable
     * model->scores is the matrix of scores
     * model->loadings is the matrix of loadings
     * model->colavg is the column average obtained from the input matrix
     * model->scaling is the scaling factor obtained from the input matrix
     */

    // Free the memory spaces
    DelPCAModel(&model);
    DelMatrix(&m);
}

