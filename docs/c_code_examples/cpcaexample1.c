#include <stdio.h>
#include <scientific.h>

int main(void)
{
    tensor *t; // Definition of the input tensor
    CPCAMODEL *model; // Definition of the CPCA model
    int i, j, k;
    int nblocks = 4;
    int nobj = 20;
    int nvars[4] = {8, 10, 5, 7}; //every block have different variables
    NewTensor(&t, nblocks);

    // Fill with random values the matrix m
    srand(nobj);
    for(k = 0; k < nblocks; k++){
        NewTensorMatrix(t, k, nobj, nvars[k]);
        for(i = 0; i < nobj; i++){
            for(j = 0; j < nvars[k]; j++){
                t->m[k]->data[i][j] = randDouble(0,20);
            }
        }
    }


    NewCPCAModel(&model); // Allocation of the CPCA model
    CPCA(t, 1, 5, model); // Calculation of the CPCA on matrix m using unit variance scaling (1) and the extraction of 5 super principal components 

    PrintCPCA(model); // Print to video the CPCA results

    /* Of course you can print in a separate way the different results contained in the model variable
     * model->super_scores is the matrix of super scores
     * model->super_weights is the matrix of super weights
     * model->block_scores is the tensor of scores for each block
     * model->block_loadings is the tensor of loadings for each block
     */

    // Free the memory spaces
    DelCPCAModel(&model);
    DelTensor(&t);
}

