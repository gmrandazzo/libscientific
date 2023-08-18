#include <stdio.h>
#include <scientific.h>
#include <sys/time.h>

int main(int argc, char **argv)
{
    tensor *t; // Definition of the input matrix 
    CPCAMODEL *model; // Definition of the PCA model
    int i, j, k;
    int nobj = atoi(argv[1]);
    int nvars = atoi(argv[2]);
    int norder = atoi(argv[3]);
    struct timeval  tv1, tv2;
    NewTensor(&t, norder);
    

    for(k = 0; k < norder; k++){
        NewTensorMatrix(t, k, nobj, nvars);
        srand_(time(0));
        for(i = 0; i < nobj; i++){
            for(j = 0; j < nvars; j++){
                t->m[k]->data[i][j] = randDouble(0,100);
            }
        }
    }
    
    gettimeofday(&tv1, NULL);
    NewCPCAModel(&model);
    CPCA(t, 1, 5, model);
    gettimeofday(&tv2, NULL);
    fprintf( stdout, "Total time (sec): %f\n",
        (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
        (double) (tv2.tv_sec - tv1.tv_sec));
    DelCPCAModel(&model);
    DelTensor(&t);
}


