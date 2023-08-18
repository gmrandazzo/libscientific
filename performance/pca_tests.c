#include <stdio.h>
#include <scientific.h>
#include <sys/time.h>

int main(int argc, char **argv)
{
    matrix *m; // Definition of the input matrix 
    PCAMODEL *model; // Definition of the PCA model
    int i, j;
    int nobj = atoi(argv[1]);
    int nvars = atoi(argv[2]);
    struct timeval  tv1, tv2;
    NewMatrix(&m, nobj, nvars);

    srand_(time(0));
    for(size_t i = 0; i < nobj; i++){
        for(size_t j = 0; j < nvars; j++){
            m->data[i][j] = randDouble(0,20);
        }
    }
    
    gettimeofday(&tv1, NULL);
    NewPCAModel(&model); // Allocation of the PCA model
    PCA(m, 1, 5, model, NULL); // Calculation of the PCA on matrix m using unit variance scaling (1) and the extraction of 5 principal components 
    gettimeofday(&tv2, NULL);
    fprintf( stdout, "Total time (sec): %f\n",
        (double) (tv2.tv_usec - tv1.tv_usec) / 1000000 +
        (double) (tv2.tv_sec - tv1.tv_sec));
    DelPCAModel(&model);
    DelMatrix(&m);
}


