#include <stdio.h>
#include <scientific.h>
#include <math.h>

int main(void)
{
    int i;
    matrix *m; // Definition of the matrix variable as a pointer
    dvector *eval; // Definition the variable to store the eigenvalues
    matrix *evect; // Definition of the variable were to store the eigenvectors

    NewMatrix(&m, 10, 10); // Allocate the matrix to invert 
    MatrixInitRandomFloat(m, -3., 3.); // Random fill the matrix with values within a range -3 < x < 3
    PrintMatrix(m); // Print to video the matrix 

    // Initialize the variables
    initDVector(&eval);
    initMatrix(&evect);

    EVectEval(m, eval, evect); // Calculate the eigenvectors and associated eigenvalues

    PrintDVector(eval); // Print to video the eigenvalues
    PrintMatrix(evect); // Print to video the eigenvectors. Each column correspond to an eingenvalue

    // Free the memory spaces
    DelDVector(&eval);
    DelMatrix(&evect);
    DelMatrix(&m);
}


