#include <stdio.h>
#include <scientific.h>
#include <math.h>

int main(void)
{
    int i;
    matrix *m; // Definition of the matrix variable as a pointer
    matrix *U; // Definition of the complex unitary matrix
    matrix *S; // Definition of the rectangular diagnonal matrix with non-negative real numbers on diagonal
    matrix Vt; // Definition of conjugate transpose of a complex unitary matrix

    NewMatrix(&m, 10, 10); // Allocate the matrix to invert 
    MatrixInitRandomFloat(m, -3., 3.); // Random fill the matrix with values within a range -3 < x < 3
    PrintMatrix(m); // Print to video the matrix 

    // Initialize the variables
    initMatrix(&U);
    initMatrix(&S);
    initMatrix(&Vt);
    SVD(m, U, S, Vt); // Internal method
    // SVDlapack(m, U, S, Vt); lapack method using dgesdd
    // Print to video the results of the factorization
    PrintMatrix(U);
    PrintMatrix(S);
    PrintMatrix(Vt);

    // Free the memory spaces
    DelMatrix(&U);
    DelMatrix(&S);
    DelMatrix(&Vt);
    DelMatrix(&m);
}


