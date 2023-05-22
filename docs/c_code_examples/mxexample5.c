#include <stdio.h>
#include <scientific.h>
#include <math.h>

int main(void)
{
    int i;
    matrix *m; // Definition of the matrix variable as a pointer
    matrix *m_inv; // Definition of the inverted matrix variable as a pointer 

    NewMatrix(&m, 10, 10); // Allocate the matrix to invert 
    MatrixInitRandomFloat(m, -3., 3.); // Random fill the matrix with values within a range -3 < x < 3
    PrintMatrix(m); // Print to video the matrix 

    initMatrix(&m_inv); // Initialize the matrix to invert
    MatrixInversion(m, m_inv); // Invert the matrix

    double det = fabs(MatrixDeterminant(m)); // Calculate the determinant

    printf("Determinant %.4f\n", det); // Print to video the matrix determinant
    PrintMatrix(m_inv); // Print to video the inverted matrix

    // Free the memory spaces
    DelMatrix(&m_inv);
    DelMatrix(&m);
}


