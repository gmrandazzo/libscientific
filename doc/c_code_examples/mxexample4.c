#include <stdio.h>
#include <scientific.h>


int main(void)
{
    int i, j;
    matrix *m, *m_T; // Definition of the pointer matrix variable
    NewMatrix(&m, 10, 15); // Create the matrix with 10 rows and 15 columns. Each value in the matrix is 0.
    NewMatrix(&m_T, m->col, m->row); // Create the transposed matrix with the flip of the columns and rows size

    for(i = 0; i < 10; i++){
        for(j = 0; j < 15; j++){
            m->data[i][j] = (float)i+j; // Fill the matrix values with these numbers
        }
    }
    MatrixTranspose(m, m_T);
    PrintMatrix(m); // Print to video the original matrix content
    puts("Transposed matrix");
    PrintMatrix(m_T); // Print to video the transposed matrix 
    // Free the memory spaces
    DelMatrix(&m);
    DelMatrix(&m_T);
}
