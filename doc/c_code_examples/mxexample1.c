#include <stdio.h>
#include <scientific.h>


int main(void)
{
    int i, j;
    matrix *m; // Definition of the pointer matrix variable
    NewMatrix(&m, 10, 15); // Create the matrix with 10 rows and 15 columns. Each value in the matrix is 0.
    for(i = 0; i < 10; i++){
        for(j = 0; j < 15; j++){
            m->data[i][j] = (float)i+j; // Fill the matrix values with these numbers
        }
    }
    PrintMatrix(m); // Print to video the matrix content
    DelMatrix(&m); // Free the memory space
}
