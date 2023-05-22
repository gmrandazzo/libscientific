#include <stdio.h>
#include <scientific.h>

int main(void)
{
    int i;
    matrix *mx; // Definition of the matrix variable as a pointer
    dvector *row; // Definition of the row variable as a pointer 

    NewDVector(&row, 15);
    for(i = 0; i < row->size; i++){
        row->data[i] = (double)i; // Fill one time the row vector 
    }

    initMatrix(&mx); // Initialize the empty matrix with rows and columns equal to 0

    for(i = 0; i < 5; i++){
        MatrixAppendRow(mx, row); // Append 5 times the row to the matrix mx
    }

    PrintMatrix(mx); // Print to video the matrix 

    DelDVector(&row); // Free the memory space for the row vector
    DelMatrix(&mx); // Free the memory space for the matrix
}


