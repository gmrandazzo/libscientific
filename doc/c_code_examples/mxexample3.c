#include <stdio.h>
#include <scientific.h>

int main(void)
{
    int i, j;
    matrix *mx; // Definition of the matrix variable as a pointer
    dvector *cvect; // Definition of the column vector as a pointer 
    dvector *result; // Definition of the result vector between the matrix and the column vector

    NewDVector(&cvect, 15);
    for(i = 0; i < cvect->size; i++){
        cvect->data[i] = (double)i; // Fill one time the row vector 
    }

    NewMatrix(&mx, 23, 15); // Initialize a matrix with 23 rows and 15 columns

    for(i = 0; i < mx->row; i++){
        for(j = 0; j < mx->col; j++){
            mx->data[i][j] = (double)i+j;
        }
    }
    NewDVector(&result, 23);

    MatrixDVectorDotProduct(mx, cvect, result);
    /*
     * or MT_MatrixDVectorDotProduct if you want to run the multitask operation.
     * This function is usefull for large matrix.
     */


    PrintDVector(result); // Print to video the result
    // Free the memory spaces
    DelDVector(&result); 
    DelDVector(&cvect);
    DelMatrix(&mx);
}


