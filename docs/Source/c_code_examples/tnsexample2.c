#include <stdio.h>
#include <scientific.h>


int main(void)
{
    size_t i, j, k;
    tensor *t; // Definition of the pointer tensor variable
    matrix *m;
    initTensor(&t); // Tensor allocation
    NewMatrix(&m, 10, 7); // Create of a matrix with 10 rows and 7 columns
    TensorAppendMatrix(t, m); // Append this matrix to the tensor t
    DelMatrix(&m); // Delete the matrix
    
    NewMatrix(&m, 10, 10); // Create of a second matrix with 10 rows and 10 columns
    TensorAppendMatrix(t, m); // Append this new matrix to the tensor t
    DelMatrix(&m);
    
    for(k = 0; k < t->order; k++){
        for(i = 0; i < t->m[k]->row; i++){
            for(j = 0; j < t->m[k]->col; j++){
                t->m[k]->data[i][j] = (float)i+j; // Fill the tensor values with these numbers
            }
        }
    }
    PrintTensor(t); // Print to video the tensor content
    DelTensor(&t); // Free the memory space
}
