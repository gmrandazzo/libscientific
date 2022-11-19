#include <stdio.h>
#include <scientific.h>


int main(void)
{
    int i, j, k;
    tensor *t; // Definition of the pointer tensor variable
    NewTensor(&t, 3); // Create the tensor with 3 blocks;
    for(k = 0; k < 3; k++){
        NewTensorMatrix(t, k, 10, 15);
        for(i = 0; i < 10; i++){
            for(j = 0; j < 15; j++){
                t->m[k]->data[i][j] = (float)i+j; // Fill the tensor values with these numbers
            }
        }
    }
    PrintTensor(t); // Print to video the tensor content
    DelTensor(&t); // Free the memory space
}
