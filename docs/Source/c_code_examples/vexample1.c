#include <stdio.h>
#include <scientific.h>


int main(void){
  int i;
  /* Define the variable vector, which in this case is a double vector.
   * If we would like to use an integer or an unsigned integer or a 
   * string vector instead this construct became
   * ivector *v; or uivector *v; or strvector *v;
   */

  dvector *v;
  NewDVector(&v, 5); // Allocate the memory space.

  // Fill the value inside the vector
  for(i = 0; i < 5; i++){
    v->data[i] = (float)i;
  }

  // Show the vector values to video
  PrintDVector(v);

  // Free the memory space
  DelDVector(&v);
}

