#include <stdio.h>
#include <scientific.h>


int main(void){
  int i;

  /*
   * We define here the double vector. 
   * If we would like to utilize another vector
   * we can change dvector with the other three possibilities:
   * - uivector
   * - ivector
   * - strvector
   *
   */
  dvector *v;

  NewDVector(&v, 5); // We intialize the vector 

  // We add 5 numbers
  for(i = 0; i < 5; i++){
    v->data[i] = (float)i;
  }
 
  // We append a new number
  DVectorAppend(v, 123.4);
  
  // Print to video the result
  PrintDvector(v);

  // Free up the memory
  DelDVector(&v);
}

