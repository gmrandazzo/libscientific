#include <stdio.h>
#include <scientific.h>


int main(void){
  int i;
  dvector *v;

  NewDVector(&v, 5);

  for(i = 0; i < 5; i++){
    v->data[i] = (float)i;
  }

  DVectorAppend(&v, 123.4);

  PrintDvector(v);

  DelDVector(&v);
}

