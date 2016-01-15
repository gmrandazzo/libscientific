#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "vector.h"

/*sort a random vector*/
void test5()
{
  size_t i;
  dvector *v;
  srand(time(0));
  NewDVector(&v, 10);
  puts("Creating vector..");
  for(i = 0; i < v->size; i++){
    setDVectorValue(v, i, rand() % 100);
  }
  
  puts("Vector before sorting...");
  PrintDVector(v);
  
  DVectorSort(v);
  puts("Vector after sorting...");
  PrintDVector(v);
  DelDVector(&v);
}

/* Extend two vector into one vector*/

void test4_1(dvector **v1)
{
  size_t i;
  srand(time(0));
  NewDVector(v1, 100);
  printf("Creating v1\n");
  for(i = 0; i < (*v1)->size; i++){
    (*v1)->data[i] = rand() % 100;
  }
}

void test4()
{
   int i;
  dvector *v1;
  test4_1(&v1);
  
  printf("Test 4\n");
  for(i = 0; i < v1->size; i++){
   printf("%f\n", v1->data[i]);
  }
  DelDVector(&v1);
}

/* Extend two vector into one vector*/
void test3()
{
   int i;
  dvector *v1;
  dvector *v2;
  dvector *v1v2;
  
  printf("Test 3\n");
  
  initDVector(&v1);
  
  printf("Creating v1\n");
  for(i = 1; i < 100; i++){
    double a = i*i;
    double val = i/(a);
    DVectorAppend(&v1, val);
  }
  

    
  printf("Creating 21\n");
  NewDVector(&v2, 100);

  for(i = 0; i < v2->size; i++){
    v2->data[i] = i;
  }
 
  printf("Appending v1 to v2\n");
  v1v2 = DVectorExtend(v1,v2);
  
  printf("Final output\n");
  for(i = 0; i < v1v2->size; i++){
   printf("%f\n", v1v2->data[i]);
  }
  
  DelDVector(&v1v2);
  DelDVector(&v2);
  DelDVector(&v1);
}

/* Initialize the dvector by using initDVector function and then append values with DVectorAppend function*/
void test2()
{
  int i;
  dvector *v;
  
  printf("Test 2\n");
  
  initDVector(&v);
  
  printf("Appending 100 value\n");
  for(i = 1; i < 100; i++){
    double a = i*i;
    double val = i/(a);
    
    DVectorAppend(&v, val);
  }
  
  printf("Final output\n");
  for(i = 0; i < v->size; i++){
   printf("%f\n", v->data[i]);
  }
  
  DelDVector(&v);

}



/* Allocate the vector by using the NewDVector function */
void test1()
{
  int i;
  dvector *v;
  
  printf("Test 1\n");
  
  NewDVector(&v, 100);
  
  for(i = 0; i < 100; i++){
    v->data[i] = i;
  }
  
  printf("Appending 123.4 to vector\n");
  DVectorAppend(&v, 123.4);
  
  printf("Final output\n");
  for(i = 0; i < v->size; i++){
   printf("%f\n", v->data[i]);
  }
  
  DelDVector(&v);

}

int main(void)
{  
  test1();
  test2();
  test3();
  test4();
  test5();
  return 0;
}
