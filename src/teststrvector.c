#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "vector.h"

int main(void)
{
  size_t i, j;
  strvector *s;
  char *letters = "QWERTYUIOPASDFGHJKLZXCVBNMqwertyuiopasdfghjklzxcvbnm123456789"; 
  int num_of_letters = strlen(letters);
  char *buffer = NULL;
  buffer = malloc(sizeof(char)*101);
  
  NewStrVector(&s, 100);
  srand(time(0));

  for(i = 0; i < s->size; i++){
    memset(buffer, 0, 101);
    for(j = 0; j < 100; j++){
      int n = rand() % num_of_letters; 
      memcpy(&buffer[j], &letters[n], 1);
    }
    buffer[j] = 0; /* 0 byte at the end of buffer */
    
    setStr(s, i, buffer);
  }
  
  for(i = 0; i < s->size; i++){
    printf("%d\t%s\n", (int)strlen(s->data[i]), getStr(s, i));
  }
  
  free(buffer);
  DelStrVector(&s);
  return 0;

}