/* teststrvector.c
*
* Copyright (C) <2016>  Giuseppe Marco Randazzo
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "vector.h"

void test2()
{
  char *string = "This;is;a;string";
  strvector *sv;
  initStrVector(&sv);
  SplitString(string, ";", sv);
  StrVectorAppendDouble(sv, 1.23456);
  PrintStrVector(sv);
  DelStrVector(&sv);
}

void test1()
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
}


int main(void)
{
 //test1();
 test2();
}
