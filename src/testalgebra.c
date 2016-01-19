/* testalgebra.c
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
#include "algebra.h"
#include "matrix.h"
#include "vector.h"

void TestPolyFit()
{
  matrix *x; /* this matrix explain all the indipendent variables x for the polinomyal equation */
  dvector *y; /* this vector represent the dependent value associated for each row of the x matrix */
  dvector *b; /* this vector represent the polynomial coefficients. */
  
  NewMatrix(&x, 19, 10); /* 10 variables with 19 experiments for CCD */
  NewDVector(&y, 19); /* 19 dependent value */
  
  initDVector(&b);
  /*NewDVector(&b, 10);  10 coefficients that are equal to the number of variables */
  
  puts(" >>>>> Polynomial Fitting Test");
  puts("Answer: Coefficients \n b0: 58.78 \n b1: 1.064179 \n b2: 0.971129 \n b3: 1.898212 \n b4: -0.953331 \n b5: -0.686538");
  puts(" b6: -1.876163 \n b7: -1.236267 \n b8: -2.174703 \n b9: -2.706770");
  
  x->data[0][0] = 1.0000; 
  x->data[0][1] = -1.0000; 
  x->data[0][2] = -1.0000;  
  x->data[0][3] =-1.0000;  
  x->data[0][4] =1.0000; 
  x->data[0][5] = 1.0000;
  x->data[0][6] = 1.0000;
  x->data[0][7] = 1.0000;
  x->data[0][8] = 1.0000;
  x->data[0][9] = 1.0000;
  y->data[0] = 45.9;
  
  x->data[1][0] = 1.0000;
  x->data[1][1] = 1.0000;
  x->data[1][2] = -1.0000;
  x->data[1][3] = -1.0000;
  x->data[1][4] = 1.0000;
  x->data[1][5] = 1.0000;
  x->data[1][6] = 1.0000;
  x->data[1][7] = -1.0000;
  x->data[1][8] = -1.0000;
  x->data[1][9] = 1.0000;
  y->data[1] = 53.3;
  
  
  x->data[2][0] = 1.0000;
  x->data[2][1] = -1.0000;
  x->data[2][2] = 1.0000;
  x->data[2][3] = -1.0000;
  x->data[2][4] = 1.0000; 
  x->data[2][5] = 1.0000;
  x->data[2][6] = 1.0000;
  x->data[2][7] = -1.0000;
  x->data[2][8] = 1.0000;
  x->data[2][9] = -1.0000;
  y->data[2] = 57.5;
  
  x->data[3][0] = 1.0000;
  x->data[3][1] = 1.0000;
  x->data[3][2] = 1.0000;
  x->data[3][3] = -1.0000;
  x->data[3][4] = 1.0000;
  x->data[3][5] = 1.0000;
  x->data[3][6] = 1.0000;
  x->data[3][7] = 1.0000;
  x->data[3][8] = -1.0000;
  x->data[3][9] = -1.0000;
  y->data[3] = 58.8;
  
  x->data[4][0] = 1.0000;
  x->data[4][1] = -1.0000;
  x->data[4][2] = -1.0000;
  x->data[4][3] = 1.0000;
  x->data[4][4] = 1.0000;
  x->data[4][5] = 1.0000;
  x->data[4][6] = 1.0000;
  x->data[4][7] = 1.0000;
  x->data[4][8] = -1.0000;
  x->data[4][9] = -1.0000;
  y->data[4] = 60.6;
  
  x->data[5][0] = 1.0000;
  x->data[5][1] = 1.0000;
  x->data[5][2] = -1.0000;
  x->data[5][3] = 1.0000;
  x->data[5][4] = 1.0000;
  x->data[5][5] = 1.0000;
  x->data[5][6] = 1.0000;
  x->data[5][7] = -1.0000;
  x->data[5][8] = 1.0000; 
  x->data[5][9] = -1.0000;
  y->data[5] = 58;
  
  x->data[6][0] = 1.0000;
  x->data[6][1] = -1.0000;
  x->data[6][2] = 1.0000;
  x->data[6][3] = 1.0000;
  x->data[6][4] = 1.0000;
  x->data[6][5] = 1.0000;
  x->data[6][6] = 1.0000;
  x->data[6][7] = -1.0000;
  x->data[6][8] = -1.0000;
  x->data[6][9] = 1.0000;
  y->data[6] = 58.6;
  
  x->data[7][0] = 1.0000;
  x->data[7][1] = 1.0000;
  x->data[7][2] = 1.0000;
  x->data[7][3] = 1.0000;
  x->data[7][4] = 1.0000;
  x->data[7][5] = 1.0000;
  x->data[7][6] = 1.0000;
  x->data[7][7] = 1.0000;
  x->data[7][8] = 1.0000;
  x->data[7][9] = 1.0000;
  y->data[7] = 52.4;
  
  x->data[8][0] = 1.0000;
  x->data[8][1] = 0.0000;
  x->data[8][2] = 0.0000;
  x->data[8][3] = 0.0000;
  x->data[8][4] = 0.0000;
  x->data[8][5] = 0.0000;
  x->data[8][6] = 0.0000;
  x->data[8][7] = 0.0000;
  x->data[8][8] = 0.0000;
  x->data[8][9] = 0.0000;
  y->data[8] = 56.9;
  
  x->data[9][0] = 1.0000;
  x->data[9][1] =  0.0000;
  x->data[9][2] =  0.0000;
  x->data[9][3] =  2.0000;
  x->data[9][4] =  0.0000;
  x->data[9][5] =  0.0000;
  x->data[9][6] = 4.0000;
  x->data[9][7] =  0.0000;
  x->data[9][8] =  0.0000;
  x->data[9][9] =  0.0000;
  y->data[9] = 55.4;
  
  x->data[10][0] = 1.0000;
  x->data[10][1] = 0.0000;
  x->data[10][2] = 0.0000;
  x->data[10][3] = -2.0000;
  x->data[10][4] = 0.0000;
  x->data[10][5] = 0.0000;
  x->data[10][6] = 4.0000;
  x->data[10][7] = 0.0000;
  x->data[10][8] = -0.0000;
  x->data[10][9] = -0.0000;
  y->data[10] = 46.9;
  
  x->data[11][0] = 1.0000;
  x->data[11][1] = 0.0000;
  x->data[11][2] = 2.0000;
  x->data[11][3] = 0.0000;
  x->data[11][4] = 0.0000;
  x->data[11][5] = 4.0000;
  x->data[11][6] = 0.0000;
  x->data[11][7] = 0.0000;
  x->data[11][8] = 0.0000;
  x->data[11][9] = 0.0000;
  y->data[11] = 57.5;
  
  x->data[12][0] = 1.0000;
  x->data[12][1] =  0.0000;
  x->data[12][2] =  -2.0000;
  x->data[12][3] =  0.0000;
  x->data[12][4] =  0.0000;
  x->data[12][5] =  4.0000;
  x->data[12][6] =  0.0000;
  x->data[12][7] =  -0.0000;
  x->data[12][8] =  0.0000;
  x->data[12][9] =  -0.0000;
  y->data[12] =  55;
  
  x->data[13][0] = 1.0000;
  x->data[13][1] = 2.0000;
  x->data[13][2] = 0.0000;
  x->data[13][3] = 0.0000;
  x->data[13][4] = 4.0000;
  x->data[13][5] = 0.0000;
  x->data[13][6] = 0.0000;
  x->data[13][7] = 0.0000;
  x->data[13][8] = 0.0000;
  x->data[13][9] = 0.000;
  y->data[13] = 58.9;
  
  
  x->data[14][0] = 1.0000;
  x->data[14][1] = -2.0000; 
  x->data[14][2] = 0.0000;
  x->data[14][3] = 0.0000;
  x->data[14][4] = 4.0000;
  x->data[14][5] = 0.0000;
  x->data[14][6] = 0.0000;
  x->data[14][7] = -0.0000;
  x->data[14][8] = -0.0000;
  x->data[14][9] = 0.0000;
  y->data[14] = 50.3;
  
  x->data[15][0] = 1.0000;
  x->data[15][1] = 0.0000;
  x->data[15][2] = -3.0000; 
  x->data[15][3] = 2.0000;
  x->data[15][4] = 0.0000;
  x->data[15][5] = 9.0000;
  x->data[15][6] = 4.0000;
  x->data[15][7] = -0.0000;
  x->data[15][8] = 0.0000;
  x->data[15][9] = -6.0000;
  y->data[15] = 61.1;
  
  x->data[16][0] = 1.0000;
  x->data[16][1] = 0.0000;
  x->data[16][2] = -3.0000;
  x->data[16][3] = 2.0000;
  x->data[16][4] = 0.0000;
  x->data[16][5] = 9.0000;
  x->data[16][6] = 4.0000;
  x->data[16][7] = -0.0000;
  x->data[16][8] = 0.0000;
  x->data[16][9] = -6.0000;
  y->data[16] = 62.9;
  
  x->data[17][0] = 1.0000;
  x->data[17][1] = 0.6667;
  x->data[17][2] = 2.6000;
  x->data[17][3] = -1.4000;
  x->data[17][4] = 0.4444;
  x->data[17][5] = 6.7600;
  x->data[17][6] = 1.9600;
  x->data[17][7] = 1.7333;
  x->data[17][8] = -0.9333;
  x->data[17][9] = -3.6400;
  y->data[17] = 60;
  
  x->data[18][0] = 1.0000;
  x->data[18][1] = 0.6667;
  x->data[18][2] = 2.6000;
  x->data[18][3] = -1.4000;
  x->data[18][4] = 0.4444;
  x->data[18][5] = 6.7600;
  x->data[18][6] = 1.9600;
  x->data[18][7] = 1.7333;
  x->data[18][8] = -0.9333;
  x->data[18][9] = -3.6400;
  y->data[18] = 60.6;
       
  puts("Calculating...");
  OrdinaryLeastSquares(x, y, b);
  
  printf("Coefficients value\n");
  PrintDVector(b);
  
  DelDVector(&b);
  DelDVector(&y);
  DelMatrix(&x);

  
}

void testLESolv4()
{
  long int i;
  matrix *equation;
  dvector *c; /*coefficients results*/
  puts(">>>>>>>> Test4");
  puts("Problem: Solve the following system:\n 5z  = 3 \n w -x + z = -4 \n x - 3y = 2\n -w -y = -1\n  Answer (0.100000, 4.700000, 0.900000, 0.600000)");
  NewMatrix(&equation, 4, 5); /* 4 row(equation) and 4 column*/
  equation->data[0][0] = 0.f; equation->data[0][1] = 0.f; equation->data[0][2] = 0.f; equation->data[0][3] = 5; equation->data[0][4] = 3;
  equation->data[1][0] = 1; equation->data[1][1] = -1; equation->data[1][2] = 0.f; equation->data[1][3] = 1; equation->data[1][4] = -4;  
  equation->data[2][0] = 0.f; equation->data[2][1] = 1; equation->data[2][2] = -3; equation->data[2][3] = 0.f; equation->data[2][4] = 2;
  equation->data[3][0] = -1; equation->data[3][1] = 0.f; equation->data[3][2] = -1; equation->data[3][3] = 0.f; equation->data[3][4] = -1;
  
  initDVector(&c);
  puts("Compute...");
  SolveLSE(equation, &c);
  
  printf("The answer is: (");
  for(i = 0; i < c->size; i++){
    if(i < c->size-1)
      printf("%f, ", c->data[i]);
    else
      printf("%f).\n", c->data[i]);
      
  }
  
  DelMatrix(&equation);
  DelDVector(&c);
  
}


void testLESolv3()
{
  long int i;
  matrix *equation;
  dvector *c; /*coefficients results*/
  puts(">>>>>>>> Test3");
  puts("Problem: Solve the following system:\n 3y + 5z  = 3 \n 1w -x - 2y + z = -4 \n 2w + 4x - 3y - 7z = 2\n w +3y = 2\n  Answer (-4.4821, 3.1428, 2.1607, 0.6964)");
  
  NewMatrix(&equation, 4, 5); /* 4 row(equation) and 4 column*/
  equation->data[0][0] = 0; equation->data[0][1] = 0; equation->data[0][2] = 3; equation->data[0][3] = 5; equation->data[0][4] = 3;
  equation->data[1][0] = 1; equation->data[1][1] = -1; equation->data[1][2] = 2; equation->data[1][3] = 1; equation->data[1][4] = -4;  
  equation->data[2][0] = 2; equation->data[2][1] = 4; equation->data[2][2] = -3; equation->data[2][3] = -7; equation->data[2][4] = 2;
  equation->data[3][0] = 1; equation->data[3][1] = 0; equation->data[3][2] = 3; equation->data[3][3] = 0.f; equation->data[3][4] = 2;
  
  initDVector(&c);
  puts("Compute...");
  SolveLSE(equation, &c);
  
  printf("The answer is: (");
  for(i = 0; i < c->size; i++){
    if(i < c->size-1)
      printf("%f, ", c->data[i]);
    else
      printf("%f).\n", c->data[i]);
      
  }
  
  DelMatrix(&equation);
  DelDVector(&c);
  
}


void testLESolv2()
{
  long int i;
  matrix *equation;
  dvector *c; /*coefficients results*/
  puts(">>>>>>>> Test2");
  puts("Problem: Solve the following system:\n 2w - 1x + 3y + 5z  = 3 \n -x - 2y + z = -4 \n 2w + 4x - 3y - 7z = 2\n w +3y = 2\n  Answer (-7.947368, 15.526316, 3.315789, 4.894737)");
  NewMatrix(&equation, 4, 5); /* 4 row(equation) and 4 column*/
  equation->data[0][0] = 2; equation->data[0][1] = -1; equation->data[0][2] = 3; equation->data[0][3] = 5; equation->data[0][4] = 3;
  equation->data[1][0] = 0; equation->data[1][1] = -1; equation->data[1][2] = 2; equation->data[1][3] = 1; equation->data[1][4] = -4;  
  equation->data[2][0] = 2; equation->data[2][1] = 4; equation->data[2][2] = -3; equation->data[2][3] = -7; equation->data[2][4] = 2;
  equation->data[3][0] = 1; equation->data[3][1] = 0; equation->data[3][2] = 3; equation->data[3][3] = 0.f; equation->data[3][4] = 2;
  
  initDVector(&c);
  puts("Compute...");
  SolveLSE(equation, &c);
  
  printf("The answer is: (");
  for(i = 0; i < c->size; i++){
    if(i < c->size-1)
      printf("%f, ", c->data[i]);
    else
      printf("%f).\n", c->data[i]);
      
  }
  
  DelMatrix(&equation);
  DelDVector(&c);
  
}

void testLESolv1()
{
  long int i;
  matrix *equation;
  dvector *c; /*coefficients results*/
  puts(">>>>>>>> Test1");
  puts("Problem: Solve the following system:\n x + y + z  = 4 \n x - 2y - z = 1 \n 2x - y - 2z = -1\n Answer (2.000000, -1.000000, 3.000000)");
  NewMatrix(&equation, 3, 4); /* 3 row and 4 column*/
  equation->data[0][0] = 1; equation->data[0][1] = 1; equation->data[0][2] = 1; equation->data[0][3] = 4;
  equation->data[1][0] = 1; equation->data[1][1] = -2; equation->data[1][2] = -1; equation->data[1][3] = 1; 
  equation->data[2][0] = 2; equation->data[2][1] = -1; equation->data[2][2] = -2; equation->data[2][3] = -1; 
  
  initDVector(&c);
  puts("Compute...");
  SolveLSE(equation, &c);
  
  printf("The answer is: (");
  for(i = 0; i < c->size; i++){
    if(i < c->size-1)
      printf("%f, ", c->data[i]);
    else
      printf("%f).\n", c->data[i]);
      
  }
  
  DelMatrix(&equation);
  DelDVector(&c);
  
}

int main(void)
{
  testLESolv1();
  testLESolv2();
  testLESolv3();
  testLESolv4();
  TestPolyFit();
  return 0;
}
