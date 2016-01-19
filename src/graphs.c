/* graphs.c
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

#include "graphs.h"

#include <math.h>
#include <string.h>

#include "memwrapper.h"
#include "numeric.h"
#include "matrix.h"
#include "scientificinfo.h"

void NewGraph(G **graph)
{
  (*graph) = xmalloc(sizeof(G));
  (*graph)->adjweight = NULL;
  (*graph)->adjmx = NULL;
  (*graph)->nodes = NULL;
  (*graph)->size_nodes = 0;
}

void DelGraph(G **graph)
{
  size_t i;
  if((*graph)->adjmx != NULL){
    for(i = 0; i < (*graph)->size_nodes; i++){
      if((*graph)->size_nodes-(i+1) > 0){
        xfree((*graph)->adjmx[i]);
        xfree((*graph)->adjweight[i]);
      }
    }
    xfree((*graph)->adjmx);
    xfree((*graph)->adjweight);
  }

  if((*graph)->nodes != NULL){
    for(i = 0; i < (*graph)->size_nodes; i++){
      xfree((*graph)->nodes[i].name);
    }
    xfree((*graph)->nodes);
  }
  xfree((*graph));
}

void GenerateAdjMX(matrix *dmx, double radius, G **graph)
{
  if(dmx->row == dmx->col){
    size_t i, j;
    (*graph)->size_nodes = dmx->row;
    printf("nodes %d\n", (*graph)->size_nodes);
    (*graph)->adjmx = xmalloc(sizeof(int*)*(*graph)->size_nodes);
    (*graph)->adjweight = xmalloc(sizeof(double*)*(*graph)->size_nodes);
    (*graph)->nodes = xmalloc(sizeof(node)*(*graph)->size_nodes);
    for(i = 0; i < dmx->row; i++){
      (*graph)->nodes[i].x = 0.;
      (*graph)->nodes[i].y = 0.;
      (*graph)->nodes[i].name = strdup("No Node Name");
      if((*graph)->size_nodes-(i+1) > 0){
        (*graph)->adjmx[i] = xmalloc(sizeof(int)*(*graph)->size_nodes-(i+1));
        (*graph)->adjweight[i] = xmalloc(sizeof(double)*(*graph)->size_nodes-(i+1));
      }
    }

    radius = radius + radius*0.01;
    for(i = 0; i < dmx->row; i++){
      for(j = i+1; j < dmx->col; j++){
        if(dmx->data[i][j] < radius){
          (*graph)->adjmx[i][j-(i+1)] = 1;
          (*graph)->adjweight[i][j-(i+1)] = dmx->data[i][j];
        }
        else{
          (*graph)->adjmx[i][j-(i+1)] = 0;
          (*graph)->adjweight[i][j-(i+1)] = 0.0;
        }
      }
    }
  }
  else{
    fprintf(stderr, "Error! The distance matrix must be squared!\n");
  }
}

void GenerateNodePositions(G *graph)
{
  
}

void FindMinimumPath(G *graph, uivector **minpath)
{

}

void PrintGraph(G *graph){
  size_t i, j;

  puts("Adjacence Matrix");
  /*for(i = 0; i < graph->size_nodes-1; i++){
    for(j = i+1; j < graph->size_nodes-1; j++){
      printf("%d\t", graph->adjmx[i][j-(i+1)]);
    }
    printf("%d\n", graph->adjmx[i][(graph->size_nodes-1)-(i+1)]);
  }*/
  for(i = 0; i < graph->size_nodes; i++){
    int inv = -1;
    for(j = 0; j < graph->size_nodes; j++){
      if(i == j){
        printf("%d\t", 0);
        inv  = 1;
      }
      else{
        if(inv == -1){
          if(j == 0){
            printf("%d\t", graph->adjmx[i+(j-(i))][j+(i-1)]);
            //printf("inv %d %d\n", i+(j-(i)),  j+(i-1));
          }
          else{
            printf("%d\t", graph->adjmx[i+(j-(i))][abs(j-(i-1))]);
            //printf("inv %d %d\n", i+(j-(i)),  (int)abs(j-(i-1)));
          }
          //printf("%d\t", graph->adjmx[j-(i+1)][i]);
        }
        else{
          //printf("%d %d\n", i, j-(i+1));
          printf("%d\t", graph->adjmx[i][j-(i+1)]);
        }
      }
    }
    printf("\n");
  }

  puts("Weights Matrix");
  for(i = 0; i < graph->size_nodes; i++){
    int inv = -1;
    for(j = 0; j < graph->size_nodes; j++){
      if(i == j){
        printf("%f\t", 0.0);
        inv  = 1;
      }
      else{
        if(inv == -1){
          if(j == 0){
            printf("%f\t", graph->adjweight[i+(j-(i))][j+(i-1)]);
          }
          else{
            printf("%f\t", graph->adjweight[i+(j-(i))][abs(j-(i-1))]);
          }
        }
        else{
          printf("%f\t", graph->adjweight[i][j-(i+1)]);
        }
      }
    }
    printf("\n");
  }

  /*
  for(i = 0; i < graph->size_nodes-1; i++){
    for(j = i+1; j < graph->size_nodes-1; j++){
      printf("%f\t", graph->adjweight[i][j-(i+1)]);
    }
    printf("%f\n", graph->adjweight[i][(graph->size_nodes-1)-(i+1)]);
  }*/


  if(graph->nodes != NULL){
    puts("Nodes");
    for(i = 0; i < graph->size_nodes; i++){
      printf("%s %f %f\n", graph->nodes[i].name, graph->nodes[i].x, graph->nodes[i].y);
    }
  }
}
