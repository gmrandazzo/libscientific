#ifndef GRAPHS_H
#define GRAPHS_H
#include "matrix.h"
#include "scientificinfo.h"


typedef struct{
  double x;
  double y;
  char *name;
} node;

typedef struct{
  int size_nodes;
  int **adjmx;
  double **adjweight;
  node *nodes;
} G;

/* Initialize a new graph in memory*/
void NewGraph(G **graph);

/* Delete a graph from memory */
void DelGraph(G **graph);

/*
 * Generate the adjacence matrix from a distance matrix.
 * If radius is specified and not NULL then the edges were generated if the
 * distance is <= of the radius.
 */
void GenerateAdjMX(matrix *dmx, double radius, G **graph);

/*
 * Convert the adjacence matrix to square matrix.

void ConvertAdjMXtoSquareMX(G *graph, matrix **adj);
 */
 
/*
 * Generate the node position according to the FRUCHTERMAN-REINGOLD Algorithm.
 *
 * Reference:
 * Graph drawing by force-directed placement
 * Thomas M. J. Fruchterman† andEdward M. Reingold
 * Software: Practice and Experience
 * Volume 21, Issue 11, pages 1129–1164, November 1991
 * DOI: 10.1002/spe.4380211102
 */
void GenerateNodePositions(G *graph);

void FindMinimumPath(G *graph, uivector **minpath);

void PrintGraph(G *graph);


#endif
