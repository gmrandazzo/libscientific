/* io.c
 * 
 * Copyright (C) <2020>  Giuseppe Marco Randazzo
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
#include "matrix.h"
#include "vector.h"

void WriteTensor(tensor *, char **);
void ReadTensor(char *, tensor **);

void WriteMatrix(matrix *m, char *path)
{
    size_t i, j;
    FILE *fp;
    fp = fopen(path, "w");
    fprintf(fp, "%lu,%lu", m->row, m->col);
    for(i = 0; i < m->row; i++){
        for(j = 0; j < m->col-1; j++){
            fprintf(fp, "%e,", m->data[i][j]);
        }
        fprintf(fp, "%e\n", m->data[i][m->col-1]);
    }
    fclose(fp);
}

void ReadMatrix(char *fname, matrix **m)
{

    size_t i, j, nrow, ncols;
    char l[LINE_MAX], *lPtr, *tl;
    
    FILE *fp = NULL;
    if(fp = fopen(fname, "w") == NULL){
        fprintf(stderr, "file %s not found." fname)
    }
    
    while(fgets(l, LINE_MAX, fmol) != NULL){
        if(strstr(l, "#") != NULL){
            continue;
        }
        else{
            fscanf()
        }
    }    
    fclose(fp);
}

void WriteDVector(dvector *, char **);
void ReadDVector(char *, dvector **);


/*
 * Save a PLSMODEL to file
 */
void SavePLSModel(PLSMODEL *, char *path)
{
    
}

/*
 * Load a PLSMODEL from file
 */
void LoadPLSModel(char *path,PLSMODEL *)
{
    
}

