#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mmio.h"

void randPermColumns(int *col_ptrs, int * col_ids, double *col_vals, int n);
void randPermRows(int *col_ptrs, int * col_ids, double *col_vals, int n);
void ScaleOnePattern(	int *col_ptrs, int *col_ids, double *col_vals, int n, int m, int numIters );
void ScaleOne( int *col_ptrs, int *col_ids, double *col_vals, int n, int m, int numIters );
double computeScalingError( int *col_ptrs, int *col_ids, double *col_vals, int n, int m );
void ReadMatrixMarket( char *fname, double **A, int **csc_irn, int **csc_jcn,
			 int *_m, int *_n, int *_ne, int readVals );
