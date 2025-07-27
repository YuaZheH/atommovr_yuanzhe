#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>

#include "mmio.h"
#include "portdefs.h"
#include "matrixUtils.h"


static double u_wseconds(void) 
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return (double) tp.tv_sec + (double) tp.tv_usec / 1000000.0;
}

int findMinDegs(int *col_ids, int *col_ptrs, int n, int m)
{
  int *rdegs = (int *) malloc(sizeof(int) * m);
  int i;
  int minDeg = m+n;
  for (i = 0; i < m; i++)
    rdegs[i] = 0;

  for (i=0; i < n; i++)
  {
    int j;
    if(col_ptrs[i+1] - col_ptrs[i] > minDeg)
      minDeg = col_ptrs[i+1] - col_ptrs[i];
    for (j = col_ptrs[i]; j < col_ptrs[i+1]; j++)
      rdegs[col_ids[j]] ++;

  }
  for (i=0; i < m; i++)
    if(rdegs[i] < minDeg) minDeg = rdegs[i];

  free(rdegs);
  return minDeg;
}




int main(int argc, char *argv[])
{

	int j, nz;

	int m, n, minmn, mdeg;
	double *col_vals;
	double t0, t1;

	int* col_ptrs, * col_ids;

	int permute = 0;
	int maxcmatching = 0;
	if ( argc != 3 ){
	  printf("Usage: %s matrix_name  permute_columns\n", argv[0]);
		printf("\tpermutation (1,2,3, 4): no, col, row, row and col\n");
		exit(12);
	}

	ReadMatrixMarket(argv[1], &col_vals, &col_ids, &col_ptrs, &m,	&n, &nz, 0);

	for(j = 0; j<nz; j++) 
	{
	  col_vals[j] = col_vals[j] >=0 ? col_vals[j] : -col_vals[j];;
	}

	permute = atoi(argv[2]);
	switch( permute ){
	case 1:
	  printf("No permutation\n");
	  break;
	case 2:
	  printf("Columns permutation\n");
	  randPermColumns(col_ptrs, col_ids, col_vals, n);
	  break;
	case 3:
	  printf("Rows permutation\n");
	  randPermRows(col_ptrs, col_ids, col_vals, m);
	  break;
	case 4:
	  printf("Columns and Rows permutation\n");
	  randPermColumns(col_ptrs, col_ids, col_vals, n);
	  randPermRows(col_ptrs, col_ids, col_vals, m);
	  break;
	default:
	  printf("Potential options for third parameter are 1 (no permutation), 2 (columns permutation), 3 (rows permutation) and 4 (both columns and rows permutation)\n");
	  exit(12);
	}
	{
		int *			tmpspace = (int *) malloc(sizeof(int) * (m+n+1));
		t0 = u_wseconds();	       
		maxcmatching =  sprank(col_ptrs, col_ids,  n, m, tmpspace);
		t1 = u_wseconds();
		free(tmpspace);
		
	}
       /* report results code */
	minmn = m > n ? n : m;
mdeg = findMinDegs(col_ids, col_ptrs, n, m);
	printf("%s %d %d %.2f %d\n", argv[1], minmn, mdeg, t1-t0, maxcmatching);

	free(col_vals);
	free(col_ptrs);
	free(col_ids);

	return 0;
}
