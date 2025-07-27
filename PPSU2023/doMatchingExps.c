#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include <float.h>

#include "mmio.h"
#include "portdefs.h"
#include "matrixUtils.h"


void  mc64id_(int *, double *);

void  mc64ad_(int *, int *, int *, int *, int *, int *,double *,
	      int *, int *, int *,int *, int *, double *,
	      int *, double *, int *);




static double u_wseconds(void) 
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return (double) tp.tv_sec + (double) tp.tv_usec / 1000000.0;
}



void mc64_wrapper(int n, int m, int nz, int* col_ptrs, int* col_ids, double* col_vals, int job) 
{

	/* Run MC64 */

	  double t0, t1;

	  int i;
	  int nicntl, ninfo;
	  int /*job,*/ num, liw, ldw;
	  int *perm, *iw;
	  int *icntl, *info;
	  double *dw = NULL, *cntl;

	  /* Fortran indices */
	  for(i = 0; i < n+1; i++){
	    col_ptrs[i] = col_ptrs[i]+1;
	  }
	  for(i = 0; i < nz; i++){
	    col_ids[i] = col_ids[i]+1;
	  }

	  /* fix mc64 control parameters */
	  nicntl = 10;
	  ninfo  = 10;
	  icntl = calloc(nicntl,sizeof(int));
	  info  = calloc(ninfo,sizeof(int));
	  cntl  = calloc(ninfo,sizeof(double));

	  mc64id_(icntl, cntl);

	  icntl[0] = -1;
	  icntl[1] = -1;
	  /* allocation of mc64 workspace */
	  if ( job == 2 ){
	    liw = 2*n+2*m;
	    ldw = m;
	  } else if(job == 3 ){
	    liw = 8*n+2*m+nz;
	    ldw = nz;
	  } else{
	    free(icntl);
	    free(info); 
	    free(cntl); free(col_vals); free(col_ids); free(col_ptrs);
	    printf ("Wrong job value\n");
	    exit(12);
	  } 
	  iw = calloc(liw, sizeof(int));
	  if (ldw > 0){
	    dw = calloc(ldw, sizeof(double));
	  }

	  perm = calloc(m, sizeof(int));

	  t0 = u_wseconds();
	  mc64ad_(&job, &m, &n, &nz,
		  col_ptrs, col_ids, col_vals,
		  &num, perm,
		  &liw, iw, &ldw, dw,
		  icntl, cntl, info);
	  t1 = u_wseconds();
	  printf("total mc64 time %.2f\n", t1-t0);

	  free(icntl);
	  free(info); 
	  free(cntl);
	  free(iw);
	  if (ldw > 0){
	    free(dw);
	  }
	  free(perm);
	  
	  /* Restore C indices */
	  for(i = 0; i < n+1; i++){
	    col_ptrs[i] = col_ptrs[i]-1;
	  }
	  for(i = 0; i < nz; i++){
	    col_ids[i] = col_ids[i]-1;
	  }


}

int main(int argc, char *argv[])
{

	int j, nz;

	int m, n;
	int *match;
	int *row_match;
	double *col_vals, *row_vals, thrshld, *org_col_vals, thrshld_puremc64j3;
	double t0, t1;

	int* col_ptrs, *row_ptrs;
	int* col_ids, *row_ids, *org_col_ids;
	int* fend_cols, *fend_rows;
	int iterBttl, iterBisect;
	int permute = 0;
	int scaling = 0;
	int maxcmatching = 0;
	if ( argc != 5 ){
		printf("Usage ./doMatchingExps matrix_name scaling permute_columns mc64-job_id\n");
		printf("\tscaling (1, 2, 3): no, pattern, value\n");
		printf("\tpermutation (1,2,3, 4): no, col, row, row and col\n");
		printf("\tmc64-job_id (2 or 3): sap, threshold-based\n");
		exit(12);
	}
	printf("Performing experiments on matrix %s with scaling %s and permutation %s\n",
	       argv[1], argv[2], argv[3]);

	/* Reading the matrix */
	t0 = u_wseconds();

	scaling = atoi(argv[2]);

	ReadMatrixMarket(argv[1], &col_vals, &col_ids, &col_ptrs, &m,	&n, &nz, !(scaling==2));
	t1 = u_wseconds();
	printf("read-in time %.2f\n", t1-t0);

	for(j = 0; j<nz; j++) 
	{
	  col_vals[j] = col_vals[j] >=0 ? col_vals[j] : -col_vals[j];;
	}

	switch( scaling ){
	case 1:
	  printf("No Scaling performed\n");
	  break;
	case 2:
	  printf("Should scale the matrix pattern\n");
	  ScaleOnePattern( col_ptrs, col_ids, col_vals, n, m, 20);
	  break;
	case 3:
	  printf("Should scale the values\n");
	  ScaleOne( col_ptrs, col_ids, col_vals, n, m, 20);
	  break;
	default:
	  printf("Potential options for second parameter are 1 (no scaling), 2 (scale the pattern) or 3 (scale the values)\n");
	  exit(12);
	}

	permute = atoi(argv[3]);
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
		maxcmatching =  sprank(col_ptrs, col_ids,  n, m, tmpspace);
		free(tmpspace);
	}
       /* run mc64 code */
	if(m < n)
	{
		printf("MC64 needs m (%d) >=n (%d)\n", m, n);
		exit(12);
	}
	int job_id = atoi(argv[4]); 
	
	org_col_ids = (int *) malloc(sizeof(int) * nz);
	org_col_vals = (double *) malloc(sizeof(double) * nz);
	memcpy(org_col_ids, col_ids, sizeof(int) * nz);
	memcpy(org_col_vals, col_vals, sizeof(double) * nz);
	if(maxcmatching == n && job_id == 2)
		mc64_wrapper(n,m,nz,col_ptrs,col_ids,col_vals,job_id); 
	else
		printf("\tmc64job-2 is not run, because no full sprank\n");

	{
	  match = (int*) calloc(sizeof(int), n);
	  row_match = (int*) calloc(sizeof(int), m);

	  row_ptrs = (int *)malloc(sizeof(int) * (m+1));
	  row_ids = (int *) malloc(sizeof(int) * nz );
	  row_vals = (double *) malloc(sizeof(double) * nz);
	  fend_cols = (int *) malloc(sizeof(int) * n);
	  fend_rows = (int *) malloc(sizeof(int) * m);

	  memcpy(col_ids, org_col_ids, sizeof(int) * nz);
	  memcpy(col_vals, org_col_vals, sizeof(double) * nz);
	  thrshld = thrshld_puremc64j3 = DBL_MAX;
	  t0 = u_wseconds();
	  iterBttl = bttlThreshold(col_ptrs, col_ids, col_vals, n, m, match, row_match, 
		      row_ptrs, 
		      row_ids,
		      row_vals,
		      fend_cols,fend_rows,
		      1, &thrshld, maxcmatching);
	  t1 = u_wseconds();
	  printf("total bottled time %.2f, threshold %.4f iters %d\n", t1-t0, thrshld, iterBttl);


	memcpy(col_ids, org_col_ids, sizeof(int) * nz);
	memcpy(col_vals, org_col_vals, sizeof(double) * nz);

	t0 = u_wseconds();
	iterBisect = bisectionBasedOnMC64J3(col_ptrs, col_ids, col_vals, n,  m, match, row_match, 
		row_ptrs, 
		row_ids,
		row_vals,
		fend_cols, fend_rows, &thrshld_puremc64j3, maxcmatching);
	t1 = u_wseconds();

	printf("total thresh time %.2f, threshold %.4f iters %d\n", t1-t0, thrshld_puremc64j3, iterBisect);

	  free(fend_rows);
	  free(fend_cols);
	  free(row_vals);
	  free(row_ids);
	  free(row_ptrs);
	  free(row_match);
	  free(match);
	}
	free(org_col_vals);
	free(org_col_ids);
	free(col_vals);
	free(col_ptrs);
	free(col_ids);

	return 0;
}
