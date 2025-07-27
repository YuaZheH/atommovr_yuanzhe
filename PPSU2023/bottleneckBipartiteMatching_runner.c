#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>

#include "portdefs.h"
#include "matrixUtils.h"

static double u_wseconds(void) 
{
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return (double) tp.tv_sec + (double) tp.tv_usec / 1000000.0;
}


#ifndef MAIN_C

#include "mex.h"
/*inputs*/
#define A_IN            ( prhs[0] )
/*outputs*/
#define CPRM_OUT       ( plhs[0] )
#define TIME_OUT       ( plhs[1] )

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int i, j, nz;

	mwSize _m, _n, _nz;
	mwIndex *irn_in;
	mwIndex *jcn_in;

	int m, n, maxcardmatching;
	int *match;
	int *row_match;
	int* col_ptrs, *row_ptrs;
	int* col_ids, *row_ids;
	int* fend_cols, *fend_rows;
	
	double *tmp_pr, *col_vals, *row_vals, *val_in;
	double etime, thrshld;

	if (!(nrhs == 1  && (nlhs == 1 || nlhs == 2))) {
		mexErrMsgTxt("not ok");
	}

	if (!mxIsSparse(A_IN)) {
		mexErrMsgTxt("bottleneckBipartiteMatching.c: First parameter must be a sparse matrix");
	}

	_n = mxGetN(A_IN);
	_m = mxGetM(A_IN);
	jcn_in = mxGetJc(A_IN);
	irn_in = mxGetIr(A_IN);
	val_in = mxGetPr(A_IN);

	_nz = jcn_in[_n];

	n = (int)_n; m = (int) _m; nz = (int) _nz;
	if(n != _n || m != _m || nz != _nz) {
		mexErrMsgTxt("bottleneckBipartiteMatching.c: problems with conversion...will abort.");
	}

	col_ptrs = (int *) mxCalloc(sizeof(int), (n+1));
	col_ids = (int *) mxCalloc(sizeof(int), nz);
	col_vals = (double*) mxCalloc(sizeof(double), nz);
	for(j = 0; j<=n; j++) {
		col_ptrs[j] = (int) jcn_in[j];
	}

	for(j = 0; j<nz; j++) {
		col_ids[j] = (int) irn_in[j];
		col_vals[j] = val_in[j] >=0.0 ? val_in[j] : -val_in[j];;
	}

	match = (int*) mxCalloc(sizeof(int), n);
	row_match = (int*) mxCalloc(sizeof(int), m);
	row_ptrs = (int *)malloc(sizeof(int) * (m+1));
	row_ids = (int *) malloc(sizeof(int) * nz );
	row_vals = (double *) malloc(sizeof(double) * nz);
	fend_cols = (int *) malloc(sizeof(int) * n);
	fend_rows = (int *) malloc(sizeof(int) * m);

	{
		int *tmpspace = (int *) malloc(sizeof(int) * (m+n+1));
		maxcardmatching = sprank(col_ptrs, col_ids,  n,  m, tmpspace);
		free(tmpspace);
	}
/*do computation*/
	double t0 = u_wseconds();


 		bttlThreshold(col_ptrs, col_ids, col_vals, n, m, match, row_match, 
		row_ptrs, 
		row_ids,
		row_vals,
		fend_cols,fend_rows,
		    1, &thrshld, maxcardmatching);

double t1 = u_wseconds();
 	printf("total time %.2f, thrshld %.4f\n", t1-t0, thrshld);


/*done computing*/
	CPRM_OUT = mxCreateDoubleMatrix(m, 1, mxREAL);
	tmp_pr = mxGetPr(CPRM_OUT);
	for (i = 0; i < m; i++) 
	{	
		if(row_match[i] != -1) 
			tmp_pr[i] = row_match[i] + 1;	
		else
			tmp_pr[i] =-1;
	}
	if(nlhs == 2)
	{
		TIME_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
		tmp_pr = mxGetPr(TIME_OUT);
		*tmp_pr = t1-t0;
	}

	free(fend_rows);
	free(fend_cols);
	free(row_vals);
	free(row_ids);
	free(row_ptrs);

	mxFree(row_match);
	mxFree(match);
	mxFree(col_vals);
	mxFree(col_ptrs);
	mxFree(col_ids);
	return;
}
#else
#include <string.h>


int main(int argc, char *argv[])
{
	int j, nz;

	int m, n;
	int *match;
	int *row_match;
	double *col_vals, *row_vals, thrshld;
	double t0, t1;

	int* col_ptrs, *row_ptrs;
	int* col_ids, *row_ids;
	int* fend_cols, *fend_rows;
	
	int optScaling, permute, maxcardmatching ;

if (argc != 4)
	{
		printf("need three args\n mtx_filename \nscale_opt(noscale 0, scale-pat 2, scale-A 3)\n");
		printf("\tpermutation (1,2,3, 4): no, col, row, row and col\n");
		exit(12);
	}

	t0 = u_wseconds();
	ReadMatrixMarket(argv[1], &col_vals, &col_ids, &col_ptrs, &m,	&n, &nz,  atoi(argv[2]) != 2);
	t1 = u_wseconds();
	printf("read-in time %.2f\n", t1-t0);

	for(j = 0; j<nz; j++) 
	{
		col_vals[j] = col_vals[j] >=0 ? col_vals[j] : -col_vals[j];;
	}

	optScaling = atoi(argv[2]);

	if(optScaling == 2)
	{
		printf("\tpattern scaling\n");
		ScaleOnePattern(	 col_ptrs,  col_ids, col_vals,  n, m, 20);
	}
	else if(optScaling == 3)
	{
		printf("\tval scaling\n");
		ScaleOne(col_ptrs,  col_ids, col_vals,  n,  m,20);
	}
	else
	{
		printf("\tNo scaling\n");
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

	match = (int*) calloc(sizeof(int), n);
	row_match = (int*) calloc(sizeof(int), m);
	row_ptrs = (int *)malloc(sizeof(int) * (m+1));
	row_ids = (int *) malloc(sizeof(int) * nz );
	row_vals = (double *) malloc(sizeof(double) * nz);
	fend_cols = (int *) malloc(sizeof(int) * n);
	fend_rows = (int *) malloc(sizeof(int) * m);
	{
		int *tmpspace = (int *) malloc(sizeof(int) * (m+n+1));
		maxcardmatching = sprank(col_ptrs, col_ids,  n,  m, tmpspace);
		free(tmpspace);
	}

/*do computation*/
	t0 = u_wseconds();

	bttlThreshold(col_ptrs, col_ids, col_vals, n, m, match, row_match, 
		row_ptrs, 
		row_ids,
		row_vals,
		fend_cols,fend_rows,
		    1, &thrshld, maxcardmatching);

	t1 = u_wseconds();
	printf("total thrshld time %.2f, thrshld %.4f\n", t1-t0, thrshld);

/*done computing*/
	free(fend_rows);
	free(fend_cols);
	free(row_vals);
	free(row_ids);
	free(row_ptrs);

	free(row_match);
	free(match);
	free(col_vals);
	free(col_ptrs);
	free(col_ids);

	return 0;
}

#endif
