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
#define ITERBTTL_OUT   ( plhs[2] )
#define ITERBSCT_OUT   ( plhs[3] )
#define NUMAUG_OUT     ( plhs[4] )

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int i, j, nz;

	mwSize _m, _n, _nz;
	mwIndex *irn_in;
	mwIndex *jcn_in;

	int m, n;
	int *match;
	int *row_match;
	int* col_ptrs, *row_ptrs;
	int* col_ids, *row_ids;
	int* fend_cols, *fend_rows;
	
	double *tmp_pr, *col_vals, *row_vals, *val_in;
	double etime;
	double thrshld, thrshld_puremc64j3, t0, t1;
int iterBttl, iterBisect, numAugs, maxcardmatching;
	if (!(nrhs == 1  && (nlhs == 1 || nlhs == 2 || nlhs == 3 || nlhs == 4 || nlhs == 5))) {
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
	


	myPrintf("total thrshld time %.2f, thrshld %.4f, org %.4f\n", t1-t0, thrshld_puremc64j3, thrshld);
	iterBttl  = bttlThreshold(col_ptrs, col_ids, col_vals, n, m, match, row_match, 
		row_ptrs, 
		row_ids,
		row_vals,
		fend_cols,fend_rows,
		    1, &thrshld, maxcardmatching);
	
			t0 = u_wseconds();
	iterBisect  = bisectionBasedOnMC64J3(col_ptrs, col_ids, col_vals, n,  m, match, row_match, 
	row_ptrs, 
	row_ids,
	row_vals,
	fend_cols, fend_rows, &thrshld_puremc64j3, maxcardmatching);

	t1 = u_wseconds();

	if(thrshld != thrshld_puremc64j3 )
	{
		myPrintf("the values are not equal %.4f %.4f diff %.6f\n", thrshld, thrshld_puremc64j3, thrshld - thrshld_puremc64j3);
		if (iterBttl>0)
			iterBttl = -iterBttl;
		iterBisect = -iterBisect;
	}
/*
	t0 = u_wseconds();
	numAugs = pureSAP(col_ptrs, col_ids, col_vals, n,  m, match, row_match, 
		row_ptrs, 
		row_ids,
		row_vals,
		fend_cols, fend_rows, &thrshld, maxcardmatching);
	t1 = u_wseconds();
	printf("pureSAP total time %.2f iters %d\n", t1-t0, numAugs);
	if(thrshld != thrshld_puremc64j3)
	{
		myPrintf("the values are not equal for pureSAP %.4f %.4f diff %.6f\n", thrshld, thrshld_puremc64j3, thrshld - thrshld_puremc64j3);
		if (iterBttl>0)
			iterBttl = -iterBttl;
	iterBisect = -iterBisect;

	}
*/

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
	if(nlhs >= 2)
	{
		TIME_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
		tmp_pr = mxGetPr(TIME_OUT);
		*tmp_pr = t1-t0;
	}
	if(nlhs >= 3)
	{
		ITERBTTL_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
		tmp_pr = mxGetPr(ITERBTTL_OUT);
		*tmp_pr = iterBttl;
	}
	if(nlhs >= 4)
	{
		ITERBSCT_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
		tmp_pr = mxGetPr(ITERBSCT_OUT);
		*tmp_pr = iterBisect;
	}
	if(nlhs >= 5)
	{
		NUMAUG_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
		tmp_pr = mxGetPr(NUMAUG_OUT);
		*tmp_pr = numAugs;
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
	int *row_match, maxcardmatching;
	double *col_vals, *row_vals, *org_col_vals, thrshld, thrshld_puremc64j3;
	double t0, t1;

	int* col_ptrs, *row_ptrs, *org_col_ids;
	int* col_ids, *row_ids;
	int* fend_cols, *fend_rows;
	
	int iterBttl, iterBisect;
	int optScaling, permute ;

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

	//randPermColumns(col_ptrs, col_ids, col_vals,  n);

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
	org_col_ids = (int *) malloc(sizeof(int) * nz);
	org_col_vals = (double *) malloc(sizeof(double) * nz);
	memcpy(org_col_ids, col_ids, sizeof(int) * nz);
	memcpy(org_col_vals, col_vals, sizeof(double) * nz);

	t0 = u_wseconds();
	iterBttl = bttlThreshold(col_ptrs, col_ids, col_vals, n, m, match, row_match, 
		row_ptrs, 
		row_ids,
		row_vals,
		fend_cols,fend_rows,
		    1, &thrshld, maxcardmatching);
	t1 = u_wseconds();
	printf("bttlThreshold total time %.2f iters %d vals %.6f\n", t1-t0, iterBttl, thrshld);

	memcpy(col_ids, org_col_ids, sizeof(int) * nz);
	memcpy(col_vals, org_col_vals, sizeof(double) * nz);

	t0 = u_wseconds();
	iterBisect = bisectionBasedOnMC64J3(col_ptrs, col_ids, col_vals, n,  m, match, row_match, 
		row_ptrs, 
		row_ids,
		row_vals,
		fend_cols, fend_rows, &thrshld_puremc64j3, maxcardmatching);
	t1 = u_wseconds();

	printf("bisection total time %.2f iters %d\n", t1-t0, iterBisect);
	if(thrshld != thrshld_puremc64j3)
	{
		printf("***********************************************************************\n");
		printf("*\t\tthe values are not equal %.4f %.4f %.10f\t\t*\n", thrshld, thrshld_puremc64j3, thrshld - thrshld_puremc64j3);
		printf("***********************************************************************\n");
	}
/*
*/
/*done computing*/
	free(org_col_vals);
	free(org_col_ids);
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
