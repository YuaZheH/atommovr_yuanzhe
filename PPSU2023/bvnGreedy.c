#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <float.h>

#include "portdefs.h"
#include "matrixUtils.h"


//#define SANITY_DEBUG 1


static double u_wseconds(void) 
{
	struct timeval tp;

	gettimeofday(&tp, NULL);

	return (double) tp.tv_sec + (double) tp.tv_usec / 1000000.0;
}

void transpose(int* sptrs, int* sids, double *svals, int* tptrs, int* tids, double *tvals, int n, int m) ;
void sortValsInsSortQS(int* ptrs, int* ids, double *vals, int n, int m);

void checkDoublyStochastic(int *col_ptrs, int *col_ids, double *col_vals, int n, double *my_zero, double *my_one)
{
	int j, i;

	double my_error = -1.0, lmyone = -1.0;
	for(j = 0; j<n; j++) 
	{
		
		double mysum = 0.0;
		for (i = col_ptrs[j]; i < col_ptrs[j+1]; i++)
		{
			col_vals[i] = col_vals[i] >=0 ? col_vals[i] : -col_vals[i];
			mysum += col_vals[i];
		}	
		if (mysum > lmyone)
			lmyone = mysum;
	}	
	for(j = 0; j<n; j++) 
	{
		
		double mysum = 0.0;
		for (i = col_ptrs[j]; i < col_ptrs[j+1]; i++)
		{
			col_vals[i] = col_vals[i] >=0 ? col_vals[i] : -col_vals[i];
			mysum += col_vals[i];
		}		



		if(fabs(mysum - lmyone) > my_error)
			my_error = fabs(mysum - lmyone);
	}
	*my_one = lmyone;	
	*my_zero = my_error * 1.1;
}

void inssortdecreades(int *ids, double *vals, int sz)
{
/*we sort until we see -1, or sz*/
/*vals already sorted in descending order; there are potentially -1.0's at the end*/
	if( sz == 1 ) return;

	int k, idsave = ids[0];
	double v = vals[0];

	for (k = 1; k < sz && vals[k] > v ; k++)
	{
		ids[k-1] = ids[k];
		vals[k-1] = vals[k];
	}
	ids[k-1] = idsave;
	vals[k-1] = v; 
}

void updateValues(int *ptrs, int *ids, double *vals, int *fends, int n, int *match, double *threshold, double my_zero)
{
	int i, j, other;
	double thnew = 1.5 * *threshold;
#ifdef SANITY_DEBUG
    int matchok;
#endif

	for(i = 0; i < n; i++)
	{
#ifdef SANITY_DEBUG
		matchok = 0;
#endif
		for( j = ptrs[i]; j <= fends[i]; j++)/*find the mate of i.*/
		{
			other = ids[j];
			
			if(match[i] == other)
			{
				vals[j] = (vals[j] - *threshold <= my_zero+1.0e-16) ?  -1.0 :  vals[j] - *threshold;
	//			vals[j] =  vals[j] - *threshold;
				/*sort from j to ptrs ptrs[i+1]*/
				inssortdecreades(ids+j, vals+j, ptrs[i+1]-j);

				/*fend can only go left by one* */
				if(vals[fends[i]] < *threshold)
					fends[i] --;
				if(fends[i]  < ptrs[i])
				{
					if (thnew > vals[ptrs[i]] && vals[ptrs[i]] >my_zero+1.0e-16)/*the largest element in each row/col is a potential threshold*/
					{
						thnew = vals[ptrs[i]];
						fends[i] = ptrs[i]	;
					}
					else if (vals[ptrs[i]] <= my_zero+1.0e-16) thnew = -1.0;
				}

#ifdef SANITY_DEBUG
				matchok = 1;
#endif								

				break;
			}
		}
#ifdef SANITY_DEBUG
		if (matchok ==0 )
		{
			myExit("bvnGreedy: match is not ok\n");
		}
#endif		
	}
	if(thnew < *threshold)
		*threshold = thnew;

}

void printMatching(int *mtch, int sz)
{
	int i;
	if(sz > 10)
		return;
	for (i = 0; i < sz; i++)
	{
		myPrintf("%d\n", mtch[i]);
	}
}

//#define _MATRIX_IS_EMPTY_ ( (currentThreshold  <= my_zero ) || fabs(totalAlpha - my_one) <= my_zero ) 
#define _MATRIX_IS_EMPTY_ ( (currentThreshold <= my_zero+1.0e-16 ) || fabs(totalAlpha - my_one) <= my_zero + 1.0e-16 ) 

void bvnGreedyInterface(int *col_ptrs, int *col_ids, double *col_vals, int n, int m, double my_zero, double my_one, 
	int **row_matches_sol, double **alphas_sol, int *num_factors)
{

	int *row_matches ;/*this is an alias to *row_matches_sol to be able to allocate space*/

	int *row_ptrs, *row_ids;
	double *row_vals;
	double currentThreshold, *alphas, totalAlpha;
	int nz = col_ptrs[n];
	int *fend_cols, *fend_rows;
	int *match, *row_match;
	int num_factors_in = *num_factors;
	int numComputed, rmsz ;

	if (num_factors_in < 1)
	{
		row_matches = (int *) malloc(sizeof(int) * n ); /*we start with 1; we will reallocate later*/
		rmsz = 1;
	}
	else
	{
		row_matches = (int *) malloc(sizeof(int) * n * num_factors_in);
		rmsz = num_factors_in;
	}
	alphas = (double *) malloc(sizeof(double) * rmsz);

	row_ptrs = (int *)malloc(sizeof(int) * (m+1));
	row_ids = (int *) malloc(sizeof(int) * nz );
	row_vals = (double *) malloc(sizeof(double) * nz);
	fend_cols = (int *) malloc(sizeof(int) * n);
	fend_rows = (int *) malloc(sizeof(int) * m);

	bttlThresholdInitializer(col_ptrs,  col_ids, col_vals,  n,  m,  nz, 
		row_ptrs, row_ids, row_vals,
		fend_cols, fend_rows,
		&currentThreshold, n);


	/*check input*/

	match = (int*) malloc(sizeof(int) * n);
	row_match = (int*) malloc(sizeof(int) * m);

	numComputed = -1;
	totalAlpha = 0.0;
	myPrintf("this is the current threshold %.5f (%.5f %.5f)\n", currentThreshold, my_zero, my_one);
	while ( (num_factors_in == 0 ||(num_factors_in  > 0 && numComputed < num_factors_in-1)) &&  ! (_MATRIX_IS_EMPTY_))
	{
		double thcols, throws;
		numComputed ++;
//		myPrintf("inside: this is the current threshold %.16f (%.16f) %d\n", currentThreshold, my_zero,currentThreshold <= my_zero  );

		if (numComputed >= rmsz)
		{
			row_matches = realloc(row_matches, n * 2 * rmsz  *sizeof(int));
			if(row_matches == NULL)
			{
				myPrintf("requested (potentially) %d factors, cannot allocate\n", rmsz * 2);
				myPrintf("suggested num_factors as an input argument is %d\n", rmsz );
				myExit("see above");
			}

			alphas = realloc(alphas, 2*rmsz*sizeof(double));
			if (alphas == NULL)
			{
				myPrintf("requested (potentially) %d coeffiecitns, cannot allocate\n", rmsz * 2);
				myPrintf("suggested num_factors as an input argument is %d\n", rmsz );
				myExit("see above");
			}
			rmsz = 2 * rmsz;
		}
		/*we can compute another perfect factor now into row_matches[numComputer * m] */
		int numIters = bttlThreshold(col_ptrs, col_ids, col_vals, n, m, match, row_matches + numComputed * m, 
			row_ptrs, row_ids, row_vals,fend_cols, fend_rows, 0, &currentThreshold, n);
 	/*	
 		printMatching(row_matches + numComputed * m, m);
 	*/
		if(currentThreshold == DBL_MAX ||  numIters < 0)
		{
			numComputed -- ;
			break;
		}			
	//	myPrintf("\tfound one factor with %.4f\n", currentThreshold);

		/*updates*/
		/* 1 : total alpha*/
		alphas[numComputed] = currentThreshold;
		totalAlpha += currentThreshold;

		/* 2 : subtract from cols, subtract from rows*/
		throws = thcols = currentThreshold;
		updateValues(col_ptrs, col_ids, col_vals, fend_cols, n, match, &thcols,  my_zero);
		updateValues(row_ptrs, row_ids, row_vals, fend_rows, m, row_matches + numComputed * m, &throws,  my_zero);

	//	myPrintf("\throws %.10f, cols %.10f\n", thcols, throws);
		if(throws >= my_zero && thcols >= my_zero)
		{
		if (throws < currentThreshold)
			currentThreshold = throws;
		if (thcols < currentThreshold)
			currentThreshold = thcols;
	//	myPrintf("\tvalue of the next thrshold %.10f, total alpha %.4f at %d (zero %.10f)\n", currentThreshold, totalAlpha, numComputed+1, my_zero);
	}
	else
		currentThreshold = -1.0;
	//	getc(stdin);
	}
	*num_factors = numComputed+1;
	*row_matches_sol = row_matches;
	*alphas_sol = alphas;
	free(row_match);
	free(match);
	free(fend_rows);
	free(fend_cols);	
	free(row_vals);
	free(row_ids);
	free(row_ptrs);
}

#ifndef MAIN_C
#include "mex.h"
/*inputs*/
#define A_IN               ( prhs[0] )
#define MY_ZERO            ( prhs[1] )
#define MY_ONE             ( prhs[2] )
#define NUM_FACTORS        ( prhs[3] )

/*outputs*/
#define CPRM_OUT         ( plhs[0] )
#define ALPHAS_OUT       ( plhs[1] )
#define TIME_OUT         ( plhs[2] )

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	int i, j, nz;

	mwSize _m, _n, _nz;
	mwIndex *irn_in;
	mwIndex *jcn_in;
	double my_zero = 1.0e-6, my_one = 1.0;
	int num_factors = 16;

	int m, n;
	int* col_ptrs;
	int* col_ids;
	double *tmp_pr, *col_vals, *val_in;
	double etime;
	double *alphas; 
	int *row_matches;
	double t0, t1;

	if (!( nrhs >= 2 &&  nrhs <=4  && (nlhs == 2 || nlhs == 3))) {
		mexErrMsgTxt("[cprms, alphas<, TIME(.s)> ]= bvnGreedy(A, myzero, myone, numfactors)");
	}

	if (!mxIsSparse(A_IN)) {
		mexErrMsgTxt("bvnGreedy: First parameter must be a sparse matrix");
	}

	_n = mxGetN(A_IN);
	_m = mxGetM(A_IN);
	jcn_in = mxGetJc(A_IN);
	irn_in = mxGetIr(A_IN);
	val_in = mxGetPr(A_IN);

	_nz = jcn_in[_n];

	n = (int)_n; m = (int) _m; nz = (int) _nz;
	if(n != _n || m != _m || nz != _nz) {
		mexErrMsgTxt("bvnGreedy: problems with conversion...will abort.");
	}
	if(n != m) {
		mexErrMsgTxt("bvnGreedy: we expect a square sparse matrix\n");
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
    

	if(nrhs == 2)
	{
		tmp_pr = mxGetPr(MY_ZERO);
		my_zero = *tmp_pr;	
	}
	else if (nrhs == 3)
	{
		tmp_pr = mxGetPr(MY_ZERO);
		my_zero = *tmp_pr;

		tmp_pr = mxGetPr(MY_ONE);
		my_one = *tmp_pr;
	}
	else if (nrhs == 4)
	{
	tmp_pr = mxGetPr(MY_ZERO);
		my_zero = *tmp_pr;

		tmp_pr = mxGetPr(MY_ONE);
		my_one = *tmp_pr;

		tmp_pr = mxGetPr(NUM_FACTORS);
		num_factors = (int) ceil(*tmp_pr);
	}



	checkDoublyStochastic(col_ptrs, col_ids, col_vals, n, my_zero, my_one);

/*do computation*/
t0 = u_wseconds();
 bvnGreedyInterface(col_ptrs, col_ids, col_vals,  n,  m, my_zero,  my_one,  	
	&row_matches, &alphas, &num_factors);
 t1 = u_wseconds();
myPrintf("num matches %d in %.2f seconds\n", num_factors, t1-t0);

/*done computing*/
	CPRM_OUT = mxCreateDoubleMatrix(m*num_factors, 1, mxREAL);
	tmp_pr = mxGetPr(CPRM_OUT);
	for (j = 0; j < num_factors; j++)
	{
		for (i = 0; i < m; i++) 
		{	
			if(row_matches[n*j + i ] >= 0 && row_matches[n*j + i ] < n) 
				tmp_pr[ n*j + i ] = row_matches[n*j + i ] + 1;	
			else
				mexErrMsgTxt("bvngreedy: not perfect matching\n");
		}
	}
	ALPHAS_OUT = mxCreateDoubleMatrix(num_factors, 1, mxREAL);
	tmp_pr = mxGetPr(ALPHAS_OUT);
	for (j = 0; j < num_factors; j++)
		tmp_pr[j] = alphas[j];

	if(nlhs == 3)
	{
		TIME_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
		tmp_pr = mxGetPr(TIME_OUT);
		*tmp_pr = t1-t0;

	}
	free(alphas);
	free(row_matches);
	mxFree(col_vals);
	mxFree(col_ptrs);
	mxFree(col_ids);
	return;
}
#else

#include <string.h>
//#include "mmio.h"


double sumalphas(double *alphas, int sz)
{
	double ans = 0.0;
	int i;
	for (i = 0; i < sz; i ++)
		ans += alphas[i];
	return ans;
}
void printlphas(double *alphas, int sz)
{
	int i;
	for (i = 0; i < sz; i ++)
		myPrintf("%d %.16f\n", i+1, alphas[i]);

}
int main(int argc, char *argv[])
{

	int j, nz;

	int m, n;
	int* col_ptrs;
	int* col_ids;
	int num_factors = 0;
	double *col_vals;
	double my_zero = 1.0e-6, my_one = 1.0;
	double t0, t1;
	double *alphas; int *row_matches;
	if (argc < 4 || argc > 7)
	{
		printf("need at least 3 additional arg (mtx filename) and at most 6\n");
		printf(" scaling (1=No, 2=Pattern, 3=val), permutation (1=No, 2=cols, 3=rows, 4=rowscols) \n");
		printf(" myZero, myOne --- these are ignored for now\n");
		printf("the number of requested permutation matrices; of not bvnGreedy finds\n");

		exit(12);
	}

	if(argc == 5)
	{
		my_zero = (double) atof(argv[4]);	
	}
	else if (argc == 6)
	{
		my_zero = (double) atof(argv[4]);	
		my_one = (double) atof(argv[5]);	
	}
	else if (argc == 7)
	{
		my_zero = (double) atof(argv[4]);	
		my_one = (double) atof(argv[5]);	
		num_factors = atoi(argv[6]);
	}


	/* Reading the matrix */
	int permute = 0;
	int scaling = 0;
	scaling = atoi(argv[2]);
	t0 = u_wseconds();
	ReadMatrixMarket(argv[1], &col_vals, &col_ids, &col_ptrs, &m,	&n, &nz, !(scaling==2));
	t1 = u_wseconds();
	printf("read-in time %.2f\n", t1-t0);
	if(m != n)
	{
		printf("we expect a square matrix\n");
		exit(12);
	}

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
	  randPermRows(col_ptrs, col_ids, col_vals, n);
	  break;
	case 4:
	  printf("Columns and Rows permutation\n");
	  randPermColumns(col_ptrs, col_ids, col_vals, n);
	  randPermRows(col_ptrs, col_ids, col_vals, n);
	  break;
	default:
	  printf("Potential options for third parameter are 1 (no permutation), 2 (columns permutation), 3 (rows permutation) and 4 (both columns and rows permutation)\n");
	  exit(12);
	}


	/*check input*/
	checkDoublyStochastic(col_ptrs, col_ids, col_vals, n, &my_zero, &my_one);

	myPrintf("the params zero, one, numf %4f %.2f %d\n", my_zero, my_one, num_factors);

/*do computation*/
	t0 = u_wseconds();
	bvnGreedyInterface(col_ptrs, col_ids, col_vals,  n,  m, my_zero,  my_one, 
		&row_matches, &alphas, &num_factors);
	t1 = u_wseconds();

	printf("total bvnGreedy time %.2f, numfactors computed %d, sum alphas %.4f\n", t1-t0, num_factors, sumalphas(alphas, num_factors));
	//printlphas(alphas, num_factors);

/*done computing*/

	free(alphas);
	free(row_matches);
	free(col_vals);
	free(col_ptrs);
	free(col_ids);

	return 0;
}

#endif
