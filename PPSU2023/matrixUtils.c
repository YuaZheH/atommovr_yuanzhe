#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "portable_time.h"
#include "mmio.h"

/*we do a lousy random permutation*/
static inline void shuffle(int *a, int n)
{
  int i;
  for (i = n-1; i >= 0; --i)
  {
    //generate a random number [0, n-1]
    int j = rand() % (i+1);

    //swap the last element with element at random index
    int temp = a[i];
    a[i] = a[j];
    a[j] = temp;
  }
}

void randPermRows(int *col_ptrs, int * col_ids, double *col_vals, int m)
{
	int i, j;
	int *prm, *invprm;
	
	prm = (int *) malloc(sizeof(int) * m);
	invprm = (int *) malloc(sizeof(int) * m);
	for (i = 0; i < m; i++)
		prm[i] = i;

	shuffle(prm, m);

	for (i = 0; i < m; i++)
		invprm[prm[i]] = i; 

	for (i = 0; i < m; i++)
	{		
		for( j = col_ptrs[i]; j < col_ptrs[i+1]; j++)
			col_ids[j] = invprm[col_ids[j]];		
	}

	free(invprm);
	free(prm);
}
void randPermColumns(int *col_ptrs, int * col_ids, double *col_vals, int n)
{
	int i, nz = col_ptrs[n];
	int *prm;
	
	int *ptrs_new, *ids_new;
	double *vals_new;

	srand((unsigned)time(NULL));
	prm = (int *) malloc(sizeof(int) * n);
	for (i = 0; i < n; i++)
		prm[i] = i;
	shuffle(prm, n);

	ptrs_new = (int *) malloc(sizeof(int) * (n+1));
	ids_new = (int *) malloc(sizeof(int) * nz);
	vals_new = (double *) malloc(sizeof(double) * nz);
	ptrs_new[0] = 0;
	for (i = 0; i < n; i++)
	{
		int c = prm[i];
		int j, pos = ptrs_new[i];
		for( j = col_ptrs[c]; j < col_ptrs[c+1]; j++, pos++)
		{
			ids_new[pos] = col_ids[j];
			vals_new[pos] = col_vals[j];
		}
		ptrs_new[i+1] = pos;
	}

	if( col_ptrs[n] != ptrs_new[n])
	{
		printf("array was not properly copied\n");
		exit(21);
	}
	/*copy to the original arrarys*/
	memcpy(col_ptrs, ptrs_new, sizeof(int) * (n+1));
	memcpy(col_ids, ids_new, sizeof(int) * nz);
	memcpy(col_vals, vals_new, sizeof(double) * nz);

	free(ptrs_new);
	free(ids_new);
	free(vals_new);
	free(prm);
}

/*
static double u_wseconds(void) 
{
	struct timeval tp;

	gettimeofday(&tp, NULL);

	return (double) tp.tv_sec + (double) tp.tv_usec / 1000000.0;
}
*/
void symbtranspose(int* sptrs, int* sids, int* tptrs, int* tids, int n, int m) 
{
	int i;
	int ptr, eptr;

	memset(tptrs, 0, sizeof(int)  * (m+1));

	for(i = 0; i < n; i++) 
	{
		eptr = sptrs[i+1];
		for(ptr = sptrs[i]; ptr < eptr; ptr++) 
			tptrs[sids[ptr] + 1]++;
	}
	for(i = 1; i < m; i++) 
		tptrs[i] += tptrs[i-1];

	for(i = 0; i < n; i++) 
	{
		eptr = sptrs[i+1];
		for(ptr = sptrs[i]; ptr < eptr; ptr++) 
			tids[tptrs[sids[ptr]]++] = i;            

	}
	for(i = m; i > 0; i--) 
		tptrs[i] = tptrs[i-1];
	tptrs[0] = 0;

}

 void transpose(int* sptrs, int* sids, double *svals, int* tptrs, int* tids, double *tvals, int n, int m) 
{
	int i;
	int ptr, eptr;
	memset(tptrs, 0, sizeof(int)  * (m+1));

	for(i = 0; i < n; i++) 
	{
		eptr = sptrs[i+1];
		for(ptr = sptrs[i]; ptr < eptr; ptr++) 
			tptrs[sids[ptr] + 1]++;
	}

	for(i = 1; i < m; i++) 
		tptrs[i] += tptrs[i-1];

	for(i = 0; i < n; i++) 
	{
		eptr = sptrs[i+1];
		for(ptr = sptrs[i]; ptr < eptr; ptr++) 
		{
			tids[tptrs[sids[ptr]]] = i;            
			tvals[tptrs[sids[ptr]]++] = svals[ptr]; 
		}
	}

	for(i = m; i > 0; i--) 
		tptrs[i] = tptrs[i-1];
	
	tptrs[0] = 0;
}

void ScaleOnePattern(	int *col_ptrs, int *col_ids, double *col_vals, int n, int m, int numIters)
{
	/* will overwrite col_vals
	 *  assume square matrix
	 */
	int i,j, it, eptr;
	double *rowSca, *colSca, sum;

	int *row_ptrs = (int *)malloc(sizeof(int) * (m+1));
	int *row_ids = (int *) malloc(sizeof(int) * col_ptrs[n] );

	symbtranspose(col_ptrs, col_ids, row_ptrs, row_ids,  n, m) ;

	rowSca = (double *) malloc(sizeof(double) * m);
	colSca = (double *) malloc(sizeof(double) * n);

	double art = ((double) m)/((double ) n);

	for (i  = 0; i < n; i++)
		colSca[i] = art;
	for (i  = 0; i < m; i++)
		rowSca[i] = 1.0;


	for(it = 0; it < numIters; it++)
	{

		for(j = 0; j < n; j++) 
		{
			sum = 0.0;
			if(col_ptrs[j+1] > col_ptrs[j])
			{
				for(eptr = col_ptrs[j]; eptr < col_ptrs[j+1]; eptr++) 
					sum += rowSca[col_ids[eptr]];

				colSca[j] = art / sum;
			}
		}
		for(i = 0; i < m; i++) 
		{
			sum = 0.0;
			if (row_ptrs[i+1] > row_ptrs[i])
			{
				for(eptr = row_ptrs[i]; eptr < row_ptrs[i+1]; eptr++) 
					sum += colSca[row_ids[eptr]];
				
				rowSca[i] = 1.0 / sum;
			}
		}

	}
	/*update*/
	for(j = 0; j < n; j++) 
	{
		double csca = colSca[j];	
		for(eptr = col_ptrs[j]; eptr < col_ptrs[j+1]; eptr++) 
			col_vals[eptr] = rowSca[col_ids[eptr]] * csca;

	}

	free(colSca);   free(rowSca);
	free(row_ids); free(row_ptrs);
}

void ScaleOne(int *col_ptrs, int *col_ids, double *col_vals, int n, int m, int numIters)
{
	/* will overwrite col_vals
	 *  assume square matrix
	 */
	int i,j, it, eptr;
	double *rowSca, *colSca, sum;
	int * row_ptrs = (int *)malloc(sizeof(int) * (m+1));
	int *row_ids = (int *) malloc(sizeof(int) * col_ptrs[n] );
	double * row_vals = (double *) malloc(sizeof(double) * col_ptrs[n]);

	transpose(col_ptrs, col_ids, col_vals, row_ptrs, row_ids,  row_vals, n, m) ;

	rowSca = (double *) malloc(sizeof(double) * m);
	colSca = (double *) malloc(sizeof(double) * n);

	double art = ((double) m)/((double ) n);

	for (i  = 0; i < n; i++)
		colSca[i] = art;
	for (i  = 0; i < m; i++)
		rowSca[i] = 1.0;


	for(it = 0; it < numIters; it++)
	{

		for(j = 0; j < n; j++) 
		{
			sum = 0.0;
			if(col_ptrs[j+1] > col_ptrs[j])
			{			
				for(eptr = col_ptrs[j]; eptr < col_ptrs[j+1]; eptr++) 
					sum += rowSca[col_ids[eptr]] * col_vals[eptr];

				colSca[j] = art / sum;
			}
		}
		for(i = 0; i < m; i++) 
		{
			sum = 0.0;
			if (row_ptrs[i+1] > row_ptrs[i])
			{
				for(eptr = row_ptrs[i]; eptr < row_ptrs[i+1]; eptr++) 
					sum += colSca[row_ids[eptr]] * row_vals[eptr];

				rowSca[i] = 1.0 / sum;
			}
		}

	}
	/*update*/
	for(j = 0; j < n; j++) 
	{
		double csca = colSca[j];	
		for(eptr = col_ptrs[j]; eptr < col_ptrs[j+1]; eptr++) 
			col_vals[eptr] *= rowSca[col_ids[eptr]] * csca ;

	}

	free(colSca);   free(rowSca);
	free(row_vals); free(row_ids); free(row_ptrs);
}

double computeScalingError(	int *col_ptrs, int *col_ids, double *col_vals, int n, int m)
{
	int i, j, eptr;
	double *rowSums = (double *) malloc(sizeof(double) * m);
	double err = -1.0, csum;
	for (i = 0; i < m; i++)
		rowSums[i] = 0.0;

	for(j = 0; j < n; j++) 
	{
		csum = 0.0;
		for(eptr = col_ptrs[j]; eptr < col_ptrs[j+1]; eptr++) 
		{
			rowSums[col_ids[eptr]] += col_vals[eptr];
			csum += col_vals[eptr];
		}
		if(fabs(1.0-csum) > err)
			err = fabs(1.0-csum) ;
	}

	for(i= 0; i < m; i++) 
	{
		if(fabs(1.0-rowSums[i]) > err)
			err = fabs(1.0-rowSums[i]);	
	}

	free(rowSums);
	return err;
}


void ReadMatrixMarket(char *fname, 
	double **A, int **csc_irn, int **csc_jcn, int *_m,
	int *_n, int *_ne, int readVals)
{
    /*START: MM stuff. */
	MM_typecode matcode;
    /*END: MM stuff. */

	FILE *ifp;
	double *val = NULL, *lA = NULL;
	int *I, *J;
	int *wrk; int *lcsc_jcn, *lcsc_irn;
	int i, hasValues;

	ifp = (FILE*) fopen(fname, "r");

	if(ifp == NULL)
	{
		printf("could not open file %s\n", fname);
		exit(12);
	}

	if (mm_read_banner(ifp, &matcode) != 0)
	{
		fclose(ifp);
		printf("Could not process Matrix Market banner.\n");
		exit(13);
	}
	if (! mm_is_matrix(matcode) )
	{
		fclose(ifp);
		printf("Sorry, this application does not support non matrices");
			exit(10);
	}
	if( !mm_is_sparse(matcode))
	{
		fclose(ifp);
		printf("Sorry, this application does not support dense matrices");
		exit(14);
	}
	if (mm_is_complex(matcode))
	{
		fclose(ifp);
		printf("In this application, the magnitudes of the complex values are used\n");
	//	exit(15);

	}
	/*read in everything.*/
	mm_read_unsymmetric_sparse_mod(fname, _m, _n, _ne, &val, &I, &J, &hasValues);
	if(hasValues == 0 && readVals != 0) 
	{
		printf("sorry no values in the file\n");
		fflush(stdout);
		exit(12);
	}
	printf("matrix is read in unsymm mod %d %d %d\n", *_m, *_n, *_ne);
	lcsc_jcn = *csc_jcn = (int *) malloc(sizeof(int) * ((*_n)+1));
	lcsc_irn = *csc_irn = (int *) malloc(sizeof(int) * (*_ne));

	printf("within this file %s we always allocate values\n", __FILE__ );
	lA = *A = (double *) malloc(sizeof(double) * (*_ne));

	wrk = (int *) malloc(sizeof(int ) * ((*_n)+1));
	memset(wrk, 0, sizeof(int)*((*_n)+1));
	for(i = 0; i < *_ne; i++)
	{
		if(J[i] < 0)
		{
			printf("col index <  0\n");
			exit(16);
		}
		if(J[i] >= *_n)
		{
			printf("col index >= n\n");
			exit(17);
		}
		wrk[J[i]]++;
	}

	lcsc_jcn[0] = 0;
	for(i = 1; i<= *_n; i++)
	{
		lcsc_jcn[i] = wrk[i-1];
		wrk[i] = wrk[i-1] + wrk[i];
	}

	if(wrk[*_n] != *_ne)
	{
		printf("numbers do not add up\n");
		exit(18);
	}

	for(i = *_ne-1; i>=0; i--)
	{
		int rn = I[i], cn = J[i];
		if(rn < 0)   
		{
			printf("row index <  0\n");
			exit(19);
		}
		if(rn >= *_m) 
		{
			printf("row index >=  m\n"); 
			exit(20);
		}  
		lcsc_irn[--(wrk[cn])] = rn;		
		if( readVals )
			lA[wrk[cn]] = val[i];

	}
	for(i=0; i<=*_n; i++)
	{
		if(lcsc_jcn[i] != wrk[i]) 
		{
			printf("arrays do not line up\n");
			exit(21);
		}
	}
	free(I);
	free(J);
	if( val )
		free(val);
	free(wrk);

	return ;
}

#ifdef MAIN_TEST_SCALING

int main(int argc, char *argv[])
{

	int j, nz;

	int m, n;
	int* col_ptrs;
	int* col_ids;
	int iters = 20, scalePattern = 1;
	double *col_vals;
	double t0, t1;
 	if (argc < 2 || argc > 4)
	{
		printf("%s needs at least two arguments filename, scalePattern, its\n", argv[0]);
		printf("\tdefult value for its = 20\n");

		exit(12);
	}

	if(argc == 4)
	{
		iters = atoi(argv[3]);	
		scalePattern = atoi(argv[2]);	
	}
	else if (argc == 3)
	{
		scalePattern = atoi(argv[2]);	
	}
	
	t0 = u_wseconds();
	printf("will call with %d (argc %d)\n", !scalePattern, argc);

	ReadMatrixMarket(argv[1], &col_vals, &col_ids, &col_ptrs, &m,	&n, &nz, !scalePattern);
	for(j = 0; j<nz; j++) 
	{
		col_vals[j] = col_vals[j] >=0 ? col_vals[j] : -col_vals[j];
	}

	t1 = u_wseconds();
	printf("read-in time %.2f\n", t1-t0);

	/*check input*/
	if(m != n)
	{
		printf("we expect a square matrix\n");
		exit(12);
	}
	
	/*do computation*/
	t0 = u_wseconds();
	if ( scalePattern ) 
	  ScaleOnePattern(	 col_ptrs,  col_ids, col_vals,  n, m, iters);
	else
		ScaleOne(col_ptrs, col_ids, col_vals, n, m, iters);
	t1 = u_wseconds();
	/*done computing*/


	printf("total scaling time %.2f err %.4f\n", t1-t0, computeScalingError(col_ptrs, col_ids, col_vals,  n, m));


	free(col_vals);
	free(col_ptrs);
	free(col_ids);

	return 0;
}

#endif
