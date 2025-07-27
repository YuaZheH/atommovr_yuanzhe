#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#include "mmio.h"
#include "matrixUtils.h"

/*
We are interested JOB=2 and JOB 3. We expect, in general, JOB=2 to be faster.

 C   2 Compute a column permutation of the matrix so that the smallest 
 C     value on the diagonal of the permuted matrix is maximized.
This is augmenting path based.

 C   3 Compute a column permutation of the matrix so that the smallest
 C     value on the diagonal of the permuted matrix is maximized
 C     The algorithm differs from the one used for JOB = 2 and may
 C     have quite a different performance. 
This is threshold based.

*/

static double u_wseconds(void) 
{
  struct timeval tp;

  gettimeofday(&tp, NULL);

  return (double) tp.tv_sec + (double) tp.tv_usec / 1000000.0;
}

#define MAX_DOUBLE 1e100 


#define JOB_DEF 5


void  mc64id_(int *, double *);

void  mc64ad_(int *, int *, int *, int *, int *, int *,double *,
	      int *, int *, int *,int *, int *, double *,
	      int *, double *, int *);



int main(int argc, char *argv[])
{ 
  int i;
  int nicntl, ncntl, ninfo;
  int job, m, n, nz, num, liw, ldw;
  
  int *col_ptrs, *col_ids, *perm, *iw;
  int *icntl, *info;
  double *col_vals;
  double *dw, *cntl;
  double t0, t1;
  int optScaling, permute;

  if (argc != 5)
  {
    printf("need at least one additional arg <mtx filename> <JOB> <Scaling> <Permutation>\n");
    printf("\tScaling options (noscale 0, scale-pat 2, scale-A 3)\n");
    printf("\tPermutation (1,2,3, 4): no, col, row, row and col\n");
    exit(12);

  }

 optScaling = atoi(argv[3]);

  /* check of input arguments */
  ReadMatrixMarket(argv[1], &col_vals, &col_ids, &col_ptrs, &m, &n, &nz, !(optScaling==2));
  
  /* fix mc64 control parameters */
  nicntl=10;
  ncntl=10;
  ninfo=10;
  icntl = calloc(nicntl,sizeof(int));
  info =  calloc(ninfo,sizeof(int));
  cntl = calloc(ninfo,sizeof(double));

  mc64id_(icntl,cntl);

  icntl[0] = -1;
  icntl[1] = -1;


  job = atoi(argv[2]);

  /*allocation of mc64 workspace*/
  if(job == 1){
    liw = 4*n+m;
    ldw =  n+2*m+nz;
  }else if(job == 2 ){
    liw = 2*n+2*m;
    ldw = m;
  }else if(job == 3 ){
    liw = 8*n+2*m+nz;
    ldw = nz;
  }else if(job == 4 ){
    liw = 3*n+2*m;
    ldw = 2*m+nz;
  }else if(job == 5 ){
    liw = 3*n+2*m;
    ldw = n+2*m+nz;
  }else if(job == 6 ){
    liw = 3*n+2*m+nz;
    ldw = n+3*m+nz;
  }else{
    free(icntl);
    free(info); 
    free(cntl); free(col_vals); free(col_ids); free(col_ptrs);
    printf ("Wrong job value\n");
    exit(12);
  } 
  iw = calloc(liw,sizeof(int));
  if(ldw > 0)  
    dw = calloc(ldw,sizeof(double));


  for(i=0;i<nz;i++)
    col_vals[i] = col_vals[i]>0 ? col_vals[i] : -col_vals[i];
  
 optScaling = atoi(argv[3]);

  if(optScaling == 2)
  {
    printf("\tpattern scaling\n");
    ScaleOnePattern(   col_ptrs,  col_ids, col_vals,  n, m, 20);
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
  
  permute = atoi(argv[4]);
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

  for(i=0;i<n+1;i++)
    col_ptrs[i] = col_ptrs[i]+1;
  
  for(i=0;i<nz;i++)
    col_ids[i] = col_ids[i]+1;

    
  perm = calloc(m,sizeof(int));

  t0 = u_wseconds();
  mc64ad_(&job,&m,&n,&nz,col_ptrs,col_ids,col_vals,&num,perm,&liw,iw,&ldw,dw,icntl,cntl,info);
  t1 = u_wseconds();
  printf("total mc64 time %.2f\n", t1-t0);

  free(icntl);
  free(info); 
  free(cntl);
  free(iw);
  if(ldw > 0)
    free(dw);
  free(col_ptrs);
  free(col_ids);
  free(perm);
  free(col_vals);
} 
