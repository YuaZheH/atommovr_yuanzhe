/*
 * File: cheap.c
 * Content: Contains jump-start heuristics
 * Author: Kamer Kaya, Johannes Langguth and Bora Ucar
 * Version: 1.0
 *
 * Please see the papers 
 * @article{klmu:13m,
 *        Author = {Kaya, Kamer and Langguth, Johannes and Manne, Fredrik and U\c{c}ar, Bora},
 *        Journal = {Computers \& Operations Research},
 *        Number = {5},
 *        Pages = {1266--1275},
 *        Title = {Push-relabel based algorithms for the maximum transversal problem},
 *        Volume = {40},
 *        URL = {https://hal.inria.fr/hal-00763920},
 *        Year = {2013}
 *}
 *
 *@article{duku:12m,
 *        Author = {Duff, Iain S. and Kaya, Kamer and U\c{c}ar, Bora},
 *        Journal = {{ACM} Transactions on Mathematical Software},
 *        Pages = {13:1--13:31},
 *        Title = {Design, implementation, and analysis of maximum transversal algorithms},
 *        Volume = {38},
 *        URL={https://hal.inria.fr/hal-00786548},
 *        Year = {2012}
 * }
 *
 *@inproceedings{pauc:20,
 *		  author = {Ioannis Panagiotas and Bora U{\c{c}}ar},
 *		  booktitle = {28th Annual European Symposium on Algorithms, {ESA} 2020, September 7-9, 2020, Pisa, Italy (Virtual Conference)},
 *		  editor = {Fabrizio Grandoni and Grzegorz Herman and Peter Sanders},
 *	      pages = {76:1--76:23},
 *	      publisher = {Schloss Dagstuhl - Leibniz-Zentrum f{\"{u}}r Informatik},
 *	      series = {LIPIcs},
 *	      title = {Engineering Fast Almost Optimal Algorithms for Bipartite Graph Matching},
 *	      volume = {173},
 *	      year = {2020}}
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

#include "matchmaker.h"

struct node {int id; int degree; struct node *next; struct node* prvs;};
typedef struct node Node;

void old_cheap(int* col_ptrs, int* col_ids, int* match, int* row_match, int n, int m) {
	int ptr;
	int i = 0;
	for(; i < n; i++) {
		int s_ptr = col_ptrs[i];
		int e_ptr = col_ptrs[i + 1];
		for(ptr = s_ptr; ptr < e_ptr; ptr++) {
			int r_id = col_ids[ptr];
			if(row_match[r_id] == -1) {
				match[i] = r_id;
				row_match[r_id] = i;
				break;
			}
		}	
	}
}

void sk_cheap(int* col_ptrs, int* col_ids, int* row_ptrs, int* row_ids,
		int* match, int* row_match, int n, int m){
	int i;

	int* col_stack = (int*)malloc(n * sizeof(int));
	int* col_degrees = (int*)malloc(n * sizeof(int));
	memset(col_degrees, 0, n * sizeof(int));

	int no_of_d1_cols = 0;
	for(i = 0; i < n; i++) {
		col_degrees[i] = col_ptrs[i+1] - col_ptrs[i];
		if(col_degrees[i] == 1) {
			col_stack[no_of_d1_cols++] = i;
		}
	}

	int* row_stack = (int*)malloc(m * sizeof(int));
	int* row_degrees = (int*)malloc(m * sizeof(int));
	memset(row_degrees, 0, m * sizeof(int));

	int no_of_d1_rows = 0;
	for(i = 0; i < m; i++) {
		row_degrees[i] = row_ptrs[i+1] - row_ptrs[i];
		if(row_degrees[i] == 1) {
			row_stack[no_of_d1_rows++] = i;
		}
	}

	int stop = 0;
	int r_id = -1, c_id, r_id2, c_id2;
	int sptr, eptr, ptr;
	int sptr2, eptr2, ptr2;

	int remain = 0;
	int c_degree = 0;

	while(!stop) {
		while(no_of_d1_rows > 0 || no_of_d1_cols > 0) {
			if(no_of_d1_rows > 0) {
				r_id = row_stack[--no_of_d1_rows];
				if(row_degrees[r_id] == 1 && row_match[r_id] == -1) {
					sptr = row_ptrs[r_id];
					eptr = row_ptrs[r_id + 1];
					for(ptr = sptr; ptr < eptr; ptr++) {
						c_id = row_ids[ptr];
						if(match[c_id] == -1) {
							match[c_id] = r_id;
							row_match[r_id] = c_id;

							sptr2 = col_ptrs[c_id];
							eptr2 = col_ptrs[c_id + 1];
							for(ptr2 = sptr2; ptr2 < eptr2; ptr2++) {
								r_id2 = col_ids[ptr2];
								if(row_match[r_id2] == -1) {
									if((--(row_degrees[r_id2])) == 1) {
										row_stack[no_of_d1_rows++] = r_id2;
									}
								}
							}
							break;
						}
					}
				}
			}

			if(no_of_d1_cols > 0) {
				c_id = col_stack[--no_of_d1_cols];
				if(col_degrees[c_id] == 1 && match[c_id] == -1) {
					sptr = col_ptrs[c_id];
					eptr = col_ptrs[c_id + 1];
					for(ptr = sptr; ptr < eptr; ptr++) {
						r_id = col_ids[ptr];
						if(row_match[r_id] == -1) {
							row_match[r_id] = c_id;
							match[c_id] = r_id;

							sptr2 = row_ptrs[r_id];
							eptr2 = row_ptrs[r_id + 1];
							for(ptr2 = sptr2; ptr2 < eptr2; ptr2++) {
								c_id2 = row_ids[ptr2];
								if( match[c_id2] == -1) {
									if((--(col_degrees[c_id2])) == 1) {
										col_stack[no_of_d1_cols++] = c_id2;
									}
								}
							}
							break;
						}
					}
				}
			}
		}

		stop = 1;
		for(i = remain; i < n; i++) {
			c_id = i;
			c_degree = col_degrees[c_id];

			if(match[c_id] == -1 && c_degree != 0) {
				sptr = col_ptrs[c_id];
				eptr = col_ptrs[c_id + 1];

				for(ptr = sptr; ptr < eptr; ptr++) {
					r_id = col_ids[ptr];
					if(row_match[r_id] == -1) {
						match[c_id] = r_id;
						row_match[r_id] = c_id;
						stop = 0;
						break;
					}
				}
				ptr++;

				for(;ptr < eptr; ptr++) {
					r_id2 = col_ids[ptr];
					if(row_match[r_id2] == -1) {
						if((--(row_degrees[r_id2])) == 1) {
							row_stack[no_of_d1_rows++] = r_id2;
						}
					}
				}

				sptr = row_ptrs[r_id];
				eptr = row_ptrs[r_id + 1];
				int count = row_degrees[r_id];
				for(ptr = sptr;ptr < eptr && count > 0; ptr++) {
					c_id2 = row_ids[ptr];
					if(match[c_id2] == -1) {
						count--;
						if((--(col_degrees[c_id2])) == 1) {
							col_stack[no_of_d1_cols++] = c_id2;
						}
					}
				}
			}

			if(no_of_d1_cols + no_of_d1_rows > 0) {
				remain = i + 1;
				break;
			}

			if(i == n-1) {
				stop = 1;
			}
		}
	}

	free(row_degrees);
	free(row_stack);
	free(col_degrees);
	free(col_stack);
}

void sk_cheap_rand(int* col_ptrs, int* col_ids, int* row_ptrs, int* row_ids,
		int* match, int* row_match, int n, int m) {
	int i;

	int* col_stack = (int*)malloc(n * sizeof(int));
	int* col_degrees = (int*)malloc(n * sizeof(int));
	memset(col_degrees, 0, n * sizeof(int));

	int no_of_d1_cols = 0;
	for(i = 0; i < n; i++) {
		col_degrees[i] = col_ptrs[i+1] - col_ptrs[i];
		if(col_degrees[i] == 1) {
			col_stack[no_of_d1_cols++] = i;
		}
	}

	int* row_stack = (int*)malloc(m * sizeof(int));
	int* row_degrees = (int*)malloc(m * sizeof(int));
	memset(row_degrees, 0, m * sizeof(int));

	int no_of_d1_rows = 0;
	for(i = 0; i < m; i++) {
		row_degrees[i] = row_ptrs[i+1] - row_ptrs[i];
		if(row_degrees[i] == 1) {
			row_stack[no_of_d1_rows++] = i;
		}
	}

	int* randarr = (int*)malloc(n * sizeof(int));
	for(i = 0; i < n; i++){randarr[i] = i;}
	int temp;
	for(i = n-1; i >= 0; i--) {
		int z = rand() % (i+1);
		temp = randarr[i]; randarr[i] = randarr[z]; randarr[z] = temp;
	}

	int stop = 0;
	int r_id = -1, c_id, r_id2, c_id2, e_id;
	int sptr, eptr, ptr;
	int sptr2, eptr2, ptr2;

	int remain = 0;
	int c_degree = 0;

	while(!stop) {
		while(no_of_d1_rows > 0 || no_of_d1_cols > 0) {
			if(no_of_d1_rows > 0) {
				r_id = row_stack[--no_of_d1_rows];
				if(row_degrees[r_id] == 1 && row_match[r_id] == -1) {
					sptr = row_ptrs[r_id];
					eptr = row_ptrs[r_id + 1];
					for(ptr = sptr; ptr < eptr; ptr++) {
						c_id = row_ids[ptr];
						if(match[c_id] == -1) {
							match[c_id] = r_id;
							row_match[r_id] = c_id;

							sptr2 = col_ptrs[c_id];
							eptr2 = col_ptrs[c_id + 1];
							for(ptr2 = sptr2; ptr2 < eptr2; ptr2++) {
								r_id2 = col_ids[ptr2];
								if(row_match[r_id2] == -1) {
									if((--(row_degrees[r_id2])) == 1) {
										row_stack[no_of_d1_rows++] = r_id2;
									}
								}
							}
							break;
						}
					}
				}
			}

			if(no_of_d1_cols > 0) {
				c_id = col_stack[--no_of_d1_cols];
				if(col_degrees[c_id] == 1 && match[c_id] == -1) {
					sptr = col_ptrs[c_id];
					eptr = col_ptrs[c_id + 1];
					for(ptr = sptr; ptr < eptr; ptr++) {
						r_id = col_ids[ptr];
						if(row_match[r_id] == -1) {
							row_match[r_id] = c_id;
							match[c_id] = r_id;

							sptr2 = row_ptrs[r_id];
							eptr2 = row_ptrs[r_id + 1];
							for(ptr2 = sptr2; ptr2 < eptr2; ptr2++) {
								c_id2 = row_ids[ptr2];
								if( match[c_id2] == -1) {
									if((--(col_degrees[c_id2])) == 1) {
										col_stack[no_of_d1_cols++] = c_id2;
									}
								}
							}
							break;
						}
					}
				}
			}
		}

		stop = 1;
		for(i = remain; i < n; i++) {
			c_id = randarr[i];
			c_degree = col_degrees[c_id];

			if(match[c_id] == -1 && c_degree != 0) {
				e_id = rand() % c_degree;

				sptr = col_ptrs[c_id];
				eptr = col_ptrs[c_id + 1];

				for(ptr = sptr; ptr < eptr; ptr++) {
					r_id = col_ids[ptr];
					if(row_match[r_id] == -1) {
						if(e_id == 0) {
							match[c_id] = r_id;
							row_match[r_id] = c_id;
							stop = 0;
							break;
						} else {
							if((--(row_degrees[r_id])) == 1) {
								row_stack[no_of_d1_rows++] = r_id;
							}
							e_id--;
						}
					}
				}
				ptr++;

				for(;ptr < eptr; ptr++) {
					r_id2 = col_ids[ptr];
					if(row_match[r_id2] == -1) {
						if((--(row_degrees[r_id2])) == 1) {
							row_stack[no_of_d1_rows++] = r_id2;
						}
					}
				}

				sptr = row_ptrs[r_id];
				eptr = row_ptrs[r_id + 1];
				int count = row_degrees[r_id];
				for(ptr = sptr;ptr < eptr && count > 0; ptr++) {
					c_id2 = row_ids[ptr];
					if(match[c_id2] == -1) {
						count--;
						if((--(col_degrees[c_id2])) == 1) {
							col_stack[no_of_d1_cols++] = c_id2;
						}
					}
				}
			}

			if(no_of_d1_cols + no_of_d1_rows > 0) {
				remain = i + 1;
				break;
			}

			if(i == n-1) {
				stop = 1;
			}
		}
	}

	free(randarr);
	free(row_degrees);
	free(row_stack);
	free(col_degrees);
	free(col_stack);
}


void mind_cheap(int* col_ptrs, int* col_ids, int* row_ptrs, int* row_ids,
		int* match, int* row_match, int n, int m) {

	Node* rnodes = (Node*)malloc(sizeof(Node) * m);
	Node* cnodes = (Node*)malloc(sizeof(Node) * n);;
	Node* tptr;

	int i, deg, maxdeg = -1, cdeg, vtx, minnbr = -1, ptr, row ,col, temp;

	for(i = 0; i < n; i++) {
		deg = col_ptrs[i+1] - col_ptrs[i];
		cnodes[i].degree = deg;
		cnodes[i].id = i;
		if(deg > maxdeg) maxdeg = deg;
	}

	for(i = 0; i < m; i++) {
		deg = row_ptrs[i+1] - row_ptrs[i];
		rnodes[i].degree = deg;
		rnodes[i].id = i + n;
		if(deg > maxdeg) maxdeg = deg;
	}

	Node* lists = (Node*)malloc(sizeof(Node) * (maxdeg + 1));
	Node* listse = (Node*)malloc(sizeof(Node) * (maxdeg + 1));

	for(i = 0; i <= maxdeg; i++) {
		lists[i].next = &(listse[i]); lists[i].prvs = (Node*)0;
		listse[i].next = (Node*)0; listse[i].prvs = &(lists[i]);
		lists[i].id = -1; listse[i].id = -1;
		lists[i].degree = i; listse[i].degree = i;
	}

	for(i = 0; i < n; i++) {
		deg = cnodes[i].degree;
		if(deg > 0) {
			tptr = lists[deg].next;
			tptr->prvs = lists[deg].next = &(cnodes[i]);
			cnodes[i].next = tptr;
			cnodes[i].prvs = &(lists[deg]);
		}
	}
	for(i = 0; i < m; i++) {
		deg = rnodes[i].degree;
		if(deg > 0) {
			tptr = lists[deg].next;
			tptr->prvs = lists[deg].next = &(rnodes[i]);
			rnodes[i].next = tptr;
			rnodes[i].prvs = &(lists[deg]);
		}
	}

	cdeg = 1;
	while(cdeg <= maxdeg) {
		if(lists[cdeg].next == &(listse[cdeg])) {cdeg++; continue;}
		tptr = lists[cdeg].next;
		lists[cdeg].next = tptr->next;
		tptr->next->prvs = &(lists[cdeg]);
		vtx = tptr->id;

		if(vtx < n) {
			for(ptr = col_ptrs[vtx]; ptr < col_ptrs[vtx+1]; ptr++) {
				if(row_match[col_ids[ptr]] == -1) {
					minnbr = col_ids[ptr];
					break;
				}
			}

			for(ptr = ptr + 1; ptr < col_ptrs[vtx+1]; ptr++) {
				row = col_ids[ptr];
				if(row_match[row] == -1) {
					if(rnodes[row].degree < rnodes[minnbr].degree) {
						minnbr = col_ids[ptr];
					}
				}
			}

			match[vtx] = minnbr; row_match[minnbr] = vtx;
			rnodes[minnbr].next->prvs = rnodes[minnbr].prvs;
			rnodes[minnbr].prvs->next = rnodes[minnbr].next;
		} else {
			vtx = vtx - n;
			for(ptr = row_ptrs[vtx]; ptr < row_ptrs[vtx+1]; ptr++) {
				if(match[row_ids[ptr]] == -1) {
					minnbr = row_ids[ptr];
					break;
				}
			}

			for(ptr = ptr + 1; ptr < row_ptrs[vtx+1]; ptr++) {
				col = row_ids[ptr];
				if(match[col] == -1) {
					if(cnodes[col].degree < cnodes[minnbr].degree) {
						minnbr = row_ids[ptr];
					}
				}
			}

			row_match[vtx] = minnbr; match[minnbr] = vtx;
			cnodes[minnbr].next->prvs = cnodes[minnbr].prvs;
			cnodes[minnbr].prvs->next = cnodes[minnbr].next;
			temp = vtx; vtx = minnbr; minnbr = temp; /* swap */
		}

		for(ptr = col_ptrs[vtx]; ptr < col_ptrs[vtx+1]; ptr++) {
			row = col_ids[ptr];
			if(row_match[row] == -1) {
				deg = --(rnodes[row].degree);
				rnodes[row].next->prvs = rnodes[row].prvs;
				rnodes[row].prvs->next = rnodes[row].next;

				if(deg > 0) {
					tptr = lists[deg].next;
					tptr->prvs = lists[deg].next = &(rnodes[row]);
					rnodes[row].next = tptr;
					rnodes[row].prvs = &(lists[deg]);
				}
			}
		}

		for(ptr = row_ptrs[minnbr]; ptr < row_ptrs[minnbr+1]; ptr++) {
			col = row_ids[ptr];
			if(match[col] == -1) {
				deg = --(cnodes[col].degree);
				cnodes[col].next->prvs = cnodes[col].prvs;
				cnodes[col].prvs->next = cnodes[col].next;

				if(deg > 0) {
					tptr = lists[deg].next;
					tptr->prvs = lists[deg].next = &(cnodes[col]);
					cnodes[col].next = tptr;
					cnodes[col].prvs = &(lists[deg]);
				}
			}
		}
		cdeg--;
	}

	free(listse);
	free(lists);
	free(cnodes);
	free(rnodes);
}

void getSinkhornKnoppScalingOfPattern(int* col_ptrs, int* col_ids, int* row_ptrs, int* row_ids, double *rowSca, double *colSca, int m, int n, int numIters)
{
/*we assume m rows and n columns
PRE: rowSca of size m is allocated
	 colSca of size n is alllocated
*/

	int i, j, it, eptr;
	double sum;

	double noverm = ((double) n) / ((double) m);

	for(i = 0; i < m; i++) 
		rowSca[i] = noverm;

	for(j = 0; j < n; j++) 
		colSca[j] = 1.0;

	for(it = 0; it < numIters; it++)
	{

		for(i = 0; i < m; i++) 
		{
			sum = 0.0;
			for(eptr = row_ptrs[i]; eptr < row_ptrs[i+1]; eptr++) 
				sum += colSca[row_ids[eptr]];
			
			rowSca[i] = noverm / sum;
		}

		for(j = 0; j < n; j++) 
		{
			sum = 0.0;
			for(eptr = col_ptrs[j]; eptr < col_ptrs[j+1]; eptr++) 
				sum += rowSca[col_ids[eptr]];
			
			colSca[j] = 1.0 / sum;
		}
	}
	
	return ;

}

/****************************************************************************
*Truncated walk version.
*****************************************************************************/
#define randv(maxv) (maxv) * (rand()/(double) RAND_MAX)

#define matchingEdgeVal(vals, offset) offset == 0 ? vals[offset] : 	vals[offset] - vals[offset - 1]

void createVisitOrder(int *visitOrder, int n)
{
	int i, tmp;
	for(i = 0; i < n; i++) 
		visitOrder[i] = i;

	for(i = n-1; i >= 0; i--) 
	{
		int z = rand() % (i+1);
		tmp = visitOrder[i]; 
		visitOrder[i] = visitOrder[z]; 
		visitOrder[z] = tmp;
	}
}


static inline void prefixSum(double *vals_col, int sz)
{
	int i;
	for (i = 1; i < sz; i++)
		vals_col[i] += vals_col[i-1];
}

static inline int binSearch(double *vals_col, int sz, double rval)
{
	int l = 0;
	int r = sz-1;

	while (l <= r)
	{
		int m = l + (r-l)/2;
		if(l == r-1)
		{
			if (vals_col [l] < rval)
				return r;
			else 
				return l;
		}
		else 
		{
			if (l == r)
				return l;
		}
		if(vals_col [m] < rval)
			l = m +1 ;
		else				
			r = m ;
	}
	return -1;/*in case the list was empty or we could not find*/
}

static inline int sampleLognTime(int *rids_col, double *vals_col, int sz, int col, int itsmate, int offset, int *rowmatch, int *lookAhead, int base, int LAend )
{
	double medgeVal = 0.0, rsum, rval;
	int newoffset = -1, rid;

	while(newoffset == -1 && lookAhead[col] < LAend)/*look ahead mechanism*/
	{
		rid = rids_col[lookAhead[col] - base];	
		if(rowmatch[rid] == -1)
			newoffset = lookAhead[col] - base;		
		lookAhead[col]++;
	}

	if(newoffset == -1)/*look a head did not find anything*/
	{
		if (itsmate != -1)	
			medgeVal = matchingEdgeVal(vals_col, offset); 

		rsum = vals_col[sz - 1] - medgeVal;
		rval = randv(rsum);         /*a random val*/

		if (itsmate == -1)/*		binary search in the whole list*/ 
			newoffset = binSearch(vals_col, sz, rval);	
		else
		{
			if(vals_col [offset] - medgeVal >= rval || offset == sz - 1)/*search 0 -- offset -1*/
			newoffset = binSearch(vals_col, offset , rval);
			else/*search vals_col[offset +1 ...sz], the value (rval + medgeVal)*/	
			newoffset = binSearch(&(vals_col[offset + 1]), sz-offset -1, rval + medgeVal) + offset + 1;	
		}
	}	

	return newoffset;

}
void goelKapralovKhannaTruncated(int *cptrs, int *rids, double *vals, int *match, int *rowmatch, int n, int m)
{
	/*we have m=n in general, but keeping it clean.*/
	/*we are going to sample-out-edge from cols.*/
	int i, j, v, offset, pathLast, col, pcol, row, rtmp;
	int *path, *rowFromWhere, *visitOrder, *medgeOffsets, *newMedgeOffsets, *lookAhead;
	int walkLength, numTrials, cardMatch; 
	double totalLength;

	rowFromWhere = (int *) malloc(sizeof(int) * m);/*if a row is already seen in a randomWalk, its position is kept here*/
	path  = (int *) malloc(sizeof(int) * n);       /*the randomWalk, where the cycles are discarded, so a path*/
	visitOrder = (int *) malloc(sizeof(int) * n);
	medgeOffsets = (int *) malloc(sizeof(int) * n);
	newMedgeOffsets = (int *) malloc(sizeof(int) * n);
	lookAhead = (int *) malloc(sizeof(int) * n);

	srand(time(NULL));

	createVisitOrder(visitOrder, n);
	for (i = 0; i < m; i++)	
		rowmatch[i] = rowFromWhere[i] = -1;

	for (j = 0; j < n; j++)
	{
		match[j] = 	medgeOffsets[j]= -1;
		prefixSum(&(vals[cptrs[j]]), cptrs[j+1] - cptrs[j]);
		lookAhead[j] = cptrs[j];
	}

	cardMatch = 0;
	totalLength = 0.0;
	for (v = 0; v < n; v ++)/*we will do n augmentations from the empty matching*/
	{
		j = visitOrder[v];

		if(cptrs[j+1] - cptrs[j] == 0)/*skip empty columns*/
			continue;
		if(cptrs[j+1] - cptrs[j] == 1 & match[j] != -1)
			continue;
		numTrials = 0;
		while(match[j] == -1 && numTrials < 1)/*success with probability 1-1/2^10, about 0.999*/
		{
			walkLength = 0;
			path[0] = j; 
			pathLast = 0;
			numTrials ++;
			while (pathLast >= 0 && walkLength <= (int) 2 *(4 + ceil ( (2.0 * n)/(n-v))))
			{
				col = path[pathLast];
				/*pick a random unmatched row from this col*/
				offset = sampleLognTime(&(rids[cptrs[col]]), &(vals[cptrs[col]]), cptrs[col+1] - cptrs[col], col,  match[col], medgeOffsets[col], rowmatch, lookAhead, cptrs[col], cptrs[col+1]);
	
				if(offset == -1) /*Guarding against vertices with degreees 0 or 1s*/							
					break;	
							
				row = rids[cptrs[col] + offset];
				newMedgeOffsets[col] = offset;

				walkLength ++;
				totalLength += 1.0;
				if(rowmatch[row] != - 1)/*not an augmenting path*/
				{	/*grow the path*/
					if(rowFromWhere[row] == -1 )/* this row is new for randomWalk at j*/
					{
						rowFromWhere[row] = pathLast+1;
						path[++pathLast] = rowmatch[row];
						walkLength++; 
					}
					else/*we have already seen this row in this randomWalk, discard the cycle*/
					{   /*from pathLast down to romFromWhere[row], we should discard the vertices; they form the cycle*/
						while(pathLast > rowFromWhere[row])
							rowFromWhere[match[path[pathLast--]]] = -1;
					}
				}	
				else/*found an augmenting path*/
				{
					cardMatch++;
					for ( ; pathLast >= 0; pathLast--)
					{
						pcol = path[pathLast];
						rtmp = match[pcol];
						medgeOffsets[pcol] = newMedgeOffsets[pcol];
						match[pcol] = row; 
						rowFromWhere[row] = -1;
						rowmatch[row] = pcol;
						row = rtmp;					
					}
					break;
				}
			}/*..."while" randomly walking*/
			for ( ; pathLast >= 0; pathLast--)/*clean up rowFromWhere*/
			{
				pcol = path[pathLast];
				row = match[pcol];	
				if(row != -1)
					rowFromWhere[row] = -1;
			}
		}
	}
	free(lookAhead);
	free(newMedgeOffsets);
	free(medgeOffsets);
	free(visitOrder);
	free(path);
	free(rowFromWhere);
}

void truncrw(int* col_ptrs, int* col_ids, int* row_ptrs, int* row_ids,
	int* match, int* row_match, int n, int m)
{
	double *vals, *rowSca, *colSca;
	int j, eptr, nnz = col_ptrs[n];
	double scaj ;
	
	rowSca = (double *) malloc(sizeof(double) * m);
	colSca = (double *) malloc(sizeof(double) * n);

	getSinkhornKnoppScalingOfPattern(col_ptrs, col_ids, row_ptrs, row_ids, rowSca, colSca, m, n, 5);

	/*apply scaling to the pattern*/
	vals = (double *) malloc(sizeof(double) * nnz);
	for(j = 0; j < n; j++) 
	{
		scaj = colSca[j];

		for(eptr = col_ptrs[j]; eptr < col_ptrs[j+1]; eptr++) 
			vals[eptr] = scaj * rowSca[col_ids[eptr]];
	}

	goelKapralovKhannaTruncated(col_ptrs, col_ids, vals, match, row_match, n, m);

	free(vals);
	free(colSca);
	free(rowSca);
}

void cheap_matching(int* col_ptrs, int* col_ids, int* row_ptrs, int* row_ids,
		int* match, int* row_match, int n, int m, int cheap_id) {
	/*All these initialization heuristics assume to start from an empty matching.*/

	if(cheap_id == do_old_cheap) {
		old_cheap(col_ptrs, col_ids, match, row_match, n, m);
	} else if(cheap_id == do_sk_cheap) {
		sk_cheap(col_ptrs, col_ids, row_ptrs, row_ids, match, row_match, n, m);
	} else if(cheap_id == do_sk_cheap_rand) {
		sk_cheap_rand(col_ptrs, col_ids, row_ptrs, row_ids, match, row_match, n, m);
	} else if(cheap_id == do_mind_cheap) {
		mind_cheap(col_ptrs, col_ids, row_ptrs, row_ids, match, row_match, n, m);
	} else if (cheap_id == do_truncrw) {
		truncrw(col_ptrs, col_ids, row_ptrs,  row_ids, match, row_match, n, m);	
	}

}
