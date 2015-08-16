/*---------------------------------------------------
 * file:    topiclib.h
 * purpose: header file for topiclib.c
 * usage:   include in *.c
 * version: 1.0
 * author:  newman@uci.edu
 * date:    12/20/05
 *-------------------------------------------------*/

#ifndef _TOPICLIB_H_
#define _TOPICLIB_H_
 
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <time.h>
#include <iostream>
#include "MersenneTwister.h"
using namespace std;
typedef unsigned long uint32;

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SQR(a)   ((a) * (a))
#define CUBE(a)  ((a) * (a) * (a))

#define irand()  randomMT()
#define drand() (randomMT()/4294967296.0)
struct doc_author{
        int count;
        int* authors;
};
static MTRand mtrand;
void seedMT(uint32 seed);
inline uint32 randomMT(void);
void randomassignment(int ntot, int D, int T, int *w, int *d, int *c, int *z, int **wp, int **dp, int **cp, int *ztot, int *cztot);

int **imat(int nr, int nc); //
void free_imat(int **x); //
double **dmat(int nr, int nc); //
void free_dmat(double **x); //
float **fmat(int nr, int nc); //
void free_fmat(float **x); //
int *ivec(int n); //
void free_ivec(int *x);
double *dvec(int n); //
float *fvec(int n); //
int imax(int n, int *x); //
double dmax(int n, double *x); //
int imin(int n, int *x); //
double dmin(int n, double *x); //
int isum(int n, int *x); //
double dsum(int n, double *x); //
double ddot(int n, double *x, double *y); //
int countlines(char *fname); //
int countntot(char *fname); //
void isort(int n, int *x, int direction, int *indx); //
int *isort2(int n, int *x); //
int *dsort(int n, double *x); //
void read_docwordcountbin(int nnz, int *w, int *d, int *c, char *fname); //
int countnnz(int nr, int nc, int **x); //
void write_sparse(int nr, int nc, int **x, char *fname); //
void write_sparsebin(int nr, int nc, int **x, char *fname); //
int **read_sparse(char *fname, int *nr_, int *nc_); //
void read_dw(char *fname, int *d, int *w, int *c, int *D, int *W); //
void read_dwc(char *fname, int *d, int *w, int *c, int *D, int *W, int *C); //
void fill_dtot(int ntot, int *d, int *dtot); //
//void read_dwc(char *fname, int *d, int *w, int *c, int *D, int *W); //
int read_nnzbin(char *fname); //
int read_nnz(char *fname); //
void read_sparsebin(char *fname, int *col1, int *col2, int *col3); //
void write_ivec (int n, int *x, char *fname); //
void read_ivec (int n, int *x, char *fname); //
int *randperm(int n); //
//void randomassignment(int ntot, int T, int *w, int *d, int *z, int **wp, int **dp, int *ztot); //
//void randomassignment_coin(int ntot, int D, int T, int *w, int *d, int *c, int *z, int *s,int **wp, int **dp, int **cp, int** doc_cite_0, int** doc_cite_1, int *ztot, int *cztot);
void randomassignment_act(int ntot, int D, int T, int *w, int *d, int *c, int *z, int *s,
    int **wp, int **dp, int **cp, int** doc_cite_0, int** doc_cite_1, int *ztot, int *cztot, doc_author* map_doc_author, int* to_sample, int* to_test, int percent);
void randomassignment_atm(int ntot, int D, int T, int *w, int *d, int *c, int *z, int *s,
    int **wp, int **dp, int **cp, int** doc_cite_0, int** doc_cite_1, int *ztot, int *cztot, doc_author* map_doc_author, int* to_sample, int* to_test);
void randomassignment_coin_tu_etal(int ntot, int D, int T, int *w, int *d, int *c,
		int *z, int *s, int **wp, int **dp, int **cp, int** doc_cite_0,
		int** doc_cite_1, int *ztot, int *cztot, doc_author* map_doc_author,  doc_author* map_doc_cited_author, 
		int* to_sample, int* to_test);
void randomassignment_alt(int ntot, int D, int T, int *w, int *d, int *c, int *z, int *s,
    int **wp, int **dp, int **cp, int** doc_cite_0, int** doc_cite_1, int *ztot, int *cztot, doc_author* map_doc_author, int* to_sample, int* to_test);
void assignment(int ntot, int *w, int *d, int *z, int **wp, int **dp, int *ztot); //
void randomassignment_coin_sact(int ntot, int D, int T, int *w, int *d, int *c, int *z, int *s,
    int **wp, int **dp, int **cp, int** doc_cite_0, int** doc_cite_1, int *ztot, int *cztot, doc_author* map_doc_author, int* to_sample, int* to_test);
void randomassignment_2layer(int ntot, int T, int S, int *w, int *d, int *z, int *y, int **wp, int **zy, int **dp, int *ztot, int *ytot); //
void randomassignment2(int ntot, int T, int *d, int *z, int **dp); //
void sample_chain (int ntot, int W, int T, double alpha, double beta, int *w, int *d, int *z, int **wp, int **dp, int *ztot, int *order); //
void sample_chain_rank (int ntot, int W, int T, double alpha, double beta, int *w, int *d, int *drank, int *z, int **wp, int **dp, int *ztot, int *order); //
void sample_chain0 (int ntot, int W, int T, double alpha, double beta, int *w, int *d, int *z, int **wp, int **dp, int *ztot); //
void sample_chain_2layer (int ntot, int W, int T, int S, double alpha, double beta, double gamma, int *w, int *d, int *z, int *y, int **wp, int **zy, int **dp, int *ztot, int *ytot); //
void resample_chain (int ntot, int W, int T, double alpha, double beta, int *w, int *d, int *z, int **wp, int **dp, int *ztot); //
void oversample_dp (int ntot, int W, int T, double alpha, double beta, int *w, int *d, int *z, int **wp, int **dp, int *ztot); //
void loglike (int ntot, int W, int D, int T, double alpha, double beta, int *w, int *d, int **wp, int **dp, int *ztot, int *dtot); //
void chksum (int n, int T, int **x, int *sumx); //
void getztot (int n, int T, int **x, int *ztot); //
double etime(); //
void quickSort(double *numbers, int *index,int array_size);
void q_sort(double *numbers, int *index, int left, int right);
#endif
