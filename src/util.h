#ifndef _UTIL_H
#define _UTIL_H
class util {
public:
	static int* ivec(int x);
	static void free_ivec(int *x);
	static int** imat(int x, int y);
	static void free_imat(int **x);
	static double* dvec(int x);
	static void free_dvec(double *x);
	static double** dmat(int x, int y);
	static void free_dmat(double **x);

};
#endif