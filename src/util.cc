
#include "util.h"

int* util::ivec(int x){
	int *ptr = new int[x];
	return ptr;
}
double* util::dvec(int x){
	double *ptr = new double[x];
	return ptr;
}
int** util::imat(int x, int y) {
	int **ptr = new int*[x];
	ptr[0] = new int [x*y];
	for (int i = 1; i < x; i++)
		ptr[i] = ptr[0] + i * y;
	return ptr;
}
double** util::dmat(int x, int y) {
	double **ptr = new double*[x];
	ptr[0] = new double [x*y];
	for (int i = 1; i < x; i++)
		ptr[i] = ptr[0] + i * y;
	return ptr;
}

void util::free_ivec(int *x) {
	delete x;
}

void util::free_dvec(double *x) {
	delete x;
}
void util::free_imat(int **x) {
	delete x[0];
	delete x;
}

void util::free_dmat(double **x) {
	
	delete x[0];delete x;
}
