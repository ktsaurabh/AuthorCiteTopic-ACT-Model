/*---------------------------------------------------
 * file:    topiclib.c
 * purpose: various utility routines for c code
 * version: 1.0
 * author:  newman@uci.edu
 * date:    12/20/05
 *-------------------------------------------------*/

// grep //$ topiclib.c

#include "topiclib.h"

/*------------------------------------------
 * private to util.c
 *------------------------------------------ */

static int icomp(const void *, const void *); /* comparison for isort */
static int dcomp(const void *, const void *); /* comparison for dsort */
static int *icomp_vec; /*  data used for isort */
static double *dcomp_vec; /*  data used for dsort */

/*------------------------------------------
 * allocation routines
 * imat
 * dmat
 * ivec
 * dvec
 *------------------------------------------ */

int **imat(int nr, int nc) //
{
	int ntot = nr * nc;
	int *tmp = (int*) calloc(ntot, sizeof(int));
	int **x = (int**) calloc(nr, sizeof(int*));
	int r;
	assert(tmp);
	assert(x);
	for (r = 0; r < nr; r++)
		x[r] = tmp + nc * r;
	return x;
}

void free_imat(int **x) //
{
	free(x[0]);
	free(x);
}

double **dmat(int nr, int nc) //
{
	int ntot = nr * nc;
	double *tmp = (double*) calloc(ntot, sizeof(double));
	double **x = (double**) calloc(nr, sizeof(double*));
	int r;
	assert(tmp);
	assert(x);
	for (r = 0; r < nr; r++)
		x[r] = tmp + nc * r;
	return x;
}

/*double **d3dmat(int nr, int nc, int h)
 {
 int 2dtot=nr*nc;
 int ntot=2dtot*h;
 double *tmp=(double*) calloc(ntot,sizeof(double));
 double **tmp1 =(double*) calloc(2dtot,sizeof(double*));
 double ***x = (double***) calloc(nr, sizeof (double*));
 int r;
 for (r=0;r<2dtot;r++) tmp1[r]=tmp + h*r;
 for (r=0;r<nr;r++) x[r]=tmp1 + nc*r;
 retunr x;

 }

 int **i3dmat(int nr, int nc, int h)
 {
 int 2dtot=nr*nc;
 int ntot=2dtot*h;
 int *tmp=(int*) calloc(ntot,sizeof(int));
 int **tmp1 =(int*) calloc(2dtot,sizeof(int*));
 int ***x = (int***) calloc(nr, sizeof (int*));
 int r;
 for (r=0;r<2dtot;r++) tmp1[r]=tmp + h*r;
 for (r=0;r<nr;r++) x[r]=tmp1 + nc*r;
 retunr x;

 }
 */

void free_dmat(double **x) //
{
	free(x[0]);
	free(x);
}

float **fmat(int nr, int nc) //
{
	int ntot = nr * nc;
	float *tmp = (float*) calloc(ntot, sizeof(float));
	float **x = (float**) calloc(nr, sizeof(float*));
	int r;
	assert(tmp);
	assert(x);
	for (r = 0; r < nr; r++)
		x[r] = tmp + nc * r;
	return x;
}

void free_fmat(float **x) //
{
	free(x[0]);
	free(x);
}

int *ivec(int n) //
{
	int *x = (int*) calloc(n, sizeof(int));
	assert(x);
	return x;
}

double *dvec(int n) //
{
	double *x = (double*) calloc(n, sizeof(double));
	assert(x);
	return x;
}

float *fvec(int n) //
{
	float *x = (float*) calloc(n, sizeof(float));
	assert(x);
	return x;
}

void free_ivec(int *x) {
	free(x);
}

/*------------------------------------------
 * vector routines
 * imax
 * dmax
 *------------------------------------------ */

int imax(int n, int *x) //
{
	int i, xmax = x[0];
	for (i = 0; i < n; i++)
		xmax = MAX(xmax, x[i]);
	return xmax;
}

double dmax(int n, double *x) //
{
	int i;
	double xmax = x[0];
	for (i = 0; i < n; i++)
		xmax = MAX(xmax, x[i]);
	return xmax;
}

int imin(int n, int *x) //
{
	int i, xmin = x[0];
	for (i = 0; i < n; i++)
		xmin = MIN(xmin, x[i]);
	return xmin;
}

double dmin(int n, double *x) //
{
	int i;
	double xmin = x[0];
	for (i = 0; i < n; i++)
		xmin = MIN(xmin, x[i]);
	return xmin;
}

int isum(int n, int *x) //
{
	int i, xsum = 0;
	for (i = 0; i < n; i++)
		xsum += x[i];
	return xsum;
}

double dsum(int n, double *x) //
{
	int i;
	double xsum = 0;
	for (i = 0; i < n; i++)
		xsum += x[i];
	return xsum;
}

double ddot(int n, double *x, double *y) //
{
	int i;
	double xdot = 0;
	for (i = 0; i < n; i++)
		xdot += x[i] * y[i];
	return xdot;
}

/*------------------------------------------
 * countlines
 * 
 *------------------------------------------ */
int countlines(char *fname) //
{
	int lines = 0;
	char buf[BUFSIZ];
	FILE *fp = fopen(fname, "r");
	assert(fp);
	while (fgets(buf, BUFSIZ, fp))
		lines++;
	fclose(fp);
	lines -= 3; // less 3 header lines
	assert(lines > 0);
	return lines;
}

int countntot(char *fname) //
{
	int i, count, ntot = 0;
	char buf[BUFSIZ];
	FILE *fp = fopen(fname, "r");
	assert(fp);
	for (i = 0; i < 3; i++)
		fgets(buf, BUFSIZ, fp); // skip 3 header lines
	while (fscanf(fp, "%*d%*d%*d%d", &count) != EOF)
		ntot += count;
	fclose(fp);
	assert(ntot > 0);
	return ntot;
}

/*------------------------------------------
 * sort: call qsort library function
 * isort
 * dsort
 *------------------------------------------ */

void isort(int n, int *x, int direction, int *indx) //
{
	int i;
	assert(direction * direction == 1);
	icomp_vec = ivec(n);
	for (i = 0; i < n; i++) {
		icomp_vec[i] = direction * x[i];
		indx[i] = i;
	}
	qsort(indx, n, sizeof(int), icomp);
	free(icomp_vec);
}
static int icomp(const void *pl, const void *p2) {
	int i = *(int *) pl;
	int j = *(int *) p2;
	return (icomp_vec[i] - icomp_vec[j]);
}

int *dsort(int n, double *x) //
{
	int *indx = ivec(n);
	int i;
	dcomp_vec = dvec(n);
	for (i = 0; i < n; i++) {
		dcomp_vec[i] = -x[i];
		indx[i] = i;
	}
	qsort(indx, n, sizeof(int), dcomp);
	free(dcomp_vec);
	return indx;
}
static int dcomp(const void *pl, const void *p2) {
	int i = *(int *) pl;
	int j = *(int *) p2;
	if (dcomp_vec[i] > dcomp_vec[j])
		return 1;
	if (dcomp_vec[i] < dcomp_vec[j])
		return -1;
	if (dcomp_vec[i] == dcomp_vec[j])
		return 0;
	return 0;
}

void read_docwordcountbin(int nnz, int *w, int *d, int *c, char *fname) //
{
	FILE *fp;
	int chk;
	fp = fopen(fname, "r");
	assert(fp);
	fscanf(fp, "%d", &chk);
	fscanf(fp, "%d", &chk);
	fscanf(fp, "%d", &chk);
	chk = fread(w, sizeof(int), nnz, fp);
	assert(chk == nnz);
	chk = fread(d, sizeof(int), nnz, fp);
	assert(chk == nnz);
	chk = fread(c, sizeof(int), nnz, fp);
	assert(chk == nnz);
	fclose(fp);
}

int countnnz(int nr, int nc, int **x) //
{
	int i, j, nnz = 0;
	for (i = 0; i < nr; i++)
		for (j = 0; j < nc; j++)
			if (x[i][j] > 0)
				nnz++;
	return nnz;
}

void write_sparse(int nr, int nc, int **x, char *fname) //
{
	FILE *fp = fopen(fname, "w");
	int i, j;
	assert(fp);
	fprintf(fp, "%d\n", nr);
	fprintf(fp, "%d\n", nc);
	fprintf(fp, "%d\n", countnnz(nr, nc, x));
	for (i = 0; i < nr; i++)
		for (j = 0; j < nc; j++)
			if (x[i][j] > 0)
				fprintf(fp, "%d %d %d\n", i + 1, j + 1, x[i][j]);
	fclose(fp);
}

void write_sparsebin(int nr, int nc, int **x, char *fname) //
{
	int i, j, k, chk;
	int nnz = countnnz(nr, nc, x);
	int *col1 = ivec(nnz);
	int *col2 = ivec(nnz);
	int *col3 = ivec(nnz);
	FILE *fp = fopen(fname, "w");
	assert(fp);
	for (i = 0, k = 0; i < nr; i++)
		for (j = 0; j < nc; j++)
			if (x[i][j] > 0) {
				col1[k] = i;
				col2[k] = j;
				col3[k] = x[i][j];
				k++;
			}
	assert(k == nnz);
	fwrite(&nr, sizeof(int), 1, fp);
	fwrite(&nc, sizeof(int), 1, fp);
	fwrite(&nnz, sizeof(int), 1, fp);
	chk = fwrite(col1, sizeof(int), nnz, fp);
	assert(chk == nnz);
	chk = fwrite(col2, sizeof(int), nnz, fp);
	assert(chk == nnz);
	chk = fwrite(col3, sizeof(int), nnz, fp);
	assert(chk == nnz);
	fclose(fp);
	free(col1);
	free(col2);
	free(col3);
}

int **read_sparse(char *fname, int *nr_, int *nc_) //
{
	FILE *fp = fopen(fname, "r");
	int i, j, c, nr, nc, nnz;
	int **x;
	assert(fp);
	fscanf(fp, "%d", &nr);
	assert(nr > 0);
	fscanf(fp, "%d", &nc);
	assert(nc > 0);
	fscanf(fp, "%d", &nnz);
	assert(nnz > 0);
	x = imat(nr, nc);
	while (fscanf(fp, "%d%d%d", &i, &j, &c) != EOF) {
		i--;
		j--;
		assert(i < nr);
		assert(j < nc);
		assert(c > 0);
		x[i][j] = c;
	}
	fclose(fp);
	*nr_ = nr;
	*nc_ = nc;
	return x;
}

/*void read_dw(char *fname, int *d, int *w, int *D, int *W) //
 {
 int i,wt,dt,ct,count,nnz;
 FILE *fp = fopen(fname ,"r"); assert(fp);
 count = 0;
 fscanf(fp,"%d", D);    assert(*D>0);
 fscanf(fp,"%d", W);    assert(*W>0);
 fscanf(fp,"%d", &nnz); assert(nnz>0);
 while (fscanf(fp, "%d%d%d", &dt, &wt, &ct) != EOF) {
 for (i = count; i < count+ct; i++) {
 w[i] = wt;
 d[i] = dt;
 }
 count += ct;
 }
 fclose(fp);
 }*/

void read_dw(char *fname, int *d, int *w, int *c, int *D, int *W) //
{
	int i, wt, dt, cite, ct, count, nnz;
	FILE *fp = fopen(fname, "r");
	assert(fp);
	count = 0;
	fscanf(fp, "%d", D);
	assert(*D > 0);
	fscanf(fp, "%d", W);
	assert(*W > 0);
	fscanf(fp, "%d", &nnz);
	assert(nnz > 0);

	while (fscanf(fp, "%d%d%d%d", &dt, &cite, &wt, &ct) != EOF) {
		for (i = count; i < count + ct; i++) {
			w[i] = wt;
			d[i] = dt;
			c[i] = cite;
		}
		count += ct;
	}
	fclose(fp);
}

void read_dwc(char *fname, int *d, int *w, int *c, int *D, int *W, int *C) //
{
	int i, wt, dt, cite, ct, count, nnz;
	FILE *fp = fopen(fname, "r");
	assert(fp);
	count = 0;
	fscanf(fp, "%d", D);
	assert(*D > 0);
	fscanf(fp, "%d", W);
	assert(*W > 0);
	fscanf(fp, "%d", C);
	assert(*C > 0);

	while (fscanf(fp, "%d%d%d%d", &dt, &cite, &wt, &ct) != EOF) {
		for (i = count; i < count + ct; i++) {
			w[i] = wt;
			d[i] = dt;
			c[i] = cite;
		}
		count += ct;
	}
	fclose(fp);
}

void fill_dtot(int ntot, int *d, int *dtot) //
{
	int i;
	for (i = 0; i < ntot; i++)
		dtot[d[i]]++;
}

/*void read_dwc(char *fname, int *d, int *w, int *c, int *D, int *W) //
 {
 FILE *fp = fopen(fname,"r");
 int i=0, dd, ww, cc, nnz;
 assert(fp);
 fscanf(fp,"%d", D);    assert(*D>0);
 fscanf(fp,"%d", W);    assert(*W>0);
 fscanf(fp,"%d", &nnz); assert(nnz>0);
 while (fscanf(fp, "%d%d%d", &dd, &ww, &cc) != EOF) {
 d[i] = --dd;
 w[i] = --ww;
 c[i] = cc;
 i++;
 }
 assert(i==nnz);
 fclose(fp);
 }
 */
int read_nnzbin(char *fname) //
{
	int nnz;
	FILE *fp = fopen(fname, "r");
	assert(fp);
	assert(fread(&nnz, sizeof(int), 1, fp)); // nr
	assert(fread(&nnz, sizeof(int), 1, fp)); // nc
	assert(fread(&nnz, sizeof(int), 1, fp)); // nnz
	fclose(fp);
	return nnz;
}

int read_nnz(char *fname) //
{
	int nnz;
	FILE *fp = fopen(fname, "r");
	assert(fp);
	fscanf(fp, "%d", &nnz); // nr
	fscanf(fp, "%d", &nnz); // nc
	fscanf(fp, "%d", &nnz); // nnz
	fclose(fp);
	return nnz;
}

void read_sparsebin(char *fname, int *col1, int *col2, int *col3) //
{
	int nr, nc, nnz, chk;
	FILE *fp = fopen(fname, "r");
	assert(fp);
	assert(fread(&nr, sizeof(int), 1, fp));
	assert(nr > 0);
	assert(fread(&nc, sizeof(int), 1, fp));
	assert(nc > 0);
	assert(fread(&nnz, sizeof(int), 1, fp));
	assert(nnz > 0);
	chk = fread(col1, sizeof(int), nnz, fp);
	assert(chk == nnz);
	chk = fread(col2, sizeof(int), nnz, fp);
	assert(chk == nnz);
	chk = fread(col3, sizeof(int), nnz, fp);
	assert(chk == nnz);
	fclose(fp);
}

void write_ivec(int n, int *x, char *fname) //
{
	FILE *fp = fopen(fname, "w");
	int i;
	assert(fp);
	for (i = 0; i < n; i++)
		fprintf(fp, "%d\n", x[i] + 1);
	fclose(fp);
}

void read_ivec(int n, int *x, char *fname) //
{
	FILE *fp = fopen(fname, "r");
	int i;
	assert(fp);
	for (i = 0; i < n; i++) {
		fscanf(fp, "%d", x + i);
		x[i]--;
	}
	fclose(fp);
}

/*------------------------------------------
 * randperm
 *------------------------------------------ */
/*
 int *randperm(int n) //
 {
 int *order = ivec(n);
 int k, nn, takeanumber, temp;
 for (k = 0; k < n; k++)
 order[k] = k;
 nn = n;
 for (k = 0; k < n; k++) {
 // take a number between 0 and nn-1
 takeanumber = (int) (nn * drand());
 temp = order[nn - 1];
 order[nn - 1] = order[takeanumber];
 order[takeanumber] = temp;
 nn--;
 }
 return order;
 }
 */
/*------------------------------------------
 * randomassignment
 *------------------------------------------ */
void randomassignment(int ntot, int D, int T, int *w, int *d, int *c, int *z,
		int **wp, int **dp, int **cp, int *ztot, int *cztot) //
{
	int i, t, rand_int;
	double rand_double;
	for (i = 0; i < ntot; i++) {
		rand_int = rand() % 1000 + 1;
		rand_double = (double) rand_int / 1000.0;
		t = (int) (T * rand_double);
		z[i] = t;
		wp[w[i]][t]++;
		dp[d[i]][t]++;
		if (c[i] >= 0) {
			cp[c[i]][t]++;
			cztot[t]++;
		}
		ztot[t]++;
	}
	/*for(i=0;i<D;i++){
	 t=(int)(T*drand());
	 cztot[t]++;
	 }*/
}
void randomassignment_atm(int ntot, int D, int T, int *w, int *d, int *c,
		int *z, int *s, int **wp, int **dp, int **cp, int** doc_cite_0,
		int** doc_cite_1, int *ztot, int *cztot, doc_author* map_doc_author,
		int* to_sample, int* to_test) {

	int i, t, coin, rand_int, random, count, count_1;
	for (i = 0; i < ntot; i++)
		if (map_doc_author[d[i]].count == 1) {
			rand_int = rand() % 4;
			if (rand_int == 0) {
				to_test[i] = 1;
				count_1++;
			}
		}
	cout << "tokens to test:" << count_1 << endl;
	double rand_double;
	int doc_seen = -1, citation_seen = -1, doc_seen_prev = -1,
			citation_seen_prev = -1;
	for (i = 0; i < ntot; i++) {
		if (map_doc_author[d[i]].count > 0) {// && to_test[i]!=1){
			to_sample[i] = 1;
			//for(int j=0;j<map_doc_author[d[i]].count;j++){
			//rand_int = rand() % 1000 + 1;
			rand_double = mtrand();//;(double) rand_int / 1000.0;
			t = (int) (T * rand_double);
			z[i] = t;
			wp[w[i]][t]++;
			count = map_doc_author[d[i]].count;
			//cout<<count<<"\t";
			random = rand() % count;
			doc_seen = d[i];
			citation_seen = c[i];
			dp[map_doc_author[d[i]].authors[random]][t]++;
			d[i] = map_doc_author[d[i]].authors[random];
			if (c[i] >= 0 && map_doc_author[c[i]].count > 0) {
				if (doc_seen_prev != doc_seen || citation_seen_prev
						!= citation_seen) {
					//rand_int = rand() % 1000 + 1;
					rand_double = mtrand();//(double) rand_int / 1000.0;
					coin = (int) (2 * rand_double);
					s[i] = coin;
					count = map_doc_author[c[i]].count;
					random = rand() % count;
					cp[map_doc_author[c[i]].authors[random]][t]++;
					c[i] = map_doc_author[c[i]].authors[random];
					cztot[t]++;
					doc_seen_prev = doc_seen;
					citation_seen_prev = citation_seen;
				}
			} else if (c[i] >= 0 && map_doc_author[c[i]].count == 0) {
				c[i] = -1;
				s[i] = -1;
			} else
				s[i] = -1;

			ztot[t]++;
		}
		//}
	}
	/*for(i=0;i<D;i++){
	 t=(int)(T*drand());
	 cztot[t]++;
	 }*/
}

void randomassignment_alt(int ntot, int D, int T, int *w, int *d, int *c,
		int *z, int *s, int **wp, int **dp, int **cp, int** doc_cite_0,
		int** doc_cite_1, int *ztot, int *cztot, doc_author* map_doc_author,
		int* to_sample, int* to_test) {

	int i = 0, t = 0, coin = 0, rand_int = 0, random = 0, count = 0, count_1 =
			0;
	cout << ntot << endl;
	for (i = 0; i < ntot; i++)
		if (map_doc_author[d[i]].count == 1) {
			coin++;
			rand_int = rand() % 4;
			if (rand_int == 0) {
				//cout<<"here"<<endl;
				to_test[i] = 1;
				count_1++;
			}
		}
	cout << "total tokens:" << coin << endl;
	cout << "tokens to test:" << count_1 << endl;
	double rand_double;
	int doc_seen = -1, citation_seen = -1, doc_seen_prev = -1,
			citation_seen_prev = -1;
	for (i = 0; i < ntot; i++) {
		if (map_doc_author[d[i]].count > 0) {// && to_test[i]!=1){
			to_sample[i] = 1;
			//for(int j=0;j<map_doc_author[d[i]].count;j++){
			//rand_int = rand() % 1000 + 1;
			rand_double = mtrand();//(double) rand_int / 1000.0;
			t = (int) (T * rand_double);
			z[i] = t;
			wp[w[i]][t]++;
			count = map_doc_author[d[i]].count;
			doc_seen = d[i];
			citation_seen = c[i];
			random = rand() % count;
			dp[map_doc_author[d[i]].authors[random]][t]++;
			d[i] = map_doc_author[d[i]].authors[random];
			if (c[i] >= 0 && map_doc_author[c[i]].count > 0) {
				if (doc_seen_prev != doc_seen || citation_seen_prev
						!= citation_seen) {
					//rand_int = rand() % 1000 + 1;
					rand_double = mtrand();//(double) rand_int / 1000.0;
					coin = (int) (2 * rand_double);
					s[i] = coin;
					count = map_doc_author[c[i]].count;
					random = rand() % count;
					cp[map_doc_author[c[i]].authors[random]][t]++;
					c[i] = map_doc_author[c[i]].authors[random];
					cztot[t]++;
					doc_seen_prev = doc_seen;
					citation_seen_prev = citation_seen;

				}
			} else if (c[i] >= 0 && map_doc_author[c[i]].count == 0) {
				c[i] = -1;
				s[i] = -1;
			} else
				s[i] = -1;

			ztot[t]++;
		}
	}
}
void randomassignment_coin_sact(int ntot, int D, int T, int *w, int *d, int *c,
		int *z, int *s, int **wp, int **dp, int **cp, int** doc_cite_0,
		int** doc_cite_1, int *ztot, int *cztot, doc_author* map_doc_author,
		int* to_sample, int* to_test) //
{

	int i, t, coin, rand_int, random, count, count_1;

	double rand_double;
	for (i = 0; i < ntot; i++) {
		if (map_doc_author[d[i]].count > 0) {// && to_test[i]!=1){
			to_sample[i] = 1;
			//for(int j=0;j<map_doc_author[d[i]].count;j++){
			rand_int = rand() % 1000 + 1;
			rand_double = (double) rand_int / 1000.0;
			t = (int) (T * rand_double);
			z[i] = t;
			wp[w[i]][t]++;
			count = map_doc_author[d[i]].count;
			//cout<<count<<"\t";
			random = rand() % count;
			dp[map_doc_author[d[i]].authors[random]][t]++;
			d[i] = map_doc_author[d[i]].authors[random];
			if (c[i] >= 0 && map_doc_author[c[i]].count > 0) {
				rand_int = rand() % 1000 + 1;
				rand_double = (double) rand_int / 1000.0;
				coin = (int) (2 * rand_double);
				s[i] = coin;
				if (coin)
					doc_cite_1[d[i]][c[i]]++;
				else
					doc_cite_0[d[i]][c[i]]++;
				//for(int k=0;k<map_doc_author[c[i]].count;k++){
				count = map_doc_author[c[i]].count;
				//cout<<count<<"\t";
				random = rand() % count;
				cp[map_doc_author[c[i]].authors[random]][t]++;
				c[i] = map_doc_author[c[i]].authors[random];
				cztot[t]++;
				//}
			} else if (c[i] >= 0 && map_doc_author[c[i]].count == 0) {
				c[i] = -1;
				s[i] = -1;
			} else
				s[i] = -1;
			ztot[t]++;
		}
		//}
	}
	/*for(i=0;i<D;i++){
	 t=(int)(T*drand());
	 cztot[t]++;
	 }*/
	int count_s1 = 0, count_s2 = 0;
	for (i = 0; i < ntot; i++) {
		if (s[i] == 1) {
			count_s1++;
		} else if (s[i] == 0) {
			count_s2++;
		}
	}
	cout << "percent on=" << (double) count_s1 / (double) (count_s1 + count_s2)
			<< endl;

}
void randomassignment_coin_tu_etal(int ntot, int D, int T, int *w, int *d,
		int *c, int *z, int *s, int **wp, int **dp, int **cp, int** doc_cite_0,
		int** doc_cite_1, int *ztot, int *cztot, doc_author* map_doc_author,
		doc_author* map_doc_cited_author, int* to_sample, int* to_test) {
	int i, t, coin, rand_int, random, count, count_1;

	for (i = 0; i < ntot; i++)
		if (map_doc_author[d[i]].count == 1) {
			rand_int = rand() % 4;
			if (rand_int == 0) {
				to_test[i] = 1;
				count_1++;
			}
		}

	double rand_double;
	for (i = 0; i < ntot; i++) {
		if (map_doc_author[d[i]].count > 0) {// && to_test[i]!=1){
			to_sample[i] = 1;
			//for(int j=0;j<map_doc_author[d[i]].count;j++){
			rand_int = rand() % 1000 + 1;
			rand_double = (double) rand_int / 1000.0;
			t = (int) (T * rand_double);
			z[i] = t;
			wp[w[i]][t]++;
			count = map_doc_author[d[i]].count;
			//cout<<count<<"\t";
			random = rand() % count;
			dp[map_doc_author[d[i]].authors[random]][t]++;
			int doc = d[i];
			d[i] = map_doc_author[d[i]].authors[random];
			if (c[i] >= 0 && map_doc_author[c[i]].count > 0) {
				for (int iter_cited = 0; iter_cited
						< map_doc_cited_author[doc].count; iter_cited++) {
					//rand_int = rand() % 1000 + 1;
					rand_double = (double) rand_int / 1000.0;
					coin = (int) (2 * rand_double);
					s[i] = coin;
					//if(coin)
					//doc_cite_1[d[i]][c[i]]++;
					//else
					//doc_cite_0[d[i]][c[i]]++;
					//for(int k=0;k<map_doc_author[c[i]].count;k++){
					//count = map_doc_author[c[i]].count;
					//cout<<count<<"\t";
					//random = rand() % count;
					cp[map_doc_author[c[i]].authors[iter_cited]][t]++;
					//c[i] = map_doc_author[c[i]].authors[iter_cited];
					cztot[t]++;
					//}
				}
			} else if (c[i] >= 0 && map_doc_author[c[i]].count == 0) {
				c[i] = -1;
				s[i] = -1;
			} else
				s[i] = -1;
			ztot[t]++;
		}
		//}
	}
	/*for(i=0;i<D;i++){
	 t=(int)(T*drand());
	 cztot[t]++;
	 }*/

}
void randomassignment_act(int ntot, int D, int T, int *w, int *d, int *c,
		int *z, int *s, int **wp, int **dp, int **cp, int** doc_cite_0,
		int** doc_cite_1, int *ztot, int *cztot, doc_author* map_doc_author,
		int* to_sample, int* to_test, int percent) //
{

	int i = 0, t = 0, coin = 0, rand_int = 0, random = 0, count = 0, count_1 = 0;
	for (i = 0; i < ntot; i++)
		if (map_doc_author[d[i]].count == 1) {
			coin++;
			rand_int = rand() % percent;
			if (rand_int == 0) {
				to_test[i] = 1;
				count_1++;
			} else
				to_test[i] = 0;
		}
	cout << "total tokens:" << coin << endl;
	cout << "tokens to test:" << count_1 << endl;
	double rand_double;
	for (i = 0; i < ntot; i++) {
		if (map_doc_author[d[i]].count > 0) {// && to_test[i]!=1){
			to_sample[i] = 1;
			//for(int j=0;j<map_doc_author[d[i]].count;j++){
			//rand_int = rand() % 1000 + 1;
			rand_double = mtrand();//(double) rand_int / 1000.0;
			t = (int) (T * rand_double);
			z[i] = t;
			wp[w[i]][t]++;
			count = map_doc_author[d[i]].count;
			//cout<<count<<"\t";
			random = rand() % count;
			dp[map_doc_author[d[i]].authors[random]][t]++;
			d[i] = map_doc_author[d[i]].authors[random];
			if (c[i] >= 0 && map_doc_author[c[i]].count > 0) {
				//rand_int = rand() % 1000 + 1;
				rand_double = mtrand();//(double) rand_int / 1000.0;
				coin = (int) (2 * rand_double);
				s[i] = coin;
				//if(coin)
				//doc_cite_1[d[i]][c[i]]++;
				//else
				//doc_cite_0[d[i]][c[i]]++;
				//for(int k=0;k<map_doc_author[c[i]].count;k++){
				count = map_doc_author[c[i]].count;
				//cout<<count<<"\t";
				random = rand() % count;
				cp[map_doc_author[c[i]].authors[random]][t]++;
				c[i] = map_doc_author[c[i]].authors[random];
				cztot[t]++;
				//}
			} else if (c[i] >= 0 && map_doc_author[c[i]].count == 0) {
				c[i] = -1;
				s[i] = -1;
			} else
				s[i] = -1;
			ztot[t]++;
		}
		//}
	}
	/*for(i=0;i<D;i++){
	 t=(int)(T*drand());
	 cztot[t]++;
	 }*/
}

/*void randomassignment_rank(int ntot, int T, int *w, int *d, int *drank, int *z,
 int **wp, int **dp, int *ztot) //
 {
 int i, t;
 for (i = 0; i < ntot; i++) {
 t = (int) (T * drand());
 z[i] = t;
 wp[w[i]][t] += drank[d[i]];
 dp[d[i]][t] += drank[d[i]];
 ztot[t] += drank[d[i]];
 }
 }
 */
void assignment(int ntot, int *w, int *d, int *z, int **wp, int **dp, int *ztot) //
{
	int i, t;
	for (i = 0; i < ntot; i++) {
		t = z[i];
		wp[w[i]][t]++;
		dp[d[i]][t]++;
		ztot[t]++;
	}
}
/*
 void randomassignment_2layer(int ntot, int T, int S, int *w, int *d, int *z,
 int *y, int **wp, int **zy, int **dp, int *ztot, int *ytot) //
 {
 int i, t, s;
 for (i = 0; i < ntot; i++) {
 t = (int) (T * drand());
 s = (int) (S * drand());
 z[i] = t;
 y[i] = s;
 wp[w[i]][t]++;
 zy[t][s]++;
 dp[d[i]][s]++;
 ztot[t]++;
 ytot[s]++;
 }
 }

 void randomassignment2(int ntot, int T, int *d, int *z, int **dp) //
 {
 int i, t;
 for (i = 0; i < ntot; i++) {
 t = (int) (T * drand());
 z[i] = t;
 dp[d[i]][t]++;
 }
 }
 */
/*------------------------------------------
 * sample_chain
 *------------------------------------------ */
/*void sample_chain(int ntot, int W, int T, double alpha, double beta, int *w,
 int *d, int *z, int **wp, int **dp, int *ztot, int *order) //
 {
 int ii, i, t;
 double totprob, maxprob, currprob;
 double *probs = dvec(T);
 double wbeta = W * beta;

 for (ii = 0; ii < ntot; ii++) {

 i = order[ii];

 t = z[i]; // take the current topic assignment to word token i
 ztot[t]--; // and substract that from the counts
 wp[w[i]][t]--;
 dp[d[i]][t]--;
 totprob = 0;

 for (t = 0; t < T; t++) {
 probs[t] = (wp[w[i]][t] + beta) / (ztot[t] + wbeta) * (dp[d[i]][t]
 + alpha);
 totprob += probs[t];
 }

 maxprob = drand() * totprob;
 currprob = probs[0];
 t = 0;

 // sample a topic t from the distribution
 while (maxprob > currprob) {
 t++;
 currprob += probs[t];
 }

 z[i] = t; // assign current word token i to topic t
 wp[w[i]][t]++; // and update counts
 dp[d[i]][t]++;
 ztot[t]++;
 }

 free(probs);
 }
 */
/*------------------------------------------
 * sample_chain_rank
 *------------------------------------------ */
/*void sample_chain_rank(int ntot, int W, int T, double alpha, double beta,
 int *w, int *d, int *drank, int *z, int **wp, int **dp, int *ztot,
 int *order) //
 {
 int ii, i, t;
 double totprob, maxprob, currprob;
 double *probs = dvec(T);
 double wbeta = W * beta;

 for (ii = 0; ii < ntot; ii++) {

 i = order[ii];

 t = z[i]; // take the current topic assignment to word token i
 ztot[t] -= drank[d[i]];
 wp[w[i]][t] -= drank[d[i]];
 dp[d[i]][t] -= drank[d[i]];
 totprob = 0;

 for (t = 0; t < T; t++) {
 probs[t] = (wp[w[i]][t] + beta) / (ztot[t] + wbeta) * (dp[d[i]][t]
 + alpha);
 totprob += probs[t];
 }

 maxprob = drand() * totprob;
 currprob = probs[0];
 t = 0;

 // sample a topic t from the distribution
 while (maxprob > currprob) {
 t++;
 currprob += probs[t];
 }

 z[i] = t; // assign current word token i to topic t
 wp[w[i]][t] += drank[d[i]];
 dp[d[i]][t] += drank[d[i]];
 ztot[t] += drank[d[i]];
 }

 free(probs);
 }

 void sample_chain0(int ntot, int W, int T, double alpha, double beta, int *w,
 int *d, int *z, int **wp, int **dp, int *ztot) //
 {
 int i, t;
 double totprob, maxprob, currprob;
 double *probs = dvec(T);
 double wbeta = W * beta;

 for (i = 0; i < ntot; i++) {

 t = z[i]; // take the current topic assignment to word token i
 ztot[t]--; // and substract that from the counts
 wp[w[i]][t]--;
 dp[d[i]][t]--;

 for (t = 0, totprob = 0.0; t < T; t++) {
 probs[t] = (dp[d[i]][t] + alpha) * (wp[w[i]][t] + beta) / (ztot[t]
 + wbeta);
 totprob += probs[t];
 }

 #if 1
 assert(totprob > 0.0);
 #endif

 maxprob = drand() * totprob;
 currprob = probs[0];
 t = 0;

 // sample a topic t from the distribution
 while (maxprob > currprob) {
 t++;
 currprob += probs[t];
 }

 z[i] = t; // assign current word token i to topic t
 wp[w[i]][t]++; // and update counts
 dp[d[i]][t]++;
 ztot[t]++;
 }

 free(probs);
 }

 void sample_chain_2layer(int ntot, int W, int T, int S, double alpha,
 double beta, double gamma, int *w, int *d, int *z, int *y, int **wp,
 int **zy, int **dp, int *ztot, int *ytot) //
 {
 int i, t, s;
 double totprob, maxprob, currprob, term1, term2, term3;
 double wbeta = W * beta;
 double tgamma = T * gamma;
 double **probs = dmat(T, S);

 for (i = 0; i < ntot; i++) {

 t = z[i];
 s = y[i];
 ztot[t]--;
 ytot[s]--;
 wp[w[i]][t]--;
 zy[t][s]--;
 dp[d[i]][s]--;

 totprob = 0;
 for (t = 0; t < T; t++) {
 for (s = 0; s < S; s++) {
 term1 = (wp[w[i]][t] + beta) / (ztot[t] + wbeta);
 term2 = (zy[t][s] + gamma) / (ytot[s] + tgamma);
 term3 = (dp[d[i]][s] + alpha);
 probs[t][s] = term1 * term2 * term3;
 totprob += probs[t][s];
 }
 }

 maxprob = drand() * totprob;
 currprob = probs[0][0];
 t = 0;
 s = 0;
 while (maxprob > currprob) {
 t++;
 if (t >= T) {
 s++;
 t = 0;
 }
 currprob += probs[t][s];
 }

 z[i] = t;
 y[i] = s;
 ztot[t]++;
 ytot[s]++;
 wp[w[i]][t]++;
 zy[t][s]++;
 dp[d[s]][t]++;
 }

 free_dmat(probs);

 }

 void resample_chain(int ntot, int W, int T, double alpha, double beta, int *w,
 int *d, int *z, int **wp, int **dp, int *ztot) //
 {
 int i, t;
 double totprob, maxprob, currprob;
 double *probs = dvec(T);
 double wbeta = W * beta;

 for (i = 0; i < ntot; i++) {

 t = z[i];
 dp[d[i]][t]--;
 totprob = 0;

 for (t = 0; t < T; t++) {
 probs[t] = (wp[w[i]][t] + beta) / (ztot[t] + wbeta) * (dp[d[i]][t]
 + alpha);
 totprob += probs[t];
 }

 maxprob = drand() * totprob;
 currprob = probs[0];
 t = 0;

 // sample a topic t from the distribution
 while (maxprob > currprob) {
 t++;
 currprob += probs[t];
 }

 z[i] = t;
 dp[d[i]][t]++;
 }

 free(probs);
 }

 void oversample_dp(int ntot, int W, int T, double alpha, double beta, int *w,
 int *d, int *z, int **wp, int **dp, int *ztot) //
 {
 int i, t, k, ntimes = 4;
 double totprob, maxprob, currprob;
 double *probs = dvec(T);
 double wbeta = W * beta;

 for (i = 0; i < ntot; i++) {

 totprob = 0;
 for (t = 0; t < T; t++) {
 probs[t] = (wp[w[i]][t] + beta) / (ztot[t] + wbeta) * (dp[d[i]][t]
 + alpha);
 totprob += probs[t];
 }

 for (k = 0; k < ntimes; k++) {
 maxprob = drand() * totprob;
 currprob = probs[0];
 t = 0;
 while (maxprob > currprob) {
 t++;
 currprob += probs[t];
 }
 dp[d[i]][t]++;
 }
 }

 free(probs);
 }
 */
void loglike(int ntot, int W, int D, int T, double alpha, double beta, int *w,
		int *d, int **wp, int **dp, int *ztot, int *dtot) //
{
	int i, j, t;
	double llike;
	static int init = 0;
	static double **prob_w_given_t;
	static double **prob_t_given_d;
	static double *dtot_;
	double ztot_;

	if (init == 0) {
		init = 1;
		prob_w_given_t = dmat(W, T);
		prob_t_given_d = dmat(D, T);
		dtot_ = dvec(D);
		for (j = 0; j < D; j++)
			dtot_[j] = dtot[j] + T * alpha;
	}

	for (t = 0; t < T; t++) {
		ztot_ = ztot[t] + W * beta;
		for (i = 0; i < W; i++)
			prob_w_given_t[i][t] = (wp[i][t] + beta) / ztot_;
		for (j = 0; j < D; j++)
			prob_t_given_d[j][t] = (dp[j][t] + alpha) / dtot_[j];
	}

	llike = 0;
	for (i = 0; i < ntot; i++)
		llike += log(ddot(T, prob_w_given_t[w[i]], prob_t_given_d[d[i]]));

	printf(">>> llike = %.6e    ", llike);
	printf("pplex = %.4f\n", exp(-llike / ntot));
}

void chksum(int n, int T, int **x, int *sumx) //
{
	int i, t, sum;
	for (t = 0; t < T; t++) {
		sum = 0;
		for (i = 0; i < n; i++)
			sum += x[i][t];
		assert(sum == sumx[t]);
	}
}

void getztot(int n, int T, int **x, int *ztot) //
{
	int i, t, sum;
	for (t = 0; t < T; t++) {
		sum = 0;
		for (i = 0; i < n; i++)
			sum += x[i][t];
		ztot[t] = sum;
	}
}

double etime() //
{
	static double last_clock = 0;
	static double now_time = 0;
	last_clock = now_time;
	now_time = (double) clock();
	return (double) (now_time - last_clock) / CLOCKS_PER_SEC;
}
void quickSort(double *numbers, int *index, int array_size) {
	q_sort(numbers, index, 0, array_size - 1);

}

void q_sort(double *numbers, int *index, int left, int right) {
	long pivot, pivot_index, l_hold, r_hold;

	l_hold = left;
	r_hold = right;
	pivot = numbers[left];
	pivot_index = index[left];
	while (left < right) {
		while ((numbers[right] <= pivot) && (left < right))
			right--;
		if (left != right) {
			numbers[left] = numbers[right];
			index[left] = index[right];
			left++;
		}
		while ((numbers[left] >= pivot) && (left < right))
			left++;
		if (left != right) {
			numbers[right] = numbers[left];
			index[right] = index[left];
			right--;
		}
	}
	numbers[left] = pivot;
	index[left] = pivot_index;
	pivot = left;
	left = l_hold;
	right = r_hold;
	if (left < pivot)
		q_sort(numbers, index, left, pivot - 1);
	if (right > pivot)
		q_sort(numbers, index, pivot + 1, right);
}
