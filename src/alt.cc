/*---------------------------------------------------
 * file:    topicfast.c
 * purpose: run standard topic model with shortcut sampling
 * usage:   topic2 -help
 * inputs:  T (number of topics), iter (number of iterations), docword.txt
 * outputs: wp.txt, dp.txt, z.txt
 * version: 1.0
 * author:  saurabh@psu.edu (modified version provided by newman@uci.edu)
 * date:    2008
 *-------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <math.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include "topiclib.h"
#include "util.h"
#include <map>
using namespace std;
/*------------------------------------------
 * global variables
 *------------------------------------------ */
static int T; // number of topics
static int W; // number of unique words
static int D; // number of docs
static int N; // number of words in corpus
static int C; //number of cited documents.
static int A;
static int CA;

doc_author* read_authors(char* file, doc_author* da, int *A) {// file format: number of authors\n doc \t count of authors \t author id seperated with tabs
	ifstream f;
	f.open(file);
	int doc, count, author;
	f >> *A;
	while (!f.eof()) {
		f >> doc >> count;
		da[doc].count = count;
		da[doc].authors = new int[count];
		for (int i = 0; i < count; i++) {
			f >> author;
			da[doc].authors[i] = author;
		}
	}
	return da;
}

struct rank_word {
	int word;
	double val;
};
struct rank_cite {
	int citation;
	double val;
};

bool SortGreater_word(const rank_word& first, const rank_word& second) {
	return first.val > second.val;
}
bool SortGreater_cite(const rank_cite& first, const rank_cite& second) {
	return first.val > second.val;
}

void print_top(char *file_to_write, char *map, double **topic_doc, int row,
		int col) {
	int count = 0;
	char s[100], d = ':', *name;
	char **names = new char*[col];
	int i = 0, id;
	for (i = 0; i < col; i++)
		names[i] = new char[20];
	ifstream in;
	in.open(map);
	while (!in.eof()) {
		in.getline(s, 100);
		char *p = strtok(s, ":");
		if (strcmp(s, "\0") == 0)
			break;
		i = 0;
		while (p != NULL) {
			if (i++ % 2 == 0)
				name = p;
			else
				id = atoi(p);
			p = strtok(NULL, ":");
		}
		strcpy(names[id], name);

	}
	ofstream fo;
	fo.open(file_to_write);
	vector<rank_cite> citation;
	vector<rank_cite>::iterator it;
	for (int i = 0; i < row; i++) {
		citation.clear();
		for (int j = 0; j < col; j++) {
			rank_cite r;
			r.citation = j;
			r.val = topic_doc[i][j];
			citation.push_back(r);
		}

		sort(citation.begin(), citation.end(), SortGreater_cite);
		it = citation.begin();
		count = 0;
		fo << "Topic-" << i << ":";
		while (count++ < 25) {
			fo << names[it->citation] << "=" << it++->val << ":";
		}
		fo << endl;
	}

}

void calculate_rank(int* test_docs, int* cited_doc, int** doc_citation,
		int* citations_2_pred, double **topic_cite, double **topic_doc, int W,
		int A, int T) {
	//for each word, make a vector and sort that vector.
	int max_rank = -1, mean_max_rank = 0, min_rank = 2 * D, mean_min_rank = 0,
			total_ret = 0, rel_ret = 0, total_ret1 = 0, rel_ret1 = 0,

			total_ret2 = 0, rel_ret2 = 0, total_ret5 = 0, rel_ret5 = 0,
			rel_ret3 = 0, rel_ret4 = 0, not_found = 0;
	double prec = 0.0, prec1 = 0.0, prec2 = 0.0, prec3 = 0.0, prec4 = 0.0,
			prec5 = 0;
	double recl = 0.0, recl1 = 0.0, recl2 = 0.0, recl3 = 0.0, recl4 = 0.0,
			recl5 = 0;
	vector<rank_cite> citation;
	vector<rank_cite>::iterator it;
	int count = 0;

	for (int j = 0; j < A; j++) {
		if (test_docs[j] == 1) {
			count++;
			citation.clear();
			for (int i = 0; i < A; i++) {
				if (cited_doc[i] == 1) {
					rank_cite r;
					r.citation = i;
					r.val = 0.0;
					for (int k = 0; k < T; k++) {
						r.val += topic_cite[k][i] * topic_doc[k][j];
						/*if (topic_cite[k][i] > 0 && topic_doc[k][j] > 0)
						 r.val += ((topic_cite[k][i] * log(topic_cite[k][i]
						 / topic_doc[k][j])) + (topic_doc[k][j]
						 * log(topic_doc[k][j] / topic_cite[k][i])));*/
					}
					citation.push_back(r);
				}
			}
			sort(citation.begin(), citation.end(), SortGreater_cite);
			//sort(citation.begin(), citation.end(), SortAscend_cite);
			/*for (it = citation.begin(); it != citation.end(); ++it){
			 cout<<it->val<<endl;
			 }
			 exit(0);*/
			rel_ret = 0;
			int k = 0;
			rel_ret1 = 0;
			rel_ret2 = 0;
			rel_ret3 = 0;
			rel_ret4 = 0;
			rel_ret5 = 0;
			max_rank = 0;
			min_rank = 2 * A;
			for (it = citation.begin(), k = 0; it != citation.end(); ++it, ++k) {
				if (doc_citation[j][it->citation] == 1) {
					if (k <= 1) {
						rel_ret1++;
						//cout<<"k:"<<k<<":"<<it->val<<endl;
					}
					if (k <= 2) {
						rel_ret2++;
						//cout<<"k:"<<k<<":"<<it->val<<endl;
					}
					if (k <= 3) {
						rel_ret3++;
						//cout<<"k:"<<k<<":"<<it->val<<endl;
					}
					if (k <= 4) {
						rel_ret4++;
						//cout<<"k:"<<k<<":"<<it->val<<endl;
					}
					if (k <= 5)
						rel_ret5++;

					if (k <= 10)
						rel_ret++;
					if (min_rank == 2 * A)
						min_rank = k;
					max_rank = k;

					//cout << "Document=" << j << ":Citation=" << it->citation<< ":Rank=" << k << ":val=" << it->val <<endl;
				}
			}
			int max_citations = 0;
			if (min_rank != 2 * A) {
				max_citations = citations_2_pred[j] > 10 ? 10
						: citations_2_pred[j];
				prec += (double) rel_ret / (double) max_citations;
				recl += (double) rel_ret / (double) citations_2_pred[j];
				prec1 += (double) rel_ret1;
				recl1 += (double) rel_ret1 / (double) citations_2_pred[j];
				max_citations = citations_2_pred[j] > 2 ? 2
						: citations_2_pred[j];
				prec2 += (double) rel_ret2 / (double) max_citations;
				recl2 += (double) rel_ret2 / (double) citations_2_pred[j];
				max_citations = citations_2_pred[j] > 3 ? 3
						: citations_2_pred[j];
				prec3 += (double) rel_ret3 / (double) max_citations;
				recl3 += (double) rel_ret3 / (double) citations_2_pred[j];
				max_citations = citations_2_pred[j] > 4 ? 4
						: citations_2_pred[j];
				prec4 += (double) rel_ret4 / (double) max_citations;
				recl4 += (double) rel_ret4 / (double) citations_2_pred[j];
				max_citations = citations_2_pred[j] > 5 ? 5
						: citations_2_pred[j];
				prec5 += (double) rel_ret5 / (double) max_citations;
				recl5 += (double) rel_ret5 / (double) citations_2_pred[j];
				mean_max_rank += max_rank;
				mean_min_rank += min_rank;
			}
			if (min_rank > 400)
				not_found++;
		}

	}
	cout << "Result:Mean Max Rank=" << mean_max_rank / (count - not_found)
			<< ":Mean Min Rank=" << mean_min_rank / (count - not_found)
			<< ":P@10=" << prec / (count - not_found) << ":P@1=" << prec1
			/ (count - not_found) << ":P@2=" << prec2 / (count - not_found)
			<< ":P@3=" << prec3 / (count - not_found) << ":P@4=" << prec4
			/ (count - not_found) << ":P@5=" << prec5 / (count - not_found)
			<< ":not found=" << not_found << ":count=" << count << endl;
	cout << "----------Precision-Recall-Start-------------" << endl;
	cout << "1:" << prec1 / (count - not_found) << ":" << recl1 / (count
			- not_found) << endl;
	cout << "2:" << prec2 / (count - not_found) << ":" << recl2 / (count
			- not_found) << endl;
	cout << "3:" << prec3 / (count - not_found) << ":" << recl3 / (count
			- not_found) << endl;
	cout << "4:" << prec4 / (count - not_found) << ":" << recl4 / (count
			- not_found) << endl;
	cout << "5:" << prec5 / (count - not_found) << ":" << recl5 / (count
			- not_found) << endl;
	cout << "10:" << prec / (count - not_found) << ":" << recl / (count
			- not_found) << endl;
	cout << "----------Precision-Recall-End-------------" << endl;

}

/*void calculate_rank(int* test_docs, int* cited_doc, int** doc_citation,
 double **topic_cite, double **topic_doc, int W, int A, int T) {
 //for each word, make a vector and sort that vector.
 int max_rank = -1, mean_max_rank = 0, min_rank = 2 * D, mean_min_rank = 0,
 total_ret = 0, rel_ret = 0, total_ret1 = 0, rel_ret1 = 0,
 total_ret2 = 0, rel_ret2 = 0, total_ret5 = 0, rel_ret5 = 0;
 double prec = 0.0, prec1 = 0.0, prec2 = 0.0, prec5 = 0;
 vector<rank_cite> citation;
 vector<rank_cite>::iterator it;
 int count = 0;

 for (int j = 0; j < A; j++) {
 if (test_docs[j] == 1) {
 count++;
 citation.clear();
 for (int i = 0; i < A; i++) {
 if (cited_doc[i] > 0) {

 rank_cite r;
 r.citation = i;
 r.val = 0.0;
 for (int k = 0; k < T; k++) {

 r.val += topic_cite[k][i] * topic_doc[k][j];

 }
 citation.push_back(r);
 }
 }
 sort(citation.begin(), citation.end(), SortGreater_cite);

 rel_ret = 0;
 int k = 0;
 rel_ret1 = 0;
 rel_ret2 = 0;
 rel_ret5 = 0;
 max_rank = 0;
 min_rank = 2 * A;
 for (it = citation.begin(), k = 0; it != citation.end(); ++it, ++k) {
 if (doc_citation[j][it->citation] == 1) {
 if (k <= 1)
 rel_ret1++;
 if (k <= 2)
 rel_ret2++;
 if (k <= 5)
 rel_ret5++;

 if (k <= 10)
 rel_ret++;
 if (min_rank == 2 * A)
 min_rank = k;
 max_rank = k;

 //cout << "Document=" << j << ":Citation=" << it->citation<< ":Rank=" << k << ":val=" << it->val <<endl;
 }
 }
 if (min_rank != 2 * A) {
 prec += (double) rel_ret / 10.0;
 prec1 += (double) rel_ret1;
 prec2 += (double) rel_ret2 / 2.0;
 prec5 += (double) rel_ret5 / 5.0;
 mean_max_rank += max_rank;
 mean_min_rank += min_rank;
 }
 }

 }
 cout << "Result:Mean Max Rank=" << mean_max_rank / count
 << ":Mean Min Rank=" << mean_min_rank / count << ":P@10=" << prec
 / count << ":P@1=" << prec1 / count << ":P@2=" << prec2 / count
 << ":P@5=" << prec5 / count << endl;

 }
 */
void print_top_words_given_zc(double **topic_cite, double **occur_topic, int W,
		int D, int T) {
	//for each word, make a vector and sort that vextor.
	vector<rank_word> word;
	vector<rank_word>::iterator it;

	for (int j = 0; j < D; j++)
		for (int k = 0; k < T; k++) {

			word.clear();
			rank_word r;
			for (int i = 0; i < W; i++) {
				r.word = i;
				r.val = topic_cite[k][j] * occur_topic[k][i];
				word.push_back(r);
			}
			sort(word.begin(), word.end(), SortGreater_word);
			it = word.begin();
			for (int i = 0; i < 10; i++) {
				cout << "Citation=" << j << ":Topic=" << k << ":Word="
						<< it->word << ":value=" << it->val << endl;
				++it;
			}
		}

}

void print_top_citation_given_d(double **topic_cite, double **topic_doc, int W,
		int D, int T) {
	//for each word, make a vector and sort that vextor.
	vector<rank_cite> cite;
	vector<rank_cite>::iterator it;

	for (int j = 0; j < D; j++)
		for (int k = 0; k < T; k++) {

			cite.clear();
			rank_cite r;
			for (int i = 0; i < W; i++) {
				r.citation = i;
				r.val = topic_cite[k][j] * topic_doc[k][i];
				cite.push_back(r);
			}
			sort(cite.begin(), cite.end(), SortGreater_cite);
			it = cite.begin();
			for (int i = 0; i < 10; i++) {
				cout << "Citation=" << j << ":Topic=" << k << ":Citation="
						<< it->citation << ":value=" << it->val << endl;
				++it;
			}
		}

}
void generate_test(int *array, int N) {
	int count = -1;
	while (++count < N) {
		array[count] = 0;
	}
	count = 0;
	while (count < 0.2 * N) {
		int r = rand() % N;
		if (array[r] == 0) {
			array[r] = 1;
			count++;
		}
	}
	cout << "count=" << count << endl;
}

int main(int argc, char* argv[]) {
	int i = 0, iter, ITER, seed = 1, t, tt, tt_, tsp, ttt = -1, ttt_chk, t_,
			Tchk, eval_mode, eval_iter, total_iters;
	int *w, *d, *c, *z, **wp, **dp, **cp, *cztot, *ztot, *dl, *indx, *nsamp,
			*nsamplda, *d_last_z, *w_last_z, *ztot_last, ztotmin = -1, d_old,
			w_old, *indx_z, *revindx_z, **indx_dp, *cited_doc, **revindx_dp,
			rand_int;
	double alpha, final_eval, beta, Wbeta, ro, Cro, coin, *probs, *Ad, *Bw, Z,
			Zp = -1, Zp_old, currprob, *a, *b, aa = -1, bb = -1, U, sumprobs,
			sumprobs_old, fac = -1, **occur_topic, **topic_cite, *word_probs,
			**topic_doc, Kalpha, **eval, *log_eval, alpha_lambda_theta,
			alpha_lambda_phi, rand_double;
	int *ldacounts, *fastcounts, *s;
	int word, doc, citation, citing_author, cited_author;
	int avg, begin_test;
	int **doc_cite_0, **doc_cite_1, **cite_topic_0, **cite_topic_1,
			**doc_topic_0, **doc_topic_1, **doc_citation;

	int *doc_array, *word_array, *citation_array, *test_docs, *to_test;
	double termGamma = 0.0, vocabGamma = 0.0, loglikelihood = 0.0,
			loglikelihood1 = 0.0, loglikelihood2 = 0.0;
	char base[20], *train, *test;
	int j, k;
	int **term_doc;
	int **citation_matrix;
	int *topic_array;
	etime();
	//srand( time(NULL));
	if (argc == 1) {
		fprintf(
				stderr,
				"usage: %s T iter alpha beta ro alpha_lambda_theta alpha_lambda_phi seed begine_test_on_iteration\n",
				argv[0]);
		exit(-1);
	}
	T = atoi(argv[1]);
	ITER = atoi(argv[2]);
	alpha = 50.0 / (double) T;
	beta = atof(argv[3]);
	ro = atof(argv[4]);
	//alpha_lambda_theta = atof(argv[5]);
	//alpha_lambda_phi = atof(argv[6]);
	seed = atoi(argv[5]);
	begin_test = atoi(argv[6]);
	int percent_test = atoi(argv[7]);
	//int model = atoi(argv[10]);
	eval_mode = 0;
	seed = 2 * seed + 1;
	assert(T > 0);
	assert(ITER > 0);
	assert(alpha > 0);
	assert(beta > 0);
	assert(seed > 0);
	char tmpbuff[10];
	int count_test = 0;
	double acc, final_val;
	int count = 0;
	cout << "Testing for ALT..." << endl;
	sprintf(tmpbuff, "test%d", eval_iter + 1);
	train = "data";
	N = countntot(train);
	w = util::ivec(N);
	d = util::ivec(N);
	z = util::ivec(N);
	c = util::ivec(N);
	s = util::ivec(N);// three possible values: -1=not in context, 0=context and citing emits the word, 1=cited emits the word
	word_probs = util::dvec(W);
	read_dwc(train, d, w, c, &D, &W, &C);
	seedMT(seed);
	Wbeta = W * beta;
	Kalpha = T * alpha;

	probs = util::dvec(T);
	doc_author *map_doc_author = new doc_author[D];
	doc_author *map_doc_cited_author = new doc_author[D];
	for (int i = 0; i < D; i++) {
		map_doc_author[i].count = 0;
		map_doc_cited_author[i].count = 0;
	}
	read_authors("author_doc", map_doc_author, &A);
	read_authors("doc_cited_authors", map_doc_cited_author, &CA);
	Cro = A * ro;

	doc_citation = util::imat(A, A);
	int *to_sample = util::ivec(N);
	test_docs = util::ivec(A);
	cited_doc = util::ivec(A);
	to_test = util::ivec(N);
	term_doc = util::imat(A, W);
	ztot = util::ivec(T);
	cztot = util::ivec(T);
	wp = util::imat(W, T);
	dp = util::imat(A, T);
	cp = util::imat(A, T);
	dl = util::ivec(A);
	c = util::ivec(N);
	d = util::ivec(N);
	w = util::ivec(N);
	z = util::ivec(N);
	int *citations_2_pred = util::ivec(A);
	read_dwc(train, d, w, c, &D, &W, &C);
	randomassignment_alt(N, D, T, w, d, c, z, s, wp, dp, cp, doc_cite_0,
			doc_cite_1, ztot, cztot, map_doc_author, to_sample, to_test);//
	int test_token_cts = 0;
	for (i = 0; i < N; i++)
		if (to_test[i] == 1)
			test_token_cts++;
	cout << test_token_cts << endl;
	occur_topic = util::dmat(T, W);
	topic_cite = util::dmat(T, A);
	topic_doc = util::dmat(T, A);
	int doc_seen = -1, citation_seen = -1;
	for (iter = 0; iter <= ITER; iter++) {
		termGamma = 0.0;
		vocabGamma = 0.0;
		loglikelihood = 0.0;
		loglikelihood1 = 0.0;
		loglikelihood2 = 0.0;
		for (i = 0; i < N; i++) {
			word = w[i];
			doc = d[i];//d[i] holds the author for this token
			citing_author = d[i];
			if (to_sample[i] == 0 || to_test[i] == 1)
				continue;
			citation = c[i];
			cited_author = c[i];
			doc_array = dp[citing_author];
			word_array = wp[word];
			if (citation >= 0) {
				citation_array = cp[cited_author];
				cited_doc[cited_author] = 1;
			}
			tt = z[i];
			word_array[tt]--;
			if (doc_array[tt] == 0) {
				cout << "here" << endl;
				exit(0);
			}
			doc_array[tt]--;
			if (citation >= 0) {
				if (doc_seen != citing_author || citation_seen != cited_author) {
					citation_array[tt]--;
					cztot[tt]--;
				}
			}
			ztot[tt]--;
			if (citation >= 0) {
				//sample s_i and topic
				if (doc_seen != citing_author || citation_seen != cited_author) {
					for (t = 0, Z = 0.0; t < T; t++) {
						probs[t] = ((citation_array[t] + ro) * (doc_array[t]
								+ alpha)) / ((cztot[t] + Cro));
						Z += probs[t];
					}
				} /*else {
				 for (t = 0, Z = 0.0; t < T; t++) {
				 probs[t] = ((doc_array[t] + alpha) * (word_array[t]
				 + beta)) / ((ztot[t] + Wbeta));
				 Z += probs[t];
				 }
				 }*/
			} //else {
			for (t = 0, Z = 0.0; t < T; t++) {
				probs[t] = ((doc_array[t] + alpha) * (word_array[t] + beta))
						/ ((ztot[t] + Wbeta));
				Z += probs[t];
			}
			//}
			//rand_int = rand() % 1000 + 1;
			rand_double = mtrand();//(double) rand_int / 1000.0;
			U = Z * rand_double;
			currprob = probs[0];
			ttt = 0;
			while (U > currprob) {
				ttt++;
				currprob += probs[ttt];
			}
			z[i] = ttt;
			word_array[ttt]++;
			doc_array[ttt]++;
			if (citation >= 0) {
				if (doc_seen != citing_author || citation_seen != cited_author) {
					citation_array[ttt]++;
					cztot[ttt]++;
				}
			}
			ztot[ttt]++;
			doc_seen = citing_author;
			citation_seen = cited_author;
		}
		if (iter >= begin_test) {
			for (i = 0; i < N; i++) {
				word = w[i];
				doc = d[i];//d[i] holds the author for this token
				citation = c[i];
				cited_author = c[i];
				if (to_test[i] == 0)
					continue;
				doc_array = dp[doc];
				word_array = wp[word];
				if (iter == begin_test) {
					term_doc[doc][word]++;
					dl[doc]++;
					//cout<<"here"<<endl;
					if (citation >= 0 && cited_doc[cited_author] == 1) {
						doc_citation[doc][citation] = 1;
						test_docs[doc] = 1;
					}
				}
				if (citation >= 0) {
					citation_array = cp[citation];
					//cited_doc[citation] = 1;//this might be changed to select citations only from cited docs

				}
				tt = z[i];
				word_array[tt]--;
				doc_array[tt]--;
				ztot[tt]--;
				for (t = 0, Z = 0.0; t < T; t++) {
					probs[t]
							= ((doc_array[t] + alpha) * (word_array[t] + beta))
									/ ((ztot[t] + Wbeta));
					Z += probs[t];
				}
				//rand_int = rand() % 1000 + 1;
				rand_double = mtrand();
				U = Z * rand_double;
				currprob = probs[0];
				ttt = 0;
				while (U > currprob) {
					ttt++;
					currprob += probs[ttt];
				}
				z[i] = ttt;
				word_array[ttt]++;
				doc_array[ttt]++;
				ztot[ttt]++;
			}

		}
		if (iter == begin_test) {

			for (i = 0; i < A; i++) {
				for (j = 0; j < A; j++) {
					if (doc_citation[i][j] == 1)
						citations_2_pred[i]++;
				}
			}
		}
		if (iter < begin_test && iter % 20 == 0) {
			cout << "iter:" << iter << ":in " << etime() << " sec." << endl;
		}
		if (iter >= begin_test) {
			for (t = 0; t < T; t++)
				for (i = 0; i < W; i++)
					occur_topic[t][i] = ((wp[i][t] + beta) / (ztot[t] + Wbeta));
			for (t = 0; t < T; t++)
				for (i = 0; i < A; i++)

					topic_cite[t][i] = ((cp[i][t] + ro) / (cztot[t] + Cro));
			for (t = 0; t < T; t++)
				for (i = 0; i < A; i++)
					if (dl[i] > 0)
						topic_doc[t][i] = ((dp[i][t] + ro) / (dl[i] + Kalpha));
			final_val = 0.0;
			for (i = 0; i < A; i++) {
				if (dl[i] > 0) {
					for (j = 0; j < W; j++) {
						acc = 0.0;
						for (t = 0; t < T; t++) {
							acc += (topic_doc[t][i] * occur_topic[t][j]);
						}
						final_val += term_doc[i][j] * log(acc);
					}
				}
			}
			double perplexity = 0.0;
			/*if (iter == begin_test)
			 for (i = 0; i < N; i++)
			 if (to_test[i] == 1)
			 test_token_cts++;*/
			perplexity = exp(-1 * final_val / test_token_cts);
			cout << "cite-lda:iter:" << iter << ":in " << etime() << " sec"
					<< ":word-ll:" << final_val << ":Perplexity:" << perplexity
					<< endl;
			calculate_rank(test_docs, cited_doc, doc_citation,
												citations_2_pred, topic_cite, topic_doc, W, A, T);
		}
	}
	util::free_ivec(w);
	util::free_ivec(d);
	util::free_ivec(z);
	util::free_ivec(c);
	util::free_ivec(dl);
	util::free_imat(dp);
	util::free_imat(wp);
	util::free_imat(term_doc);
	util::free_dmat(occur_topic);
	util::free_dmat(topic_cite);
	util::free_dmat(topic_doc);

	return 0;
}
