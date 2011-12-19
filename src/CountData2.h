/*
 * CountData2.h
 *
 *  Created on: Jun 16, 2011
 *      Author: pickrell
 */

#ifndef COUNTDATA2_H_
#define COUNTDATA2_H_
#include "Settings.hpp"


class CountData2{
public:
	CountData2(string, string, gsl_rng*);
	void read_counts(string);
	void set_alfreqs();
	void scale_alfreqs();
	void center_alfreqs();
	void randomize_order(gsl_rng*);
	void set_scatter();
	void set_cov();
	void process_scatter();
	void print_scatter(string);
	void print_cov(string);

	string meanpop;

	map<string, int> pop2id;
	vector<string> popsnomean;
	gsl_vector *means;
	gsl_matrix *alfreqs, *scatter, *cov;
	double scatter_det, scatter_gamma, cov_var;
	vector<vector<pair<int, int> > > allele_counts;
	int npop, nsnp;
};

#endif /* COUNTDATA2_H_ */
