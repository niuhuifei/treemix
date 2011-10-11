/*
 * CountData.h
 *
 *  Created on: Apr 1, 2011
 *      Author: pickrell
 */

#ifndef COUNTDATA_H_
#define COUNTDATA_H_
#include "Settings.hpp"
#include "PhyloPop_params.h"


class CountData{
public:
	CountData(string, PhyloPop_params*);
	CountData(string); //no rescaling (ie. just the allele frequencies)
	void read_counts(string);
	map<string, int> pop2id;
	map<int, double> mean_hzy;
	map<int, double> mean_ninds;
	map<int, int> id2nsnp;
	vector<vector<pair<int, int> > > allele_counts;
	int npop, nsnp, ef_nsnp;
	string get_pops(); //in Newick format
	vector<string> list_pops(); //simple list
	gsl_matrix *alfreqs, *scatter, *cov, *cov_var, *cov_var2;
	void set_alfreqs();
	void scale_alfreqs();
	void set_scatter();
	void set_cov();
	void set_cov2();
	void process_scatter();
	void process_cov();
	void read_scatter(string); //for debugging
	void read_alfreqs(string); //for debugging, can input simulated allele frequencies
	void print_scatter(string);
	void print_alfreqs(string);
	void print_fst(string);
	double get_cov(string, string);
	double get_cov_var(string, string);
	double get_scatter(string, string);
	void print_cov(string);
	void print_cov_var(string);
	void print_cov_var2(string);
	double scatter_det, scatter_gamma;
	PhyloPop_params* params;
	gsl_matrix *U, *scatter_prime;
	int ne;
	void set_ne();
};

#endif /* COUNTDATA_H_ */
