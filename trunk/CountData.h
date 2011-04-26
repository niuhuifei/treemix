/*
 * CountData.h
 *
 *  Created on: Apr 1, 2011
 *      Author: pickrell
 */

#ifndef COUNTDATA_H_
#define COUNTDATA_H_
#include "Settings.hpp"



class CountData{
public:
	CountData(string);
	void read_counts(string);
	map<string, int> pop2id;
	vector<vector<pair<int, int> > > allele_counts;
	int npop, nsnp;
	string get_pops();
	gsl_matrix *alfreqs, *scatter;
	void set_alfreqs();
	void scale_alfreqs();
	void set_scatter();
	void process_scatter();
	void read_scatter(string); //for debugging
	void print_scatter(string);
	double scatter_det, scatter_gamma;
};

#endif /* COUNTDATA_H_ */
