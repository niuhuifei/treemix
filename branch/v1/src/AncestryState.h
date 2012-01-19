/*
 * AncestryState.h
 *
 *  Created on: Nov 11, 2011
 *      Author: pickrell
 */

#ifndef ANCESTRYSTATE_H_
#define ANCESTRYSTATE_H_

#include "Settings.hpp"
#include "CountData.h"
#include "PhyloPop_params.h";

class AncestryState{
public:
	AncestryState();
	AncestryState(CountData*, CountData*, PhyloPop_params*);
	vector<gsl_matrix *> covs;
	void input_covs(string);
	CountData* counts;
	CountData* counts2;
	PhyloPop_params* params;
	double c;
	double llik();

	void set_sigmacor_from_sigma();
	void set_sigma(int);
	void print_ancestry_llik(string);
	gsl_matrix * sigma;
	gsl_matrix * sigma_cor;
};


#endif /* ANCESTRYSTATE_H_ */
