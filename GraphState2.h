/*
 * GraphState2.h
 *
 *  Created on: Jun 28, 2011
 *      Author: pickrell
 */

#ifndef GRAPHSTATE2_H_
#define GRAPHSTATE2_H_

#include "Settings.hpp"
#include "PopGraph.h"
#include "CountData.h"


class GraphState2{
public:
	GraphState2(CountData*);

	PopGraph* tree;
	PopGraph* tree_bk;

	gsl_matrix *sigma;
	CountData* countdata; //pointer to the data
	vector<string> allpopnames;
	int current_npops;
	double current_llik;

	// covariance matrix
	void compute_sigma();
	void print_sigma();

	//local hill=climbing
	void local_hillclimb(int);

	//add a new population
	void add_pop();

	//under normal model, get the max lik branch lengths
	// for a given topology by least squares
	void set_branches_ls();

	//likelihoods
	double llik_normal();
};


#endif /* GRAPHSTATE2_H_ */
