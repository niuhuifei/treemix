/*
 * GraphState.h
 *
 *  Created on: Apr 25, 2011
 *      Author: pickrell
 */

#ifndef GRAPHSTATE_H_
#define GRAPHSTATE_H_
#include "Settings.hpp"
#include "PopGraph.h"
#include "CountData.h"
#include "MCMC_params.h"

class GraphState{
public:

	/*
	 *  Initializes tree to newick string, count data)
	 */
	GraphState(string, CountData*, MCMC_params*);


	// pointer to the tree data structure
	PopGraph* tree;

	gsl_matrix *sigma;
	double current_lik; //store the current likelihood
	CountData* countdata; //pointer to the data
	MCMC_params* params; //pointer to the parameters for updates

	//initialize the likelihood
	void init();

	//update the tree
	void update_tree(gsl_rng*);
	vector<Graph::vertex_descriptor> propose_tree(gsl_rng*);

	// compute the covariance matrix from the tree
	void compute_sigma();
	void print_sigma();
	void read_sigma(string);


	//compute the log-likelihood of the data given the tree
	double llik();
	double dens_wishart();
	void print_state(ogzstream&);
};


#endif /* GRAPHSTATE_H_ */
