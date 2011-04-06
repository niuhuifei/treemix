/*
 * State.h
 *
 *  Created on: Mar 31, 2011
 *      Author: pickrell
 */

#ifndef STATE_H_
#define STATE_H_


#include "BinaryTree.h"
#include "Settings.hpp"
#include "CountData.h"
#include "MCMC_params.h"

class State{
public:

	/*
	 *  Initializes tree to newick string, count data)
	 */
	State(string, CountData*, MCMC_params*);

	/*
	 * Copy constructor
	 */

	State(const State&);

	/*
	 * read in thetas instead of estimating them (for testing from simulations)
	 */
	void read_thetas(string, int, int);


	// pointer to the tree data structure
	PhyloPop_Tree::Tree<PhyloPop_Tree::NodeData>* tree;

	gsl_matrix* thetas;
	gsl_vector* means;
	gsl_matrix* sigma;
	CountData* countdata; //pointer to the data
	MCMC_params* params; //pointer to the parameters for updates

	// workhorse function is update(), does all the Metropolis updates using the functions that follow it
	void update(gsl_rng*);
	void update_tree(gsl_rng*);
	vector<PhyloPop_Tree::iterator<PhyloPop_Tree::NodeData> > propose_tree(gsl_rng*);
	void update_means(gsl_rng*);
	void update_mean(gsl_rng*, int);

	// compute the covariance matrix from the tree
	void compute_sigma();
	void print_sigma();

	//compute the log-likelihood of the data given the tree
	double llik();
	double llik_snp(int);
	void print_state(ogzstream&, ogzstream&);

};
#endif /* STATE_H_ */
