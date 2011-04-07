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
	gsl_vector* snp_liks; //store the log likelihoods of individual snps
	double current_lik; //store the sum of the individual snp log(liks) ie. the total likelihood
	CountData* countdata; //pointer to the data
	MCMC_params* params; //pointer to the parameters for updates

	//initialize the likelihood
	void init_liks();

	// workhorse function is update(), does all the Metropolis updates using the functions that follow it
	void update(gsl_rng*);
	void update_tree(gsl_rng*);
	vector<PhyloPop_Tree::iterator<PhyloPop_Tree::NodeData> > propose_tree(gsl_rng*);
	void update_means(gsl_rng*);
	double update_mean(gsl_rng*, int);

	// compute the covariance matrix from the tree
	void compute_sigma();
	void print_sigma();
	void set_sigma_inv(); //set the inverse

	//compute the log-likelihood of the data given the tree
	double llik();
	double llik_snp(int);
	double dens_mvnorm(const gsl_vector*, const gsl_vector*);

	void print_state(ogzstream&, ogzstream&);
private:
	//store the inversion of sigma and the determinant of sigma
	gsl_matrix* winv;
	double ax;
};
#endif /* STATE_H_ */
