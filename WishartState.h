/*
 * WishartState.h
 *
 *  Created on: Apr 19, 2011
 *      Author: pickrell
 */

#ifndef WISHARTSTATE_H_
#define WISHARTSTATE_H_

#include "BinaryTree.h"
#include "Settings.hpp"
#include "CountData.h"
#include "MCMC_params.h"

class WishartState{
public:

	/*
	 *  Initializes tree to newick string, count data)
	 */
	WishartState(string, CountData*, MCMC_params*);

	/*
	 * Copy constructor
	 */

	WishartState(const WishartState&);


	// pointer to the tree data structure
	PhyloPop_Tree::Tree<PhyloPop_Tree::NodeData>* tree;

	gsl_matrix* sigma;
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
	void update_thetas(gsl_rng*);
	double update_theta_snp(gsl_rng*, int);

	// compute the covariance matrix from the tree
	void compute_sigma();
	void print_sigma();
	void set_sigma_inv(); //set the inverse

	//compute the log-likelihood of the data given the tree
	double llik();
	double llik_snp(int);
	double llik_snp_old(int);
	double dens_mvnorm(const gsl_vector*, const gsl_vector*);

	void print_state(ogzstream&, ogzstream&);
private:

};


#endif /* WISHARTSTATE_H_ */
