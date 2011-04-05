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


	PhyloPop_Tree::Tree<PhyloPop_Tree::NodeData>* tree;

	gsl_matrix* thetas;
	gsl_vector* means;
	gsl_matrix* sigma;
	CountData* countdata;
	MCMC_params* params;

	void update(gsl_rng*);
	void compute_sigma();
	void print_sigma();
	double llik();
	double llik_snp(int);
	void update_tree(gsl_rng*);
	vector<PhyloPop_Tree::iterator<PhyloPop_Tree::NodeData> > propose_tree(gsl_rng*);
	void update_means(gsl_rng*);
	void update_mean(gsl_rng*, int);

	//double epsilon; //for proposal of new trees
	//double lambda; //parameter for beta prior on ancestral allele frequences

};
#endif /* STATE_H_ */
