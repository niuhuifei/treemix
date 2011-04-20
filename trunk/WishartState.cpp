/*
 * WishartState.cpp
 *
 *  Created on: Apr 19, 2011
 *      Author: pickrell
 */

#include "WishartState.h"


WishartState::WishartState(string newick, CountData* counts, MCMC_params* p){
	countdata = counts;
	params =p;
	tree = new PhyloPop_Tree::BinaryTree<PhyloPop_Tree::NodeData>(newick, countdata->pop2id);
	sigma = gsl_matrix_alloc(counts->npop, counts->npop);

	gsl_matrix_set_zero(sigma);

}



WishartState::WishartState(const WishartState& oldstate){
	tree = oldstate.tree->copy();
	countdata = oldstate.countdata;
	params = oldstate.params;
	sigma = gsl_matrix_alloc(oldstate.countdata->npop, oldstate.countdata->npop);
	gsl_matrix_memcpy(sigma, oldstate.sigma);


}


