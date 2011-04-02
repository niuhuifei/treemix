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

class State{
public:

	/*
	 *  Initializes tree to newick string, count data)
	 */
	State(string, CountData*);

	/*
	 * Copy constructor
	 */

	State(const State&);


	PhyloPop_Tree::Tree<PhyloPop_Tree::NodeData>* tree;
	gsl_matrix* thetas;
	gsl_matrix* means;
	gsl_matrix* sigma;
	CountData* countdata;
	void compute_sigma();
	void print_sigma();
	double llik();

};
#endif /* STATE_H_ */
