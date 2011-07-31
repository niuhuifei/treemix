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
#include "PhyloPop_params.h"

class GraphState2{
public:
	GraphState2(CountData*, PhyloPop_params*);

	PhyloPop_params* params; //paramters for run
	PopGraph* tree;
	PopGraph* tree_bk;

	gsl_matrix *sigma;
	gsl_matrix *sigma_cor;
	CountData* countdata; //pointer to the data
	vector<string> allpopnames; //names of populations, will be added one at a time after the first 3
	int current_npops; //current total number of populations
	double current_llik;
	double scatter_det, scatter_gamma;
	gsl_matrix *scatter; //current scatter matrix

	//set the graph structure to a Newick string
	void set_graph(string);

	//covariance matrix
	void compute_sigma();
	void print_sigma();
	void print_sigma_cor(string);

	//local hill-climbing
	int local_hillclimb(int);
	int many_local_hillclimb();
	void iterate_hillclimb();

	//add a new population
	void add_pop();
	void process_scatter();

	//under normal model, get the max lik branch lengths
	// for a given topology by least squares
	void set_branches_ls();
	void set_branches_ls_old();

	//likelihoods
	double llik();
	double llik_normal();
	double llik_wishart();

	//migration
	void add_mig();

	//get newick string with trimmed terminal branch lengths
	string get_trimmed_newick();
};


#endif /* GRAPHSTATE2_H_ */
