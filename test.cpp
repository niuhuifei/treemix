/*
 * test.cpp
 *
 *  Created on: Mar 29, 2011
 *      Author: pickrell
 */

//#include "Tree.h"
//#include "BinaryTree.h"
#include "MCMC.h"
#include "CountData.h"
#include "WishartState.h"
#include "GraphState.h"

int main(){
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_ranlxs2;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, 100);
	string testnewick = "((YRI:0.8,LWD:0.4):0.1,(CEU:0.4,(JPT:0.2,CHB:0.3):0.2):0.5);";
	ogzstream treeout("treeout.gz");
	ogzstream mout("meanout.gz");
	CountData counts("testin_counts.gz");
	MCMC_params p;

	GraphState teststate(testnewick, &counts, &p);
	teststate.tree->print();
	//teststate.tree->randomize_tree(r);
	teststate.compute_sigma();

	for(int i = 0; i < counts.npop; i++){
		for (int j = 0; j < counts.npop; j++){
			cout << gsl_matrix_get(teststate.sigma, i, j) << " ";
		}
		cout << "\n";
	}
	cout << "llik: "<< teststate.llik() << "\n";

	//vector<vector<double> > freqs = counts.get_alfreqs();
	//for (int i = 0; i < counts.nsnp ; i++){
	//	for (int j = 0; j < counts.npop; j++){
	//		std::cerr << gsl_matrix_get(counts.alfreqs, i, j) << " ";
	//	}
	//	std::cerr << "\n";
	//}


};
