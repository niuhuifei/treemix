/*
 * MCMC_graph.cpp
 *
 *  Created on: Apr 25, 2011
 *      Author: pickrell
 */


#include "MCMC_graph.h"
#include "mvnorm.h"

MCMC_graph::MCMC_graph(GraphState* s, MCMC_params* p){
	state = s;
	param =p;
}


void MCMC_graph::run(gsl_rng* r, ogzstream& treefile){
	cout <<"Burn in of "<< param->burnin << "\n";
	for (int i = 0; i < param->burnin; i++){
		//cout << i  <<"\n";
		if (i % param->samp == 0)			state->print_state(treefile);
		if (i % param->psamp == 0) cout <<  "Burn-in step "<< i << " llik: "<< state->current_lik << "\n";
		state->update_tree(r);
	}
	param->epsilon = param->epsilon/10;
	cout << "Sampling\n";
	for (int i = 0; i < param->total; i++){
		if (i % param->samp == 0) state->print_state(treefile);
		if (i % param->psamp == 0) cout <<  "Sampling step "<< i << " llik: "<< state->current_lik << "\n";
		state->update_tree(r);
	}
}
