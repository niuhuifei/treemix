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
		if (i % param->psamp == 0) cout <<  "Burn-in step "<< i << " llik: "<< state->current_lik << "\n"; cout.flush();
		run_one(r);
	}
	param->epsilon = param->epsilon/5;
	cout << "Sampling\n";
	for (int i = 0; i < param->total; i++){
		if (i % param->samp == 0) state->print_state(treefile);
		if (i % param->psamp == 0) cout <<  "Sampling step "<< i << " llik: "<< state->current_lik << "\n";
		run_one(r);
	}
}

void MCMC_graph::run_one(gsl_rng* r){
	double ran = gsl_rng_uniform(r);
	if (ran < param->f)	state->update_tree(r);
	else state->local_update_tree(r);
}
