/*
 * MCMC.cpp
 *
 *  Created on: Apr 5, 2011
 *      Author: pickrell
 */

#include "MCMC.h"

MCMC::MCMC(State* s, MCMC_params* p){
	state = s;
	param =p;
}


void MCMC::run(gsl_rng* r){
	cout <<"Burn in of "<< param->burnin << "\n";
	for (int i = 0; i < param->burnin; i++){
		state->update(r);
	}
	cout << "Sampling\n";
	for (int i = 0; i < param->total; i++){
		if (i % param->samp == 0){
			//do something to sample the state
		}
		state->update(r);
	}
}
