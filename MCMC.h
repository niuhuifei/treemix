/*
 * MCMC.h
 *
 *  Created on: Apr 5, 2011
 *      Author: pickrell
 */

#ifndef MCMC_H_
#define MCMC_H_
#include "Settings.hpp"
#include "State.h"
class MCMC{
public:
	MCMC(State*, MCMC_params*);
	State* state;
	MCMC_params* param;
	void run(gsl_rng*, ogzstream&, ogzstream&);
};

#endif /* MCMC_H_ */
