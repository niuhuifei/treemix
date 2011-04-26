/*
 * MCMC.h
 *
 *  Created on: Apr 5, 2011
 *      Author: pickrell
 */

#ifndef MCMC_H_
#define MCMC_H_

#include "Settings.hpp"
#include "WishartState.h"
class MCMC{
public:
	MCMC(WishartState*, MCMC_params*);
	WishartState* state;
	MCMC_params* param;
	void run(gsl_rng*, ogzstream&);
};

#endif /* MCMC_H_ */
