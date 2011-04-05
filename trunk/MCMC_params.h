/*
 * MCMC_params.h
 *
 *  Created on: Apr 5, 2011
 *      Author: pickrell
 */

#ifndef MCMC_PARAMS_H_
#define MCMC_PARAMS_H_

class MCMC_params{
public:
	MCMC_params();
	double epsilon; //for proposal of new trees
	double lambda; //paramters for beta prior on ancestral allele frequencies
	double s2; //parameter for normal distribution for proposal of new ancestral frequencies
	double B; //maximum distance from any tip to the root
	int burnin; //number of iterations to burn in
	int total; //total number of iterations
	int samp; //sample from the chain every samp iterations
};

#endif /* MCMC_PARAMS_H_ */
