/*
 * MCMC_params.h
 *
 *  Created on: Apr 5, 2011
 *      Author: pickrell
 */

#ifndef MCMC_PARAMS_H_
#define MCMC_PARAMS_H_
#include "Settings.hpp"

class MCMC_params{
public:
	MCMC_params();
	double epsilon; //for proposal of new trees
	double lambda; //paramters for beta prior on ancestral allele frequencies
	double s2; //parameter for normal distribution for proposal of new ancestral frequencies
	double s3; //parameter for normal distribution for proposal of new thetas
	double B; //maximum distance from any tip to the root
	double f; //fraction of iterations that are global updates
	int burnin; //number of iterations to burn in
	int total; //total number of iterations
	int samp; //sample from the chain every samp iterations
	int psamp; //print the likelihood every psamp iterations
	int nthread; //for multithreading; how many threads?

};

#endif /* MCMC_PARAMS_H_ */
