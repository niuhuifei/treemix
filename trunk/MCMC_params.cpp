/*
 * MCMC_params.cpp
 *
 *  Created on: Apr 5, 2011
 *      Author: pickrell
 */
#include "MCMC_params.h"

MCMC_params::MCMC_params(){
	lambda = 1;
	s2 = 0.05;
	epsilon = 0.01;
	B = 10;
	burnin = 10000;
	total = 20000;
	samp = 10;
	psamp = 200;
}
/*

MCMC_params::MCMC_params(ogzstream& tfile, ogzstream& mfile){
	lambda = 1;
	s2 = 0.05;
	epsilon = 0.1;
	B = 10;
	burnin = 2000;
	total = 20000;
	samp = 10;
	treeoutfile = tfile;
	moutfile = mfile;
}
*/
