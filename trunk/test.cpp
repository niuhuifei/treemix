/*
 * test.cpp
 *
 *  Created on: Mar 29, 2011
 *      Author: pickrell
 */


#include "CountData.h"
#include "PopGraph.h"
#include "GraphState2.h"
#include "PhyloPop_params.h"

int main(){
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_ranlxs2;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, 100);


};
