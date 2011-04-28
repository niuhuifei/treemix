/*
 * MCMC_graph.h
 *
 *  Created on: Apr 25, 2011
 *      Author: pickrell
 */

#ifndef MCMC_GRAPH_H_
#define MCMC_GRAPH_H_

#include "Settings.hpp"
#include "GraphState.h"
class MCMC_graph{
public:
	MCMC_graph(GraphState*, MCMC_params*);
	GraphState* state;
	MCMC_params* param;
	void run(gsl_rng*, ogzstream&);
	void run_one(gsl_rng*);
};


#endif /* MCMC_GRAPH_H_ */
