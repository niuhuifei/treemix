/*
 * test2.cpp
 *
 *  Created on: Apr 7, 2011
 *      Author: pickrell
 */

#include "GraphState2.h"
#include "PhyloPop_params.h"

int main (void)
{

	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_ranlxs2;
	r = gsl_rng_alloc(T);
	int seed = (int) time(0);
	gsl_rng_set(r, seed);


	PhyloPop_params p;
	p.window_size = 100;
	p.smooth_scale = 0;
	p.f2 = false;
	p.alfreq_scaling = 0;

	CountData counts("/Users/pickrell/projects/treemix/sims/test_20pop.gz", &p);
	GraphState2 state(&counts, &p);
	state.add_pop();
	state.add_pop();


	state.set_branches_ls();
	state.set_branches_ls_wmig();
	return 0;
}
