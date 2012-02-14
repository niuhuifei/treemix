/*
 * test_easia.cpp
 *
 *  Created on: Nov 14, 2011
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
	p.window_size = 500;
	CountData counts("/Users/pickrell/projects/treemix/easia/data/easia+africa+russian+french.treemix.gz", &p);
	counts.print_cov("easia_wmigs.cov.gz");
	counts.print_cov_var("easia_wmigs.covse.gz");
	GraphState2 state(&counts, &p);
	state.set_graph("/Users/pickrell/projects/treemix/easia/data/treemix_runs/migration_runs/easia+africa+russian+french2_2mig.vertices.gz", "/Users/pickrell/projects/treemix/easia/data/treemix_runs/migration_runs/easia+africa+russian+french2_2mig.edges.gz");
	state.add_mig(1616, 535);
	state.add_mig(1616, 1095);
	state.add_mig(355, 2655);
	state.add_mig(355, 603);
	state.optimize_fracs();
	state.optimize_weights();
	state.tree->print("easia_wmigs");



    return 0;
}
