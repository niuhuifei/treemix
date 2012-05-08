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
	p.window_size = 500;
	//p.smooth_scale = 0;

	//CountData counts("/Users/pickrell/projects/treemix/hgdp_reich_asc/harvard_hgdp_yoruba_asc_noprimate_phylopop.gz", &p);
	CountData counts("/Users/pickrell/projects/treemix/sims/test_4pop.gz", &p);
	GraphState2 state(&counts, &p);
	state.add_pop();
	cout << "\n\n\n";
	//for (int i = 0;i < 17; i++) state.add_pop();
	//state.set_graph("/Users/pickrell/projects/treemix/hgdp_reich_asc/trees/tree_replicates_nnls/test.vertices.gz","/Users/pickrell/projects/treemix/hgdp_reich_asc/trees/tree_replicates_nnls/test.edges.gz");
	//state.set_branches_ls();

	/*while (state.current_llik <= -DBL_MAX){
   		cout << "RESCALING\n"; cout.flush();
   		p.smooth_scale = p.smooth_scale *2;
   		state.current_llik = state.llik();
   	}
	double llk = state.llik();
	map<string, Graph::vertex_descriptor> pop2tip = state.tree->get_tips(state.tree->root);
	for (map<string, Graph::vertex_descriptor>::iterator it = pop2tip.begin(); it != pop2tip.end(); it++){
		int i1 = state.tree->g[it->second].index;
		for (map<string, Graph::vertex_descriptor>::iterator it2 = pop2tip.begin(); it2 != pop2tip.end(); it2++){
			int i2 = state.tree->g[it2->second].index;
			pair<bool, Graph::edge_descriptor> test = state.add_mig(i1, i2);
			if (test.first){
				Graph::edge_descriptor e = test.second;
				cout << i1 <<  " "<< i2 << " "<< it->first << " "<< it2->first << " "<< state.tree->g[e].weight << " "<< llk << " "<< state.llik() << " "<< p.smooth_scale << "\n";
				state.tree->remove_mig_edge(e);
			}
		}
	}
*/
	cout << "here\n"; cout.flush();
	state.set_branches_ls_wmig();

	cout << "here2\n"; cout.flush();
	state.set_branches_ls();
	return 0;
}
