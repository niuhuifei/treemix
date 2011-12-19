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
	p.smooth_scale = 1;

	//CountData tmpcount("20pop_testin.gz", &p);
	//CountData tmpcount("/Users/pickrell/projects/treemix/sims/rosenburg_model_20pop2_phylopop_in.gz", &p);
	CountData tmpcount("/Users/pickrell/projects/treemix/hgdp_reich_asc/harvard_hgdp_yoruba_asc_noprimate_phylopop.gz", &p);
	GraphState2 state(&tmpcount, &p);
	state.set_graph("/Users/pickrell/projects/treemix/hgdp_reich_asc/trees/tree_replicates_nnls/mig_replicates/yoruba3_6mig.vertices.gz", "/Users/pickrell/projects/treemix/hgdp_reich_asc/trees/tree_replicates_nnls/mig_replicates/yoruba3_6mig.edges.gz");
	cout << state.llik() << "\n";
	state.iterate_local_hillclimb_wmig_all();
	cout << state.llik() << "\n";
	state.tree->print("test3");
	//state.set_graph("/Users/pickrell/projects/treemix/hgdp_reich_asc/trees/tree_replicates_nnls/yoruba1.vertices.gz", "/Users/pickrell/projects/treemix/hgdp_reich_asc/trees/tree_replicates_nnls/yoruba1.edges.gz");
	//state.initialize_migupdate();
	//cout << state.llik() << "\n";
	//state.set_branches_ls_f2();
	//state.tree->print("test");
	//cout << state.llik() << "\n";
	//state.initialize_migupdate();
	//state.set_branches_ls_f2_precompute();
	//state.tree->print("test2");
	//cout << state.llik() << "\n";
	//state.print_X();
	//state.tree->print();
	//state.all_try_movemig();
	//cout << state.llik() << "\n";
	//state.tree->print("test2");
	//state.optimize_weights();


	//Graph::edge_descriptor e = state.add_mig(171, 211);
	//pair<double, double> se = state.calculate_se(e);
	//double test = se.first/ se.second;
	//double pval = 1-gsl_cdf_gaussian_P(test, 1);
	//cout << se.first << " "<< se.second << " "<< test << " " << pval << "\n";


	//vector<Graph::edge_descriptor> mig_edges = state.tree->get_mig_edges();
	//Graph::edge_descriptor e = mig_edges[0];

	//state.initialize_migupdate();
	//state.optimize_weight_quick(e);

	//cout << state.llik()<< "\n"; cout.flush();

	//Graph::vertex_descriptor t = target(e, state.tree->g);
	//Graph::vertex_descriptor s = source(e, state.tree->g);
	//Graph::vertex_descriptor pa = state.tree->get_parent_node(s).first;
	//state.tree->remove_mig_edge(e);
	//state.add_mig(916,state.tree->g[t].index );
	//state.add_mig(state.tree->g[ pa ].index, state.tree->g[t].index);

	//state.optimize_weights_quick();
	//cout << state.all_try_movemig() <<"\n";
	//cout << state.llik()<< "\n"; cout.flush();
	//cout << state.all_try_changedir() <<"\n" ;
	//cout << state.llik() << "\n"; cout.flush();
	//state.tree->print("test");


/*
	CountData tmpcount("4pop_testin.gz", &p);
	GraphState2 state(&tmpcount, &p);

	state.add_mig(4, 1);
	state.add_mig(3, 1);
	state.tree->print("test0");
  */
    return 0;
}
