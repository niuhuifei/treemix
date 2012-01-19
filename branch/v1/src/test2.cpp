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
	CountData tmpcount("/Users/pickrell/projects/treemix/hgdp_reich_asc/remove_admix/yoruba_asc_filter_admix.treemix_in.gz", &p);
	GraphState2 state(&tmpcount, &p);
	state.set_graph("/Users/pickrell/projects/treemix/hgdp_reich_asc/remove_admix/migration/y_filter1_1mig.vertices.gz", "/Users/pickrell/projects/treemix/hgdp_reich_asc/remove_admix/migration/y_filter1_1mig.edges.gz");
	cout << state.llik() << "\n";
	state.print_sigma_cor("test1.modelcov.gz");
	//state.add_mig_targeted_f2();
	Graph::edge_descriptor e = state.add_mig(2655, 2);
	state.print_sigma_cor("test2.modelcov.gz");
	state.tree->print("test");
	state.initialize_migupdate();
	set<pair<double, set<Graph::edge_descriptor> > > t = state.popname2paths["Cambodian"];
	for (set<pair<double, set<Graph::edge_descriptor> > >::iterator it = t.begin(); it != t.end(); it++){
		cout << it->first << " w\n";
		for (set<Graph::edge_descriptor>::iterator it2 = it->second.begin(); it2!= it->second.end(); it2++){
			cout << state.tree->g[source(*it2, state.tree->g)].index << " "<< state.tree->g[target(*it2, state.tree->g)].index << " "<< state.tree->g[*it2].is_mig << " "<< state.tree->g[*it2].weight << "\n";
		}
	}
	//state.tree->g[e].weight = 0.0458838;
	//state.set_branches_ls_f2();
	//for (int i = 0; i < 50; i+=5){
	//	double tmp = (double) i/ 100.0;
	//	state.tree->g[e].weight = tmp;
	//	state.set_branches_ls_f2();
	//	cout << i << " "<< state.llik() << "\n";
	//}

	//state.tree->print("test3");
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


    return 0;
}
