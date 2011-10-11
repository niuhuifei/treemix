/*
 * test2.cpp
 *
 *  Created on: Apr 7, 2011
 *      Author: pickrell
 */

#include "CountData.h"
#include "PopGraph.h"
#include "GraphState2.h"
#include "GraphState3.h"
#include "PhyloPop_params.h"

int main (void)
{

	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_ranlxs2;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, 100);

	PhyloPop_params p;
	p.window_size = 1000;
	p.alfreq_scaling = 0;
	CountData tmpcount("/Users/pickrell/projects/phylopop/rosenburg_model_20pop2_phylopop_in.gz", &p);
	tmpcount.print_alfreqs("alfreqs.gz");
	tmpcount.print_cov_var("covvar.gz");
	GraphState2 state(&tmpcount, &p);
	state.llik_wishart();
	tmpcount.print_cov_var2("covvar2.gz");
	//state.set_graph_from_file("phylopop_treein");
	//state.tree->print();

	/*
	CountData tmpcount("/Users/pickrell/projects/phylopop/data/hgdp_freqs.phylopop_in.gz", &p);
	GraphState2 state(&tmpcount, &p);
	state.set_graph("/Users/pickrell/projects/phylopop/data/hgdp/noscale/hgdp_tree_global_noscale_afroot.vertices.gz", "/Users/pickrell/projects/phylopop/data/hgdp/noscale/hgdp_tree_global_noscale_afroot.edges.gz");
	//state.place_root("San,BiakaPygmy,MbutiPygmy,Yoruba,Mandenka,BantuSouthAfrica,BantuKenya");
	cout << state.llik() <<"\n"; cout.flush();
	//state.add_mig(4, 4892);
	//cout << state.llik() << "\n";
	//state.add_mig(536, 2);
	//cout << state.llik() << "\n";
	state.add_mig(2, 751);
	cout << state.llik() << "\n";
	state.tree->print("before_flip");
	state.flip_mig();
	//state.iterate_mig_hillclimb_and_optimweight(make_pair(2, 751));
	state.tree->print("after_flip");
	cout << state.tree->get_newick_format() << "\n";
*/


/*

	CountData tmpcount("testin_counts.gz", &p);
//	for( map<string, int>::iterator it = tmpcount.pop2id.begin(); it != tmpcount.pop2id.end() ; it++){
	//	cout << it->first << " "<< tmpcount.mean_hzy[it->second]<< " "<< tmpcount.mean_ninds[it->second] << " "<< tmpcount.id2nsnp[it->second]<< "\n";
//	}
	GraphState2 state(&tmpcount, &p);
//	cout.precision(8);
	state.add_pop();
	state.iterate_hillclimb();
	state.add_pop();
	state.iterate_hillclimb();
	state.add_pop();
	state.iterate_hillclimb();
	Graph::edge_descriptor e = state.add_mig(3, 51);
	state.tree->set_mig_frac(e, 0.3);
	e = state.add_mig(3, 31);
	state.tree->set_mig_frac(e, 0.7);

	state.tree->set_mig_frac(e, 0.1);

	state.tree->set_mig_frac(e, 0.7);
	//state.tree->g[source(e, state.tree->g)].mig_frac = 0.1;
	e = state.add_mig(3, 15);

	state.tree->set_mig_frac(e, 0.05);
	state.tree->print("before");
	cout << "here\n";
	state.tree->set_mig_frac(e, 0.95);
	state.tree->print("test");
//	vector<Graph::vertex_descriptor> inorder = state.tree->get_inorder_traversal(6);
//	for (int i = 0; i < inorder.size(); i++)	cout << i << " "<< state.tree->g[ inorder[i]].index<<"\n";
//	state.add_mig_targeted();
//	state.tree->print();
//	state.add_mig(4, 31);
//	state.local_hillclimb_wmig(31);
	//state.tree->place_root("pop6,pop3");
	//state.tree->print();
	//string newick =  state.tree->get_newick_subtrees( inorder[0], inorder[5] );
	//state.set_graph_from_string(newick);
	//cout << state.tree->get_newick_format() << "\n";

	/*

	cout << state.llik()<<"\n";
	state.print_sigma();

	pair<Graph::edge_iterator, Graph::edge_iterator> eds = edges(state.tree->g);
	Graph::edge_iterator it = eds.first;
	while (it != eds.second){
		if ( state.tree->g[*it].is_mig){
			cout << state.tree->g[*it].weight<< " "<< state.tree->g[*it].len << " ";
			Graph::vertex_descriptor p1 = source( *it, state.tree->g);
			p1 = state.tree->get_child_node_mig(p1);
			Graph::vertex_descriptor p2 = target(*it, state.tree->g);
			cout << state.tree->get_newick_format(p1) << " ";
			cout << state.tree->get_newick_format(p2) << "\n";
		}
		it++;
	}
	*/

    return 0;
}
