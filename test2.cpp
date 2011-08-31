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
	p.window_size = 200;
	///CountData tmpcount("/Users/pickrell/Desktop/rosenburg_model_20pop2_phylopop_in.gz", &p);
	//GraphState2 state(&tmpcount, &p);
	//state.set_graph_from_file("phylopop_treein");
	CountData tmpcount("/Users/pickrell/projects/phylopop/data/harvard_hgdp_yoruba_asc_noprimate_phylopop.gz", &p);
	GraphState2 state(&tmpcount, &p);
	state.set_graph("/Users/pickrell/projects/phylopop/data/yoruba_asc_4mig_afroot.vertices.gz", "/Users/pickrell/projects/phylopop/data/yoruba_asc_4mig_afroot.edges.gz");
	//state.add_mig_targeted();
	//state.tree->place_root("Denisova,Neander,San,BiakaPygmy,MbutiPygmy,BantuKenya,BantuSouthAfrica,Mandenka,Yoruba");
	//state.tree->print("test");
	//cout << state.llik() <<"\n";
	//state.rearrange(42, 58);
	//cout << state.llik() <<"\n";
	//state.add_mig(25, 42);
	//cout << state.llik() <<"\n";
	//state.add_mig(58, 42);
	//cout << state.llik() <<"\n";
	//state.optimize_weights();
	//state.add_mig(70, 42);
	//cout << state.llik() <<"\n";

	//state.rearrange(42, 55);
	//cout << state.llik() << "\n";
	//state.add_mig(25, 42);
	//cout << state.llik() << "\n";
	//state.add_mig(25, 54);
	//cout << state.llik() << "\n";

	//state.rearrange(42, 70);
	//cout << state.llik() << "\n";
	//state.add_mig(25, 42);
	//cout << state.llik() << "\n";
	//state.add_mig(42, 54);
	//cout << state.llik() << "\n";

	//state.add_mig(42, 54);
	//cout << state.llik() << "\n";
	//state.add_mig(70, 42);
	//cout << state.llik() << "\n";

	//state.add_mig(13, 25);
	//state.add_mig(34, 13);
	//cout << state.llik() << "\n";
	state.optimize_weights();
	cout << state.llik() << "\n";
	cout << state.tree->get_newick_format() << "\n";


	//CountData tmpcount("testin_counts.gz", &p);
//	for( map<string, int>::iterator it = tmpcount.pop2id.begin(); it != tmpcount.pop2id.end() ; it++){
	//	cout << it->first << " "<< tmpcount.mean_hzy[it->second]<< " "<< tmpcount.mean_ninds[it->second] << " "<< tmpcount.id2nsnp[it->second]<< "\n";
//	}
//	GraphState2 state(&tmpcount, &p);
//	cout.precision(8);
//	state.add_pop();
//	state.iterate_hillclimb();
//	state.add_pop();
//	state.iterate_hillclimb();
//	state.add_pop();
//	state.iterate_hillclimb();

//	state.tree->print();
//	vector<Graph::vertex_descriptor> inorder = state.tree->get_inorder_traversal(6);
//	for (int i = 0; i < inorder.size(); i++)	cout << i << " "<< state.tree->g[ inorder[i]].index<<"\n";
//	state.add_mig_targeted();
//	state.tree->print();
//	state.local_hillclimb_wmig(51);
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
