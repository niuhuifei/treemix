/*
 * GraphState3.cpp
 *
 *  Created on: Aug 12, 2011
 *      Author: pickrell
 */

#include "GraphState3.h"

GraphState3::GraphState3( GraphState2 * g, PhyloPop_params * p){
	state = g;
	params = p;
	scratch = new GraphState2( g->countdata, p);
}


bool GraphState3::add_mig_targeted(){
	pair<string, string> target = state->get_max_resid();
	cout << "Targeting migration to "<< target.first << " "<< target.second << "\n"; cout.flush();

	//Identify neighborhood
	map<string, Graph::vertex_descriptor> tips = state->tree->get_tips(state->tree->root);
	Graph::vertex_descriptor n1 = state->get_neighborhood( tips[ target.first ] );
	Graph::vertex_descriptor n2 = state->get_neighborhood( tips[ target.second ] );
	//state->tree->print();
	//cout << state->tree->g[n1].index << " "<< state->tree->g[n2].index << "\n";
	string subtree = state->tree->get_newick_subtrees(n1, n2);
	cout << "Using subtree:\n";
	cout << subtree <<"\n";
	scratch->set_graph_from_string(subtree);

	// get best migration
	pair< pair<bool, bool>, pair<double, pair<int, int> > > m = scratch->add_mig_targeted( target.first, target.second);


	// add it to the main tree if it's legal
	if ( m.first.first == true){
		int nup1 = m.second.second.first;
		int nup2 = m.second.second.second;
		bool rev = m.first.second;
		Graph::vertex_descriptor a1 = tips[ target.first];
		Graph::vertex_descriptor a2 = tips[target.second];
		int i =0;
		while (i < nup1){
			a1 = state->tree->get_parent_node(a1).first;
			i++;
		}
		i = 0;
		while (i < nup2){
			a2 = state->tree->get_parent_node(a2).first;
			i++;
		}
		Graph::edge_descriptor e;
		if (rev ) e= state->tree->add_mig_edge(a2, a1);
		else e= state->tree->add_mig_edge(a1, a2);
		state->tree->g[e].weight = m.second.first;
		state->set_branches_ls_wmig();
	}


}
