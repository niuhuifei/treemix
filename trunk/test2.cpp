/*
 * test2.cpp
 *
 *  Created on: Apr 7, 2011
 *      Author: pickrell
 */

#include "CountData.h"
#include "PopGraph.h"
#include "GraphState2.h"
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
	//p.bias_correct = false;
	CountData tmpcount("testin_counts.gz", 1);
	tmpcount.read_alfreqs("testin_freqs3.gz");
	tmpcount.print_cov("test_cov.gz");
	GraphState2 state(&tmpcount, &p);
	//state.print_sigma();
	state.add_pop();
	state.add_pop();
	state.add_pop();
	//cout << state.tree->get_newick_format() <<"\n";

	//state.tree->print();
	//set<Graph::vertex_descriptor> tset = state.tree->get_root_adj();
	//for (set<Graph::vertex_descriptor>::iterator it = tset.begin(); it != tset.end(); it++){
	//	cout << state.tree->g[*it].index << "\n";
	//}

	//state.set_graph("(((((pop0:0.1,pop1:0.1):0.1,pop2:0.1):0.1,pop3:0.1):0.1,pop4:0.1):0.1,pop5:0.1);");
	//state.set_graph("((pop0:0.1,(pop1:0.1,pop2:0.1):0.1):0.1,(pop5:0.1,(pop3:0.1,pop4:0.1):0.1));");
	//state.tree->print();
	//state.set_branches_ls();
	cout << state.tree->get_newick_format() <<"\n";
	/*
	vector<string> tmpnames;
	tmpnames.push_back("pop1"); tmpnames.push_back("pop2");tmpnames.push_back("pop3");
	PopGraph tmp(tmpnames);

	tmp.print();

	vector<Graph::vertex_descriptor> inord = tmp.get_inorder_traversal_noroot(3);
	cout << "here\n";
	for (int i = 0; i < inord.size(); i++){
		cout << tmp.g[ inord[i]].index << " "<< tmp.g[ inord[i]].name << "\n";
	}

	Graph::vertex_descriptor n = tmp.add_tip(inord[1], "pop4");
	tmp.print();

	tmp.remove_tip(n);
	cout << "\n";
	tmp.print();
	*/
    return 0;
}
