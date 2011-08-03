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
	CountData tmpcount("/Users/pickrell/Desktop/rosenburg_model_20pop2_phylopop_in.gz", &p);
	GraphState2 state(&tmpcount, &p);
	state.set_graph_from_file("phylopop_treein");

	/*
	CountData tmpcount("testin_counts.gz", &p);
	for( map<string, int>::iterator it = tmpcount.pop2id.begin(); it != tmpcount.pop2id.end() ; it++){
		cout << it->first << " "<< tmpcount.mean_hzy[it->second]<< " "<< tmpcount.mean_ninds[it->second] << " "<< tmpcount.id2nsnp[it->second]<< "\n";
	}
	GraphState2 state(&tmpcount, &p);

	state.add_pop();
	state.iterate_hillclimb();
	state.add_pop();
	state.iterate_hillclimb();
	state.add_pop();
	state.iterate_hillclimb();

	state.print_sigma();
	cout << "llik: "<< state.llik()<<"\n";
	cout << state.tree->get_newick_format() <<"\n";


	vector<Graph::vertex_descriptor> inorder = state.tree->get_inorder_traversal(6);
	state.tree->print();

	int i = 0;
	for(vector<Graph::vertex_descriptor>::iterator it = inorder.begin(); it != inorder.end(); it++){
		cout << i << " "<< state.tree->g[*it].index << "\n";
		i++;
	}
	state.many_global_hillclimb();
	//state.tree->global_rearrange(inorder[10], inorder[8]);
	//state.tree->print();
*/

	/*
	cout << state.tree->does_mig_exist(inorder[10], inorder[6]) <<"\n";
	cout << "legal? "<< state.tree->is_legal_migration(inorder[10], inorder[6]) <<"\n";
	state.tree->add_mig_edge(inorder[10], inorder[6]);
	cout << state.tree->does_mig_exist(inorder[10], inorder[6]) <<"\n";
	state.tree->print();
	state.set_branches_ls();
	state.tree->print();
	cout << "llik: "<< state.llik()<<"\n";
	state.print_sigma();
*/
	//state.set_graph("(((((pop0:0.1,pop1:0.1):0.1,pop2:0.1):0.1,pop3:0.1):0.1,pop4:0.1):0.1,pop5:0.1);");
	//state.tree->print();
	//state.set_branches_ls();
	//cout << state.tree->get_newick_format() <<"\n";
	//state.print_sigma();
	//cout << state.llik() << "\n";

	//state.set_graph("((pop0:0.1,(pop1:0.1,pop2:0.1):0.1):0.1,(pop5:0.1,(pop3:0.1,pop4:0.1):0.1));");
	//state.set_branches_ls();
	//cout << "\n";
	//state.print_sigma();
	//cout << state.llik() << "\n";
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
