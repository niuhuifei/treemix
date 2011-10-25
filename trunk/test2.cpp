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
	p.window_size = 1000;
	p.alfreq_scaling = 0;
	ofstream tmpfile("llikout");
	CountData tmpcount("/Users/pickrell/projects/treemix/sims/rosenburg_model_20pop2_phylopop_in.gz", &p);
	GraphState2 state(&tmpcount, &p);

	state.set_graph("((pop5:0.00259392,((pop3:0.00303462,(pop1:0.0036395,pop2:0.00312606):0.000436326):0.000279838,pop4:0.00283215):0.000274779):0.000135986,((pop7:0.00230621,(pop8:0.002108,((((((pop14:0.00111402,((pop16:0.000744062,(pop17:0.00052396,((pop20:0.000360078,pop19:0.000223979):0.000297038,pop18:0.000416602):0.000294024):0.000276486):0.000264836,pop15:0.000971684):0.000265374):0.000280506,pop13:0.00132408):0.000285774,pop12:0.00150098):0.000263486,pop11:0.00155312):0.000291052,pop10:0.00176539):0.000316614,pop9:0.00193484):0.000247762):0.000325461):0.000386744,pop6:0.00246966):0.000135986);");
	state.print_sigma_cor("incov.gz");
	vector<string> pops = state.allpopnames;
	cout << "here\n"; cout.flush();
	CountData tmpcount2(&tmpcount, pops, state.sigma_cor, &p, r);
	GraphState2 randstate(&tmpcount2, &p);
	randstate.set_graph(&state);
	/*
	for (int i = 0; i < 90; i++){
		stringstream ss;
		ss<< "/Users/pickrell/projects/treemix/sims/sim_mvn/simfreqs" << i << ".gz";
		string file  = ss.str();
		CountData tmpcount(file, &p);
	//	tmpcount.print_alfreqs("alfreqs.gz");
	//	tmpcount.print_cov_var("covvar.gz");
		GraphState2 state(&tmpcount, &p);
		//cout << "setting graph\n";cout.flush();
		state.set_graph("(pop1:0.00532549,((pop3:0.0183318,(pop4:0.0178017,(pop5:0.0163895,(pop6:0.0155497,(pop7:0.0146479,(pop8:0.0134163,(pop9:0.0125101,(pop10:0.0115314,(pop11:0.0105312,(pop12:0.00968735,(pop13:0.0087489,(pop15:0.00765831,pop14:0.00734194):0.00104139):0.00111911):0.000921365):0.00100286):0.00106475):0.000900554):0.00114427):0.000932313):0.00108868):0.000831278):0.000982886):0.00122326,pop2:0.0193423):0.0162424);");

		state.set_branches_ls_wmig();
		//cout << "done\n"; cout.flush();
		double l0 = state.llik_wishart();
		map<string, Graph::vertex_descriptor> tips = state.tree->get_tips(state.tree->root);
		int i1 = state.tree->g[ tips["pop2"] ].index;
		int i2 = state.tree->g[ tips["pop6"] ].index;
		state.add_mig(i1, i2);
		tmpfile << l0 << " "<< state.llik_wishart() << "\n";
		//state.compute_sigma();
		//state.print_sigma();
		//state.print_sigma_cor("testcor.gz");
		//state.tree->print("testtree");
	}
*/
	//state.set_graph("/Users/pickrell/projects/treemix/hgdp+new/san_asc_tree_archaic_root.vertices.gz", "/Users/pickrell/projects/treemix/hgdp+new/san_asc_tree_archaic_root.edges.gz");
	//cout << state.llik() << " "<< state.llik_wishart() << "\n";
	//state.rearrange(8055, 536);
	//cout << state.llik() << " "<< state.llik_wishart() << "\n";
	//state.rearrange(535, 8055);
	//cout << state.llik() << " "<< state.llik_wishart() << "\n";
//	state.llik_wishart();
//	tmpcount.print_cov_var2("covvar2.gz");
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
