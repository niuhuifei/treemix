/*
 * mig_est_4pop.cpp
 *
 *  Created on: Oct 25, 2011
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
	p.restrict_pop = true;
	p.pops2use = 3;
	p.smooth_lik = false;
	CountData tmpcount("/Users/pickrell/projects/treemix/hgdp/sa_fr_pap_bed.gz", &p);
	GraphState2 state(&tmpcount, &p);
	state.compute_sigma();
	//state.print_sigma();
	cout << state.tree->get_newick_format() << "\n";
	PhyloPop_params p2;
	p2.window_size = 1000;
	p2.alfreq_scaling = 0;
	p2.restrict_pop = true;
	p2.pops2use = 4;
	p2.smooth_lik = true;
	p2.smooth_scale = 1000;
	p2.tau = 0.001;
	CountData tmpcount2("/Users/pickrell/projects/treemix/hgdp/sa_fr_pap_bed.gz", &p2);
	string outgroup = tmpcount2.get_pop_in_index(0);
	string ingroup = tmpcount2.get_pop_in_index(1);
	string ingroup2 = tmpcount2.get_pop_in_index(2);
	string testgroup = tmpcount2.get_pop_in_index(3);
	p.smooth_scale = 1000;
	p.smooth_lik = true;
	state.set_countdata( &tmpcount2 );
	state.place_root(outgroup);
	cout << state.tree->get_newick_format() << "\n";
	state.add_pop(ingroup, testgroup);
	state.compute_sigma();
	state.set_sigmacor_from_sigma();

	map<string, Graph::vertex_descriptor> p2v = state.tree->get_tips(state.tree->root);
	Graph::vertex_descriptor t = p2v[testgroup];
	Graph::vertex_descriptor o = p2v[outgroup];
	Graph::vertex_descriptor c = p2v[ingroup];
	Graph::vertex_descriptor o2 = state.tree->get_parent_node(c).first;
	state.tree->print("test");
	double l1 = state.tree->get_parent_node(c).second;
	Graph::in_edge_iterator it = in_edges(t, state.tree->g).first;
	Graph::edge_descriptor e = *it;
	state.golden_section_edge(e, -20, log(state.tree->g[e].len), 1, p2.tau);
	p.smooth_scale = 100;
	state.golden_section_edge(e, -20, log(state.tree->g[e].len), 1, p2.tau);
	p.smooth_scale = 10;
	state.golden_section_edge(e, -20, log(state.tree->g[e].len), 1, p2.tau);
	state.print_sigma();
	Graph::edge_descriptor m = state.add_mig_noopt( state.tree->g[o2].index, state.tree->g[ p2v[testgroup] ].index );

	Graph::in_edge_iterator ei = in_edges(c, state.tree->g).first;
	Graph::edge_descriptor c1 = *ei;
	Graph::in_edge_iterator ei2 = in_edges(state.tree->get_parent_node(c).first, state.tree->g).first;
	Graph::edge_descriptor c2 = *ei2;
	ofstream outfile("outfile");
	for(double f = 0; f < 1; f+=0.01){
		state.tree->g[m].weight = f;
		for(double f2 = 0; f2 < 1; f2+= 0.1){
			state.tree->set_mig_frac(m, f2);
			for(double f3 = 0; f3 < 1; f3+= 0.1){
				state.tree->g[c1].len = l1*f3;
				state.tree->g[c2].len = l1*(1-f3);

				state.compute_sigma();
				state.set_sigmacor_from_sigma();
				state.golden_section_edge(e, -20, log(state.tree->g[e].len), 1, p2.tau);
				outfile << f << " " << f2 << " " << f3 << " "<< state.tree->g[e].len << " "<< state.llik() << "\n";
			}
		}
	}
	state.print_sigma();
    return 0;
}
