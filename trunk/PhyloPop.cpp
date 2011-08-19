/*
 * PhyloPop.cpp
 *
 *  Created on: Apr 12, 2011
 *      Author: pickrell
 */

#include "Settings.hpp"
#include "GraphState3.h"
#include "PhyloPop_params.h"

string infile;
string outstem = "PhyloPop";
int seed = 200;
int nthread = 1;
int norm_type = 3;
void printopts(){
	cout << "\nPhyloPop v.0.0 \n by JKP\n\n";
	cout << "Options:\n";
	cout << "-i [file name] input file\n";
	cout << "-o [stem] output stem (will be [stem].treeout.gz, [stem].cov.gz), [stem].modelcov.gz\n";
	cout << "-arcsin (perform the arcsin square root transformation on the allele frequencies before centering)\n";
	cout << "-k [int] number of SNPs per block for estimation of covariance matrix (1)\n";
	cout << "-global (Do a round of global rearrangements after adding all populations)\n";
	cout << "-tf [file name] Read the tree topology from a file, rather than estimating it\n";
	cout << "-m [int] number of migration edges to add (0)\n";
	cout << "-mn [int] depth of subtree to use in migration optimization (3)\n";
}


int main(int argc, char *argv[]){
    CCmdLine cmdline;
    PhyloPop_params p;
    if (cmdline.SplitLine(argc, argv) < 1){
    	printopts();
    	exit(1);
    }
    if (cmdline.HasSwitch("-i")) infile = cmdline.GetArgument("-i", 0).c_str();
    else{
    	printopts();
    	exit(1);
    }
    if (cmdline.HasSwitch("-o"))	outstem = cmdline.GetArgument("-o", 0).c_str();
    if (cmdline.HasSwitch("-tf"))	{
    	p.treefile = cmdline.GetArgument("-tf", 0).c_str();
    	p.readtree = true;
    }
    if (cmdline.HasSwitch("-arcsin")) p.alfreq_scaling = 1;
    if (cmdline.HasSwitch("-global")) p.global = true;
    if (cmdline.HasSwitch("-k"))	p.window_size = atoi(cmdline.GetArgument("-k", 0).c_str());
    if (cmdline.HasSwitch("-m"))	p.nmig = atoi(cmdline.GetArgument("-m", 0).c_str());
    if (cmdline.HasSwitch("-mn"))	p.m_neigh = atoi(cmdline.GetArgument("-mn", 0).c_str());
    string treefile = outstem+".treeout.gz";
    string covfile = outstem+".cov.gz";
    string modelcovfile = outstem+".modelcov.gz";
    string cov_sefile = outstem+".covse.gz";

    //p.bias_correct = false;
    ogzstream treeout(treefile.c_str());
    CountData counts(infile, &p);
    counts.print_cov(covfile);
    counts.print_cov_var(cov_sefile);
    GraphState2 state(&counts, &p);
    //GraphState3 s3(&state, &p);
    cout.precision(10);
    if (p.readtree) {
    	state.set_graph_from_file(p.treefile);
    	//state.iterate_hillclimb();
    }

    while (!p.readtree && counts.npop > state.current_npops){
    	//cout << counts.npop << " "<< state.current_npops << "\n";
    	state.add_pop();
    	state.iterate_hillclimb();
    	cout << "ln(likelihood): "<< state.current_llik << "\n";
    	cout << state.tree->get_newick_format() << "\n";
    }
    if (p.global){
    	cout << "Testing global rearrangements\n"; cout.flush();
    	state.iterate_global_hillclimb();
    }
    for (int i = 0; i < p.nmig; i++){
    	pair<bool, int> add = state.add_mig_targeted();
    	//cout << "Added? "<< add.first << "\n";
    	if (add.first == true) state.iterate_mig_hillclimb_and_optimweight(add.second);
    	cout << "ln(likelihood):" << state.current_llik << "\n";
    }

    //cout << state.tree->get_newick_format() << "\n";
    //state.set_branches_ls_wmig();
    //cout << state.tree->get_newick_format() << "\n";
    treeout << state.tree->get_newick_format() << "\n";
    treeout << state.get_trimmed_newick() << "\n";

    pair<Graph::edge_iterator, Graph::edge_iterator> eds = edges(state.tree->g);
    Graph::edge_iterator it = eds.first;
    while (it != eds.second){
    	if ( state.tree->g[*it].is_mig){
    		treeout << state.tree->g[*it].weight<< " "<< state.tree->g[*it].len << " ";
    		Graph::vertex_descriptor p1 = source( *it, state.tree->g);
    		p1 = state.tree->get_child_node_mig(p1);
    		Graph::vertex_descriptor p2 = target(*it, state.tree->g);
    		treeout << state.tree->get_newick_format(p1) << " ";
    		treeout << state.tree->get_newick_format(p2) << "\n";
    	}
		it++;
    }

    state.print_sigma_cor(modelcovfile);
	return 0;
}
