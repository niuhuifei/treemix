/*
 * TreeMix.cpp
 *
 *  Created on: Apr 12, 2011
 *      Author: pickrell
 */

#include "Settings.hpp"
#include "GraphState3.h"
#include "PhyloPop_params.h"

string infile;
string outstem = "TreeMix";
int nthread = 1;
int norm_type = 0;
void printopts(){
	cout << "\nTreeMix v.0.1 \n by JKP\n\n";
	cout << "Options:\n";
	cout << "-i [file name] input file\n";
	cout << "-o [stem] output stem (will be [stem].treeout.gz, [stem].cov.gz, [stem].modelcov.gz)\n";
	cout << "-scale Rescale the allele frequencies by sqrt( p (1-p) )\n";
	cout << "-k [int] number of SNPs per block for estimation of covariance matrix (1)\n";
	cout << "-global Do a round of global rearrangements after adding all populations\n";
	cout << "-tf [file name] Read the tree topology from a file, rather than estimating it\n";
	cout << "-m [int] number of migration edges to add (0)\n";
	cout << "-mn [int] depth of subtree to use in migration optimization (3)\n";
	cout << "-root [string] comma-delimited list of populations to set on one side of the root (for migration)\n";
	cout << "-g [vertices file name] [edges file name] read the graph from a previous TreeMix run\n";
	cout << "-quick do fast optimization of migration weights\n";
	cout << "\n";
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
    if (cmdline.HasSwitch("-g"))	{
      	p.vfile = cmdline.GetArgument("-g", 0);
      	p.efile = cmdline.GetArgument("-g", 1);
      	p.read_graph = true;
      }
    if (cmdline.HasSwitch("-arcsin")) p.alfreq_scaling = 1;
    if (cmdline.HasSwitch("-scale")) p.alfreq_scaling = 3;
    if (cmdline.HasSwitch("-quick")) p.quick = true;
    if (cmdline.HasSwitch("-global")) p.global = true;
    if (cmdline.HasSwitch("-k"))	p.window_size = atoi(cmdline.GetArgument("-k", 0).c_str());
    if (cmdline.HasSwitch("-m"))	p.nmig = atoi(cmdline.GetArgument("-m", 0).c_str());
    if (cmdline.HasSwitch("-mn"))	p.m_neigh = atoi(cmdline.GetArgument("-mn", 0).c_str());
    if (cmdline.HasSwitch("-root")) {
    	p.set_root = true;
    	p.root = cmdline.GetArgument("-root", 0);
    }
    string treefile = outstem+".treeout.gz";
    string covfile = outstem+".cov.gz";
    string modelcovfile = outstem+".modelcov.gz";
    string cov_sefile = outstem+".covse.gz";
    string llikfile = outstem+".llik";
    ofstream likout(llikfile.c_str());

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
    else if (p.read_graph){
    	state.set_graph(p.vfile, p.efile);
    	cout << "Set tree to: "<< state.tree->get_newick_format() << "\n";
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
    	state.set_branches_ls_wmig();
    }
    if (p.set_root) state.place_root(p.root);
    likout << "Tree likelihood: "<< state.llik() << "\n";
    for (int i = 0; i < p.nmig; i++){
    	double st = state.llik();
    	pair<bool, pair<int, int> > add = state.add_mig_targeted();
    	//cout << "Added? "<< add.first << "\n";
    	if (add.first == true) state.iterate_mig_hillclimb_and_optimweight(add.second);
    	//cout << "Optim 1\n";
    	if (p.quick){
    		p.quick = false;
    		state.optimize_weights();
    		p.quick = true;
    	}
    	else{
    		state.optimize_weights();
    	}
    	//cout << "Optim 2\n";
    	state.trim_mig();
    	state.flip_mig();
    	state.trim_mig();
    	//likout << st << " "<< state.current_llik << "\n";
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
     		double w = state.tree->g[*it].weight;

        		double lk = state.llik();
        		treeout << state.tree->g[*it].weight<< " ";
        		state.tree->g[*it].weight = 0;
        		state.set_branches_ls_wmig();

        		double lk0 = state.llik();
        		treeout << lk0 << " "<< lk << " ";
        		state.tree->g[*it].weight = w;
        		state.set_branches_ls_wmig();

        		Graph::vertex_descriptor p1 = source( *it, state.tree->g);
        		p1 = state.tree->get_child_node_mig(p1);
        		Graph::vertex_descriptor p2 = target(*it, state.tree->g);
        		treeout << state.tree->get_newick_format(p1) << " ";
        		treeout << state.tree->get_newick_format(p2) << "\n";
    	}
		it++;
    }
   // cout << state.llik_wishart() << "\n";
    //cout << state.tree->get_newick_format();
   // state.set_branches_ls_wmig();
    //cout << state.tree->get_newick_format();
    //state.tree->print();
    state.print_trimmed(outstem);
    state.print_sigma_cor(modelcovfile);
  //  state.tree->print();
  //  map<int, Graph::vertex_descriptor> test = state.tree->index2vertex();
  //  Graph::vertex_descriptor moz = test[28];
  //  set<Graph::edge_descriptor> mig = state.tree->get_in_mig_edges(moz);
  //  Graph::edge_descriptor mige = *(mig.begin());
  //  for(double m = 0; m < 0.95; m+=0.05){
  //  	state.tree->g[mige].weight = m;
  //  	state.set_branches_ls_wmig();
   // 	cout << m << "\n";
  //  	state.tree->print();
  //  	cout << state.llik() << "\n";
  // }

	return 0;
}
