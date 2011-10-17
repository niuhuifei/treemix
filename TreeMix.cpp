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
	cout << "\nTreeMix v.0.2 \n\n";
	cout << "Options:\n";
	cout << "-i [file name] input file\n";
	cout << "-o [stem] output stem (will be [stem].treeout.gz, [stem].cov.gz, [stem].modelcov.gz)\n";
	cout << "-scale Rescale the allele frequencies by sqrt( p (1-p) )\n";
	cout << "-k [int] number of SNPs per block for estimation of covariance matrix (1)\n";
	cout << "-global Do a round of global rearrangements after adding all populations\n";
	cout << "-tf [file name] Read the tree topology from a file, rather than estimating it\n";
	cout << "-m [int] number of migration edges to add (0)\n";
	cout << "-root [string] comma-delimited list of populations to set on one side of the root (for migration)\n";
	cout << "-g [vertices file name] [edges file name] read the graph from a previous TreeMix run\n";
	cout << "-quick do fast optimization of migration weights\n";
	cout << "-nofrac Force migration nodes to fall in the middle of their branches\n";
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
    if (cmdline.HasSwitch("-nofrac")) p.nofrac = true;
    if (cmdline.HasSwitch("-scale")) p.alfreq_scaling = 3;
    if (cmdline.HasSwitch("-quick")) p.quick = true;
    if (cmdline.HasSwitch("-global")) p.global = true;
    if (cmdline.HasSwitch("-k"))	p.window_size = atoi(cmdline.GetArgument("-k", 0).c_str());
    if (cmdline.HasSwitch("-m"))	p.nmig = atoi(cmdline.GetArgument("-m", 0).c_str());
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
    cout.precision(10);
    if (p.readtree) {
    	state.set_graph_from_file(p.treefile);

    	//state.iterate_hillclimb();
    }
    else if (p.read_graph){
    	state.set_graph(p.vfile, p.efile);
    	cout << "Set tree to: "<< state.tree->get_newick_format() << "\n";
    	cout << "ln(lk) = " << state.current_llik <<"\n";
    }

    while (!p.readtree && counts.npop > state.current_npops){
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

    	if (add.first == true) state.iterate_mig_hillclimb_and_optimweight(add.second);

    	//cout << "Optim 1\n";
    	if (p.quick){
    		p.quick = false;
    		if (!p.nofrac) state.optimize_fracs();
    		state.optimize_weights();

    		p.quick = true;
    	}
    	else{

    		if (!p.nofrac) state.optimize_fracs();
    		state.optimize_weights();

    	}
    	//state.tree->print("before_trim1");
    	state.trim_mig();
    	//state.tree->print("before_flip");
    	state.flip_mig();
    	//state.tree->print("before_trim2");
    	state.trim_mig();
    	state.set_branches_ls_wmig();
    	cout << "ln(likelihood):" << state.current_llik << "\n";
    }

    treeout << state.tree->get_newick_format() << "\n";
    treeout << state.get_trimmed_newick() << "\n";
    state.current_llik_w = state.llik_wishart();
    state.optimize_fracs_wish();
    state.optimize_weights_wish();
    pair<Graph::edge_iterator, Graph::edge_iterator> eds = edges(state.tree->g);
    Graph::edge_iterator it = eds.first;
    p.smooth_lik = false;
    while (it != eds.second){
    	if ( state.tree->g[*it].is_mig){
     		double w = state.tree->g[*it].weight;

        		double lk = state.llik();
        		double wlk = state.llik_wishart();
        		treeout << state.tree->g[*it].weight<< " ";
        		state.tree->g[*it].weight = 0;
        		state.set_branches_ls_wmig();
        		//state.optimize_fracs(*it);
        		//state.optimize_weights(*it);
        		double lk0 = state.llik();
        		double wlk0 = state.llik_wishart();
        		treeout << lk0 << " "<< lk << " ";
        		treeout << wlk0 << " "<< wlk << " ";
        		double pval = 1-gsl_cdf_chisq_P(-2 *(wlk0-wlk), 2 );
        		treeout << pval << " ";
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

    state.print_trimmed(outstem);
    state.print_sigma_cor(modelcovfile);


	return 0;
}
