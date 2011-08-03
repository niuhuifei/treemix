/*
 * PhyloPop.cpp
 *
 *  Created on: Apr 12, 2011
 *      Author: pickrell
 */

#include "Settings.hpp"
#include "GraphState2.h"
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
    if (p.readtree) state.set_graph_from_file(p.treefile);


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
    treeout << state.tree->get_newick_format() << "\n";
    treeout << state.get_trimmed_newick() << "\n";
    state.print_sigma_cor(modelcovfile);
	return 0;
}
