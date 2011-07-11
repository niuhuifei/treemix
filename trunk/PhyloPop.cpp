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
	cout << "-i input file\n";
	cout << "-o output stem (will be [stem].treeout.gz, [stem].cov.gz), [stem].modelcov.gz\n";
	cout << "-arcsin (perform the arcsin square root transformation on the allele frequencies before centering)\n";
}


int main(int argc, char *argv[]){
    CCmdLine cmdline;
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
    if (cmdline.HasSwitch("-arcsin")) norm_type = 1;
    string treefile = outstem+".treeout.gz";
    string covfile = outstem+".cov.gz";
    string modelcovfile = outstem+".modelcov.gz";
    PhyloPop_params p;
    //p.bias_correct = false;
    ogzstream treeout(treefile.c_str());
    CountData counts(infile, norm_type);
    counts.print_cov(covfile);
    GraphState2 state(&counts, &p);
    //while (5 > state.current_npops){
    while (counts.npop > state.current_npops){
    	//cout << counts.npop << " "<< state.current_npops << "\n";
    	state.add_pop();
    	state.iterate_hillclimb();
    	cout << "ln(likelihood): "<< state.current_llik << "\n";
    	cout << state.tree->get_newick_format() << "\n";
    }
    treeout << state.tree->get_newick_format() << "\n";
    state.print_sigma_cor(modelcovfile);
	return 0;
}
