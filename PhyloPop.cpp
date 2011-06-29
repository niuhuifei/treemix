/*
 * PhyloPop.cpp
 *
 *  Created on: Apr 12, 2011
 *      Author: pickrell
 */

#include "Settings.hpp"
#include "GraphState2.h"

string infile;
string outstem = "PhyloPop";
int seed = 200;
int nthread = 1;

void printopts(){
	cout << "\nPhyloPop v.0.0 \n by JKP\n\n";
	cout << "Options:\n";
	cout << "-i input file\n";
	cout << "-o output stem (will be [stem].treeout.gz, [stem].cov.gz)\n";
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

    string treefile = outstem+".treeout.gz";
    string covfile = outstem+".cov.gz";

    ogzstream treeout(treefile.c_str());
    CountData counts(infile, 1);
    counts.print_cov(covfile);
    GraphState2 state(&counts);
    while (counts.npop > state.current_npops){
    	//cout << counts.npop << " "<< state.current_npops << "\n";
    	state.add_pop();

    }
    cout << state.tree->get_newick_format() << "\n";
	return 0;
}
