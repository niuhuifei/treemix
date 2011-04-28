/*
 * PhyloPop_graph.cpp
 *
 *  Created on: Apr 25, 2011
 *      Author: pickrell
 */


#include "Settings.hpp"
#include "MCMC_graph.h"

string infile;
string outstem = "PhyloPop";
int seed = 200;
int nthread = 1;
bool printscatter = false;
int scaling = 1;
int burn = 20000;
int sample = 40000;
double epsilon = 0.002;

void printopts(){
	cout << "\nPhyloPop v.0.0 \n by JKP\n\n";
	cout << "Options:\n";
	cout << "-i input file\n";
	cout << "-o output stem (will be [stem].treeout.gz and [stem].meanout.gz)\n";
	cout << "-s print scatter matrix (will be [stem].scatter.gz)\n";
	cout << "-burn number of burn-in iterations of the MCMC (default 20,000)\n";
	cout << "-samp number of sampling iterations of the MCMC (default 40,000)\n";
	cout << "-ep parameter for changing topology during burn-in (default 0.002)\n";
	cout << "-sqrt use square root transformation of allele frequencies (default is arcsin)\n";
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
    if (cmdline.HasSwitch("-s"))	printscatter = true;
    if (cmdline.HasSwitch("-burn"))	burn = atoi(cmdline.GetArgument("-burn", 0).c_str());
    if (cmdline.HasSwitch("-samp"))	sample = atoi(cmdline.GetArgument("-samp", 0).c_str());
    if (cmdline.HasSwitch("-ep"))	epsilon = atof(cmdline.GetArgument("-ep", 0).c_str());
    if (cmdline.HasSwitch("-sqrt"))	scaling = 2;

    string treefile = outstem+".treeout.gz";
    ogzstream treeout(treefile.c_str());

    //cout << "Reading data and calculating scatter matrix.."; cout.flush();
    CountData counts(infile, scaling);
    //cout << "done\n"; cout.flush();




    // MCMC parameters
    MCMC_params p;
    p.epsilon = epsilon;
    p.burnin = burn;
    p.total = sample;
    p.nthread = nthread;

    string pops = counts.get_pops();
    GraphState initstate(pops,  &counts, &p);


    //start random number generator
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_ranlxs2;
    r = gsl_rng_alloc(T);
    if (seed > 0)  gsl_rng_set(r,(int)seed);
    else{
    	seed = (int)time(0);
    	gsl_rng_set(r, seed);
    }

    // initialization
    initstate.init_tree(r);
    //initstate.tree->randomize_tree(r);
    cout << initstate.tree->get_newick_format() <<"\n";
    initstate.compute_sigma();
    initstate.init();
    //MCMC
    MCMC_graph mcmc(&initstate, &p);
    mcmc.run(r, treeout);
    if (printscatter) counts.print_scatter(outstem+".scatter.gz");
	return 0;
}
