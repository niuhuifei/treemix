/*
 * scan_snp.cpp
 *
 *  Created on: Dec 12, 2011
 *      Author: pickrell
 */

#include "Settings.hpp"
#include "SNP.h"
#include "CmdLine.h"
#include "CountData.h"

string infile;
string ffile;
string outfile = "scanout.gz";
int main(int argc, char *argv[]){
    CCmdLine cmdline;
    if (cmdline.SplitLine(argc, argv) < 1){
    	exit(1);
    }
    if (cmdline.HasSwitch("-i")) infile = cmdline.GetArgument("-i", 0).c_str();
    else{
    	exit(1);
    }
    if (cmdline.HasSwitch("-f")) ffile = cmdline.GetArgument("-f", 0).c_str();
     else{
     	exit(1);
    }
    if (cmdline.HasSwitch("-o")) outfile = cmdline.GetArgument("-o", 0).c_str();

    ogzstream out(outfile.c_str());

    PhyloPop_params p;
    p.snpinfo = true;
    p.window_size = 500;
    CountData counts(infile, &p);
    F2_matrix f2(ffile);

	map<string,double> trim;
	double sumtrim = 0;
	for ( map<string, int>::iterator it = counts.pop2id.begin(); it!= counts.pop2id.end(); it++){
		int id = it->second;
		string pop = it->first;
		double meanhzy = counts.mean_hzy.find(id)->second;
		double mean_n = counts.mean_ninds.find(id)->second;
		double t = meanhzy / (4.0* mean_n);
		sumtrim+= t;
		//cout << pop  << " "<< t << " "<< meanhzy << " "<< mean_n << "\n";
		trim.insert(make_pair(pop, t));
	}

    for (int i = 0; i < counts.nsnp; i++){
    	pair<vector<string>, vector<double> > freqs = counts.get_freqs(i);
    	SNP tmpsnp(freqs.second, freqs.first, &f2, &trim);
    	tmpsnp.optimize_lambda();
    	out << counts.rss.at(i) << " "<< counts.chr.at(i) << " "<< counts.pos.at(i) << " "<< counts.a1.at(i) << " "<< counts.a2.at(i) << " ";
    	out << tmpsnp.lambda << " "<< tmpsnp.ss() << "\n"; out.flush();
    }




}