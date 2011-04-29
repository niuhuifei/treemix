/*
 * mean_fst.cpp
 *
 *  Created on: Apr 28, 2011
 *      Author: pickrell
 */
#include "CountData.h"
#include "gzstream.h"
#include "CmdLine.h"

string infile;
string outfile = "fst.out.gz";

int main(int argc, char *argv[]){
    CCmdLine cmdline;
    if (cmdline.SplitLine(argc, argv) < 1){
    	exit(1);
    }
    if (cmdline.HasSwitch("-i")) infile = cmdline.GetArgument("-i", 0).c_str();
    else{
    	exit(1);
    }
    if (cmdline.HasSwitch("-o"))	outfile = cmdline.GetArgument("-o", 0).c_str();
    CountData c(infile);
    c.print_fst(outfile);
    return 0;
}
