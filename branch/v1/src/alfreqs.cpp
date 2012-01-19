/*
 * alfreqs.cpp
 *
 *  Created on: Nov 10, 2011
 *      Author: pickrell
 */







#include "GraphState2.h"
#include "PhyloPop_params.h"

int main (void)
{

	PhyloPop_params p;
	p.window_size = 500;
	p.alfreq_scaling = 4;

	ofstream out("outfile");
	CountData counts("/Users/pickrell/projects/treemix/hgdp+new/san_asc.treemix_in.gz", &p);
	counts.print_alfreqs("alfreqs.gz");
	return 0;
}


