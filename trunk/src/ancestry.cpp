/*
 * ancestry.cpp
 *
 *  Created on: Nov 11, 2011
 *      Author: pickrell
 */





#include "AncestryState.h"
#include "PhyloPop_params.h"

int main (void)
{

	PhyloPop_params p;
	p.window_size = 500;
	p.alfreq_scaling = 4;

	CountData counts1("/Users/pickrell/projects/ancestry/afam/data/single_afam_nopop.treemix.gz", &p);
	double f2 = counts1.calculate_f2(1,2);

	PhyloPop_params p2;
	p2.window_size = 1000;
	p2.smooth_scale = sqrt( (double) counts1.nsnp / (double) p2.window_size);
	CountData counts2("/Users/pickrell/projects/ancestry/afam/data/single_afam_nopop.treemix.gz", &p2);
	AncestryState state(&counts2, &counts1, &p2);
	state.print_ancestry_llik("test");
	return 0;
}

