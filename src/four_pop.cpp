/*
 * four_pop.cpp
 *
 *  Created on: Nov 1, 2011
 *      Author: pickrell
 */



#include "GraphState2.h"
#include "PhyloPop_params.h"

int main (void)
{

	PhyloPop_params p;
	p.window_size = 500;
	p.alfreq_scaling = 4;

	CountData counts("/Users/pickrell/projects/treemix/dog/data/box_bull_bm_mass.gz", &p);
	pair<double, double> f4 = counts.calculate_f4();
	cout << f4.first << " "<< f4.second << "\n";
	return 0;
}

