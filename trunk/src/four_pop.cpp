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

	CountData counts("/Users/pickrell/projects/southaf/data/test.gz", &p);
	set<pair<string, pair<double, double> > >  f4 = counts.calculate_f4();
	for (set<pair<string, pair<double, double> > >::iterator it = f4.begin(); it != f4.end(); it++){
		cout << it->first << " "<< it->second.first << " "<< it->second.second << " "<< it->second.first/ it->second.second  << "\n";
		//
	}
	return 0;
}

