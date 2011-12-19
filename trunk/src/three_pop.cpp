/*
 * three_pop.cpp
 *
 *  Created on: Nov 9, 2011
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
	map<string, map<string, map<string, double> > > f3 = counts.calculate_f3s();
	for(map<string, map<string, map<string, double> > >::iterator it = f3.begin(); it != f3.end(); it++){
		for(map<string, map<string, double> >::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++){
			for (map<string, double>::iterator it3 = it2->second.begin(); it3 != it2->second.end(); it3++){
				out << it->first << " "<< it2->first << " "<< it3->first << " "<< it3->second << "\n";
			}
		}
	}
	return 0;
}


