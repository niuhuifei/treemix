/*
 * PhyloPop_params.h
 *
 *  Created on: Jul 1, 2011
 *      Author: pickrell
 */

#ifndef PHYLOPOP_PARAMS_H_
#define PHYLOPOP_PARAMS_H_

#include "Settings.hpp"

class PhyloPop_params{
public:
	PhyloPop_params();
	bool bias_correct, global, readtree;
	string treefile;
	int window_size;
	int alfreq_scaling; // 0 = no scaling, 1 = asin(sqrt(f))

	//optimization of weights
	double tau; //stopping criterion
	double minweight, maxweight;
	int maxit; //maximum number of iterations when maximizing weights

	int nmig; //number of migration edges to add

	int m_neigh; //"neighborhood" size for addition of migration events

	// setting root
	bool set_root;
	string root;

	// read from previous run
	bool read_graph;
	string vfile;
	string efile;

	//quick optimization of weights
	bool quick;

};

#endif /* PHYLOPOP_PARAMS_H_ */
