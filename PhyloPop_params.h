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

	int nmig; //number of migration edges to add

};

#endif /* PHYLOPOP_PARAMS_H_ */
