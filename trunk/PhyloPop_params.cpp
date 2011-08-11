/*
 * PhyloPop_params.cpp
 *
 *  Created on: Jul 1, 2011
 *      Author: pickrell
 */

#include "PhyloPop_params.h"

PhyloPop_params::PhyloPop_params(){
	bias_correct = true;
	window_size = 1;
	alfreq_scaling = 0;
	global = false;
	readtree = false;
	treefile = "NA";
	tau = 0.5;
	minweight = -10;
	maxweight = 10;
	nmig = 0;
}
