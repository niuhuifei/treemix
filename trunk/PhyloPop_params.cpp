
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
	alfreq_scaling = 4;
	global = false;
	readtree = false;
	treefile = "NA";
	tau = 0.1;
	minweight = -10;
	maxweight = 10;
	nmig = 0;
	m_neigh = 5;
	maxit = 10;
	set_root = false;
	root = "NA";
	read_graph = false;
	vfile = "NA";
	efile = "NA";
	quick = false;
	min_migw = 0.001;
	nofrac = false;
	smooth_lik = true;
	smooth_scale = 1;
	nrand = 0;
	restrict_pop = false;
	pops2use = 0;
	sample_size_correct = true;
	calc_se = false;
	f2 = true;
	neg_penalty = 10;
}
