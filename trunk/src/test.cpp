/*
 * test.cpp
 *
 *  Created on: Mar 29, 2011
 *      Author: pickrell
 */


#include "CountData.h"
#include "PopGraph.h"
#include "GraphState2.h"
#include "PhyloPop_params.h"
#include "SNP.h"
#include "nnls.h"
int main(){
	const gsl_rng_type * T;

	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_ranlxs2;
	r = gsl_rng_alloc(T);
	int seed = (int) time(0);
	gsl_rng_set(r, seed);

	double d = 5.0;
	double se = 1.0;
	//GraphState2 g();
	double t1 = gsl_ran_gaussian_pdf(d, se);
	double t2 = lndgauss(d, se);
	cout << t1 << " "<< t2 << "\n";
	return 0;

}
