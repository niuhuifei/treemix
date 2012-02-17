/*
 * SNP.h
 *
 *  Created on: Dec 12, 2011
 *      Author: pickrell
 */

#ifndef SNP_H_
#define SNP_H_
#include "F2_matrix.h"
#include "Settings.hpp"

class SNP{
public:
	SNP(vector<double>, vector<string>, F2_matrix *, map<string, double> *);
	vector<double> freqs;
	vector<string> names;
	double lambda;
	F2_matrix * f2;
	double ss();
	int optimize_lambda();
	int golden_section_lambda(double, double, double, double, int *);
	double phi, resphi;
	int maxit;
	map<string, double> * trim;
	double sumtrim;
	int npop;
};


#endif /* SNP_H_ */
