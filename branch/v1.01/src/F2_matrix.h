/*
 * F2_matrix.h
 *
 *  Created on: Dec 9, 2011
 *      Author: pickrell
 */

#ifndef F2_MATRIX_H_
#define F2_MATRIX_H_
#include "Settings.hpp"

class F2_matrix{
public:
	F2_matrix(string);
	void set_matrix(string);
	gsl_matrix * f2;
	int npop;
	map<string, int> name2index;
	map<int, string> index2name;
	void print();
	double get(string, string);
};


#endif /* F2_MATRIX_H_ */
