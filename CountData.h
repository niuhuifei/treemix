/*
 * CountData.h
 *
 *  Created on: Apr 1, 2011
 *      Author: pickrell
 */

#ifndef COUNTDATA_H_
#define COUNTDATA_H_
#include "Settings.hpp"


class CountData{
public:
	CountData(string);
	void read_counts(string);
	map<string, int> pop2id;
	vector<vector<pair<int, int> > > allele_counts;
	int npop, nsnp;
};

#endif /* COUNTDATA_H_ */
