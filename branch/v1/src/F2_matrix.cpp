/*
 * F2_matrix.cpp
 *
 *  Created on: Dec 9, 2011
 *      Author: pickrell
 */

#include "F2_matrix.h"

F2_matrix::F2_matrix(string infile){
	f2 = gsl_matrix_alloc(1, 1);
	set_matrix(infile);
}

void F2_matrix::set_matrix(string infile){
	name2index.clear();
	string ext = infile.substr(infile.size()-3, 3);
	if (ext != ".gz"){
		std::cerr << infile << " is not gzipped (only .gz files accepted)\n";
		exit(1);
	}
	igzstream in(infile.c_str()); //only gzipped files
	vector<string> line;
	struct stat stFileInfo;
	int intStat;
	string st, buf;

	intStat = stat(infile.c_str(), &stFileInfo);
	if (intStat !=0){
		std::cerr<< "ERROR: cannot open file " << infile << "\n";
		exit(1);
	   }

	/*
	 * header contains population names
	 */
	getline(in, st);
	stringstream ss(st);
	line.clear();
	while (ss>> buf){
		line.push_back(buf);
	}
	/*
	 * make map from header, number populations according to order
	 */
	for(int i = 0; i < line.size(); i++) {
		name2index.insert(make_pair(line[i], i));
		index2name.insert(make_pair(i, line[i]));
	}
	gsl_matrix_free(f2);
	npop = name2index.size();
	f2 = gsl_matrix_alloc(npop, npop);

	/*
	 * read counts, store in allele_counts
	 */
	int i =0;
	while(getline(in, st)){
		buf.clear();
		stringstream ss(st);
		line.clear();
		while (ss>> buf){
			line.push_back(buf);
		}
		for (int j = 1; j < line.size(); j++){
			float d = atof(line[j].c_str());
			gsl_matrix_set(f2, i, j-1, d);
		}
		i++;
	}
}

void F2_matrix::print(){
	for (int i = 0; i < npop; i++)cout << index2name[i] << " ";
	cout << "\n";
	for (int i = 0; i < npop; i++){
		cout << index2name[i] << " ";
		for (int j = 0; j < npop; j++) cout << gsl_matrix_get(f2, i, j)<< " ";
		cout <<  "\n";
	}
}

double F2_matrix::get(string p1, string p2){
	int i1 = name2index[p1];
	int i2 = name2index[p2];
	return gsl_matrix_get(f2, i1, i2);
}
