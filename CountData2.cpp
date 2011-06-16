/*
 * CountData2.cpp
 *
 *  Created on: Jun 16, 2011
 *      Author: pickrell
 */

#include "CountData2.h"

CountData2::CountData2(string infile, string mpop, gsl_rng *r){
	read_counts(infile);
	meanpop = mpop;
	cout << "npop:"<< npop<< " nsnp:"<<nsnp<< "\n";
	alfreqs = gsl_matrix_alloc(nsnp, npop);
	scatter = gsl_matrix_alloc(npop-1, npop-1);
	cov = gsl_matrix_alloc(npop-1, npop-1);
	means = gsl_vector_alloc(nsnp);

	set_alfreqs();
	scale_alfreqs();
	center_alfreqs();
	randomize_order(r);
	set_scatter();
	set_cov();

}

void CountData2::read_counts(string infile){
	    allele_counts.clear();
	    pop2id.clear();
	    npop = 0;
	    nsnp = 0;
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
	            std::cerr<< "ERROR: cannot open file " << in << "\n";
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
	    	pop2id.insert(make_pair(line[i], i));
	    	npop ++;
	    }

	    /*
	     * read counts, store in allele_counts
	     */
	    while(getline(in, st)){
	            buf.clear();
	            stringstream ss(st);
	            line.clear();
	            while (ss>> buf){
	                    line.push_back(buf);
	            }
	            vector<pair<int, int> > topush;

	            for ( vector<string>::iterator it = line.begin(); it != line.end(); it++){
	                typedef boost::tokenizer<boost::char_separator<char> >
	                tokenizer;
	                boost::char_separator<char> sep(",");
	                tokenizer tokens(*it, sep);
	                vector<int> tmpcounts;
	                for (tokenizer::iterator tok_iter = tokens.begin();  tok_iter != tokens.end(); ++tok_iter){
	                        int tmp = atoi(tok_iter->c_str());
	                        tmpcounts.push_back(tmp);
	                }
	                if (tmpcounts.size() != 2){
	                	std::cerr << *it << " does not have two alleles\n";
	                	exit(1);
	                }
	                topush.push_back(make_pair(tmpcounts[0], tmpcounts[1]));
	            }
	            allele_counts.push_back(topush);
	            nsnp++;
	    }

}


void CountData2::set_alfreqs(){
	for (int i = 0; i < nsnp; i++){
		for (int j = 0; j < npop; j++){
			int c1 = allele_counts[i][j].first;
			int c2 = allele_counts[i][j].second;
			double f = (double) c1 / ( (double) c1 + (double) c2 );
			gsl_matrix_set(alfreqs, i, j, f);
		}
	}
}


void CountData2::scale_alfreqs(){
	for (int i = 0; i < nsnp; i++){
		for (int j = 0; j < npop; j++){
			double f = gsl_matrix_get(alfreqs, i, j);
			double scaled = asin(sqrt(f));
			gsl_matrix_set(alfreqs, i, j, scaled);
		}
	}
}


void CountData2::center_alfreqs(){
	map<string, int>::iterator centerpop = pop2id.find(meanpop);
	if (centerpop == pop2id.end()){
		cerr << "No such population:" << meanpop << "\n";
		exit(1);
	}
	int which = centerpop->second;
	for (int i = 0; i < nsnp; i++){
		float m = gsl_matrix_get(alfreqs, i, which);
		gsl_vector_set(means, i, m);
		for (int j = 0; j < npop; j++){
			float f = gsl_matrix_get(alfreqs, i, j);
			gsl_matrix_set(alfreqs, i, j, f - m);
		}
	}
}

void CountData2::randomize_order(gsl_rng *r){

	// randomize the order in which the (non-mean) populations appear in the scatter and covariance matrices

	popsnomean.clear();
	vector<string> tmppop;
	for (map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++){
		string tmp = it->first;
		if (tmp == meanpop) continue;
		tmppop.push_back(tmp);
	}
	int src[npop-1];
	for (int i = 0; i < npop-1; i++) src[i] = i;
	gsl_ran_shuffle(r, src, npop-1, sizeof(int));
	for( int i = 0; i < npop-1; i++){
		int which = src[i];
		popsnomean.push_back(tmppop.at(which));
		//cout << which << " "<< tmppop.at(which) << "\n";
	}
}


void CountData2::set_scatter(){
	for (int i = 0; i < popsnomean.size(); i++){
		for (int j = i; j < popsnomean.size(); j++){
			double c = 0;
			int which1 = pop2id[popsnomean[i]];
			int which2 = pop2id[popsnomean[j]];
			for (int k = 0; k < nsnp; k++){
				double toadd = gsl_matrix_get(alfreqs, k, which1) * gsl_matrix_get(alfreqs, k, which2);
				c += toadd;
			}
			gsl_matrix_set(scatter, i, j, c);
			gsl_matrix_set(scatter, j, i, c);
		}
	}
}


void CountData2::process_scatter(){
	// get the log of the determinant of the covariance matrix
	// also calculate the log of the multivariate gamma function
	scatter_det = 0;
	scatter_gamma = 0;
	int n = nsnp-1;
	size_t pop = npop;
	int s;
	gsl_matrix * work = gsl_matrix_alloc(pop,pop);
	gsl_matrix_memcpy( work, scatter );
	gsl_permutation * p = gsl_permutation_alloc(pop);

	//do LU decomposition and get determinant
	gsl_linalg_LU_decomp( work, p, &s );
	scatter_det = gsl_linalg_LU_lndet( work );
	//get the log sum of the gammas
	scatter_gamma = ( (double) npop * ( (double)  npop-1.0) /4.0) * log (M_PI);
	for (int i = 1; i <= npop; i++) scatter_gamma+= gsl_sf_lngamma( (double) n/2.0 + (1.0- (double) i)/2.0);
}


void CountData2::set_cov(){
	for (int i = 0; i < npop-1; i++){
		for (int j = i; j < npop-1; j++){
			double sc = gsl_matrix_get(scatter, i, j);
			double c = sc/ ( (double) nsnp - 1);
			gsl_matrix_set(cov, i, j, c);
			gsl_matrix_set(cov, j, i, c);
		}
	}
}


void CountData2::print_scatter(string outfile){
	ogzstream out(outfile.c_str());
	for (int i = 0; i < popsnomean.size(); i++) out << popsnomean.at(i) << " ";
	out << "\n";
	for (int i = 0; i < popsnomean.size(); i++){
		out << popsnomean.at(i);
		for(int j = 0; j < popsnomean.size(); j++)	 out << " "<< gsl_matrix_get(scatter, i, j);
		out << "\n";
	}

}


void CountData2::print_cov(string outfile){
	ogzstream out(outfile.c_str());
	for (int i = 0; i < popsnomean.size(); i++) out << popsnomean.at(i) << " ";
	out << "\n";
	for (int i = 0; i < popsnomean.size(); i++){
		out << popsnomean.at(i);
		for(int j = 0; j < popsnomean.size(); j++)	 out << " "<< gsl_matrix_get(cov, i, j);
		out << "\n";
	}

}
