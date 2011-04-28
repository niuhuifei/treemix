/*
 * CountData.cpp
 *
 *  Created on: Apr 1, 2011
 *      Author: pickrell
 */
#include "CountData.h"

CountData::CountData(string infile, int which){
	read_counts(infile);
	cout << "npop:"<< npop<< " nsnp:"<<nsnp<< "\n";
	alfreqs = gsl_matrix_alloc(nsnp, npop);
	scatter = gsl_matrix_alloc(npop, npop);
	set_alfreqs();
	scale_alfreqs(which);
	set_scatter();
	process_scatter();
}

string CountData::get_pops(){
	string toreturn = "(";
	map<string , int>::iterator it = pop2id.begin();
	map<string , int>::iterator it2 = pop2id.end();
	it2--;
	toreturn+= it->first +":0.1";
	it++;
	while (it != it2){
		toreturn+=",("+it->first+":0.1";
		it++;
	}
	toreturn+=","+it->first +":0.1";
	for (int i = 0; i < npop; i++)	toreturn+= "):0.1";


	toreturn = toreturn.substr(0, toreturn.size()-4);
	toreturn+= ";";
	return toreturn;
}

void CountData::read_scatter(string infile){
	npop =0;
	gsl_matrix_free(scatter);
	vector<vector<double> > tmpcov;
	ifstream in(infile.c_str()); //only gzipped files
    vector<string> line;
    struct stat stFileInfo;
    int intStat;
    string st, buf;

    intStat = stat(infile.c_str(), &stFileInfo);
    if (intStat !=0){
            std::cerr<< "ERROR: cannot open file " << in << "\n";
            exit(1);
    }
    while(getline(in, st)){
             buf.clear();
             stringstream ss(st);
             line.clear();
             while (ss>> buf){
                     line.push_back(buf);
             }
             vector<double> tmp;
             for (vector<string>::iterator it = line.begin(); it != line.end(); it++) tmp.push_back(atof(it->c_str()));
             tmpcov.push_back(tmp);
    }
    npop = tmpcov.size();
    scatter = gsl_matrix_alloc(npop, npop);
    for (int i = 0; i < npop; i++){
    	for (int j = 0; j< npop; j++) gsl_matrix_set(scatter, i, j, tmpcov[i][j]);
    }

}
void CountData::read_counts(string infile){
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

void CountData::set_alfreqs(){
	for (int i = 0; i < nsnp; i++){
		for (int j = 0; j < npop; j++){
			int c1 = allele_counts[i][j].first;
			int c2 = allele_counts[i][j].second;
			double f = (double) c1 / ( (double) c1 + (double) c2 );
			gsl_matrix_set(alfreqs, i, j, f);
		}
	}
}

void CountData::scale_alfreqs(int which){
	for (int i = 0; i < nsnp; i++){
		double total = 0;
		for (int j = 0; j < npop; j++){
			double f = gsl_matrix_get(alfreqs, i, j);
			double scaled;
			if (which ==1) scaled = asin(sqrt(f));
			else if (which ==2) scaled = f;
			total = total+scaled;
			gsl_matrix_set(alfreqs, i, j, scaled);
		}
		double m = total/ (double) npop;
		for (int j = 0; j < npop; j++){
			double f = gsl_matrix_get(alfreqs, i, j);
			if (which ==1 )gsl_matrix_set(alfreqs, i, j, f-m);
			else if (which == 2) gsl_matrix_set(alfreqs, i, j, (f-m)/ sqrt(m*(1-m)));
		}
	}
}

void CountData::set_scatter(){
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			double c = 0;
			for (int k = 0; k < nsnp; k++){
				double toadd = gsl_matrix_get(alfreqs, k, i) * gsl_matrix_get(alfreqs, k, j);
				c += toadd;
			}
			gsl_matrix_set(scatter, i, j, c);
			gsl_matrix_set(scatter, j, i, c);
		}
	}
}

void CountData::print_scatter(string outfile){
	ogzstream out(outfile.c_str());
	for(map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++)	out << it->first << " ";
	out << "\n";
	for(map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++){
		out << it->first;
		for(map<string, int>::iterator it2 = pop2id.begin(); it2 != pop2id.end(); it2++)	 out << " "<< gsl_matrix_get(scatter, it->second, it2->second);
		out << "\n";
	}

}

void CountData::process_scatter(){
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
	//cout << "scatter_gamma "<< scatter_gamma << "\n";
}
