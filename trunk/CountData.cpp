/*
 * CountData.cpp
 *
 *  Created on: Apr 1, 2011
 *      Author: pickrell
 */
#include "CountData.h"



CountData::CountData(string infile, PhyloPop_params* p){
	params = p;
	read_counts(infile);
	if (p->restrict_pop) npop = p->pops2use;
	cout << "npop:"<< npop<< " nsnp:"<<nsnp<< "\n";
	alfreqs = gsl_matrix_alloc(nsnp, npop);
	scatter = gsl_matrix_alloc(npop, npop);
	cov = gsl_matrix_alloc(npop, npop);
	cov_var = gsl_matrix_alloc(npop, npop);
	cov_var2 = gsl_matrix_alloc(npop, npop);
	U = gsl_matrix_alloc(npop-1, npop);
	scatter_prime = gsl_matrix_alloc(npop-1, npop-1);
	nblock = nsnp/ p->window_size;
	ncomp = (npop* (npop-1))/2+ npop;
	//cout << nwind << " "<< ncomp << "\n";
	//cov_samp = gsl_matrix_alloc(nblock, ncomp);
	//gsl_matrix_set_zero(cov_samp);
	cov_cov = gsl_matrix_alloc(ncomp, ncomp);
	set_alfreqs();
	scale_alfreqs();
	set_scatter();
	//process_scatter();
	set_cov();
	///set_cov2();
	//set_ne();
	//set_ne2();
	//set_ncomp_ef();
	//cout << "Effective number of SNPs: "<< ne << "\n";
	//process_cov();
}


CountData::CountData(string infile){
	cerr << "Do not use this constructor for CountData!\n"; exit(1);
	read_counts(infile);
	cout << "npop:"<< npop<< " nsnp:"<<nsnp<< "\n";
	alfreqs = gsl_matrix_alloc(nsnp, npop);
	scatter = gsl_matrix_alloc(npop, npop);
	set_alfreqs();
	//scale_alfreqs(which);
	set_scatter();
	process_scatter();
}


CountData::CountData(CountData * c, vector<string> names, gsl_matrix* model, PhyloPop_params* p, gsl_rng *r){
	// set covariance matrix to a random Wishart with covariance from a model, copy over the standard errors
	params= p;
	npop = names.size();
	nsnp = c->nsnp;
	cov = gsl_matrix_alloc(npop, npop);
	cov_var = gsl_matrix_alloc(npop, npop);
	ne = c->ne;
	for (int i = 0; i < names.size(); i++){
		pop2id.insert(make_pair(names[i], i));
	}
	set_cov_ran(model, r);
	//for(int i = 0; i < names.size(); i++){
	//	for (int j = 0; j < names.size(); j++){
	//		gsl_matrix_set(cov_var, i, j, c->get_cov_var(names[i], names[j] ));
	//	}
	//}
}

void CountData::set_cov_ran(gsl_matrix* model,gsl_rng* r){

	// 1. Take SVD of model
	gsl_matrix * U2 = gsl_matrix_alloc(npop-1,npop);
	gsl_matrix * A = gsl_matrix_alloc(npop,npop);
	gsl_matrix * VT = gsl_matrix_alloc(npop,npop);
	gsl_matrix * model_prime = gsl_matrix_alloc(npop-1,npop-1);
	gsl_vector * S = gsl_vector_alloc(npop);
	gsl_vector * work = gsl_vector_alloc(npop);
	gsl_matrix_memcpy( A, model);
	gsl_matrix_set_zero(cov);
	gsl_matrix_set_zero(cov_var);
	gsl_linalg_SV_decomp(A, VT, S, work);


	// Now copy the first npop-1 eigenvectors to U

	for (int i = 0; i < npop-1; i++){
		for(int j = 0; j < npop; j++){
			gsl_matrix_set(U2, i, j, gsl_matrix_get(A, i, j));
		}
	}


	// multiply U model U^T
	gsl_matrix * US = gsl_matrix_alloc(npop-1, npop);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, U2, model, 0.0, US);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, US, U2, 0.0, model_prime);
	int nblock = nsnp/ params->window_size;

	//initialize blocks
	vector< vector< vector<double> > > allsamp;
	for (int i = 0; i < npop; i++){
		vector<vector<double> > tmp1;
		for(int j = 0; j < npop; j++){
			vector<double> tmp;
			tmp1.push_back(tmp);
		}
		allsamp.push_back(tmp1);
	}

	//sample nblock Wisharts
	for (int i = 0; i < nblock; i++){
		gsl_matrix * cov_prime = gsl_matrix_alloc(npop-1,npop-1);
		gsl_matrix * tmpcov = gsl_matrix_alloc(npop, npop);
		//get random wishart
		//const gsl_rng *r, const int n, const int dof, const gsl_matrix *scale, gsl_matrix *result){
		rwishart(r, npop-1, ne , model_prime, cov_prime );
		//multiply U^T cov_prime U
		gsl_matrix * UTC = gsl_matrix_alloc(npop, npop-1);
		gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, U2, cov_prime, 0.0, UTC);
		gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, UTC, U2, 0.0, tmpcov);

		//cout << gsl_matrix_get(tmpcov, 1, 3) / ne << "\n";
		for(int i = 0; i < npop; i++){
			for (int j = i; j < npop; j++){
				gsl_matrix_set(cov, i, j, gsl_matrix_get(cov, i, j)+gsl_matrix_get(tmpcov, i, j)/ (ne/nblock));
				gsl_matrix_set(cov, j, i, gsl_matrix_get(cov, j, i)+gsl_matrix_get(tmpcov, j, i)/ (ne/nblock));
				allsamp[i][j].push_back(gsl_matrix_get(tmpcov, i, j)/ (ne/nblock));
			}
		}

		gsl_matrix_free(cov_prime);
		gsl_matrix_free(tmpcov);
	}
	//cout << gsl_matrix_get(model, 1, 3) << " model\n";
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			double mean =  gsl_matrix_get(cov, i, j)/ (double) nblock;
			gsl_matrix_set(cov, i, j, mean);
			gsl_matrix_set(cov, j, i, mean);
			double sum = 0;
			vector<double> all_covs = allsamp[i][j];
			for (vector<double>::iterator it = all_covs.begin(); it != all_covs.end(); it++) sum+= (*it-mean)*(*it-mean);
			//double sd = sqrt(sum/ (double) nblock);
			double c = sqrt(sum) /(double) nblock;
			gsl_matrix_set(cov_var, i, j, c);
			gsl_matrix_set(cov_var, j, i, c);

		}
	}
	//print_cov("rancov.gz");
	//print_cov_var("rancovvar.gz");
	gsl_matrix_free(U2);
	gsl_matrix_free(A);
	gsl_matrix_free(VT);
	gsl_vector_free(S);
	gsl_vector_free(work);
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


vector<string> CountData::list_pops(){
	vector<string> toreturn;
	for (map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++) {
		if ( params->restrict_pop == true ){
			if (it->second < params->pops2use ) toreturn.push_back( it->first);
		}
		else toreturn.push_back( it->first);
	}
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



void CountData::read_alfreqs(string infile){
	npop =0;
	gsl_matrix_free(alfreqs);
	gsl_matrix_free(scatter);
	pop2id.clear();
	vector<vector<double> > tmpalfreqs;
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
    while(getline(in, st)){
             buf.clear();
             stringstream ss(st);
             line.clear();
             while (ss>> buf){
                     line.push_back(buf);
             }
             vector<double> tmp;
             for (vector<string>::iterator it = line.begin(); it != line.end(); it++) tmp.push_back(atof(it->c_str()));
             tmpalfreqs.push_back(tmp);
    }
    nsnp = tmpalfreqs.size();
    npop = tmpalfreqs[0].size();
    scatter = gsl_matrix_alloc(npop, npop);
    alfreqs = gsl_matrix_alloc(nsnp, npop);
    for (int i = 0; i < nsnp; i++){
    	for (int j = 0; j< npop; j++) gsl_matrix_set(alfreqs, i, j, tmpalfreqs[i][j]);
    }
    for(int i = 0; i < npop; i++){
    	stringstream ss;
    	ss << "pop";
    	ss << i;
    	string name = ss.str();
    	pop2id.insert(make_pair(name, i));
    }
    scale_alfreqs();
    set_scatter();

    //print_scatter("testout_scatter.gz");
	//cout << "here\n";
    set_cov();
    process_cov();

}

void CountData::read_counts(string infile){
    allele_counts.clear();
    pop2id.clear();
    id2pop.clear();
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
    	pop2id.insert(make_pair(line[i], i));
    	id2pop.insert(make_pair(i, line[i]));
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
	mean_ninds.clear();
	mean_hzy.clear();
	id2nsnp.clear();

	for (int i = 0; i < npop; i++){
		mean_ninds.insert(make_pair(i, 0.0));
		mean_hzy.insert(make_pair(i, 0.0));
		id2nsnp.insert(make_pair(i, 0));
	}
	for (int i = 0; i < nsnp; i++){
		for (int j = 0; j < npop; j++){
			int c1 = allele_counts[i][j].first;
			int c2 = allele_counts[i][j].second;
			double f = (double) c1 / ( (double) c1 + (double) c2 );
			if ( c1+c2 < 1){
				cerr << "Warning: no counts at SNP "<< i << " population "<< j <<"\n";
				continue;
			}
			gsl_matrix_set(alfreqs, i, j, f);
			mean_ninds[j] += ((double) c1+ (double) c2)/2.0;
			mean_hzy[j] += 2*f*(1-f);
			id2nsnp[j]++;
		}
	}
	for (int i = 0; i < npop; i++){
		mean_ninds[i] = mean_ninds[i]/ id2nsnp[i];
		mean_hzy[i] = mean_hzy[i]/ id2nsnp[i];
	}
}

void CountData::scale_alfreqs(){
	for (int i = 0; i < nsnp; i++){
		double total = 0;
		for (int j = 0; j < npop; j++){
			double f = gsl_matrix_get(alfreqs, i, j);
			double scaled;
			if (params->alfreq_scaling == 1) scaled = asin(sqrt(f));
			else scaled = f;
			total = total+scaled;
			gsl_matrix_set(alfreqs, i, j, scaled);
		}

		double m = total/ (double) npop;
		for (int j = 0; j < npop; j++){
			double f = gsl_matrix_get(alfreqs, i, j);
			if (params->alfreq_scaling == 3) {
				gsl_matrix_set(alfreqs, i, j, (f-m)/sqrt(m *(1-m)) );
				if (m < 1e-8) gsl_matrix_set(alfreqs, i, j, 0);
			}
			else if (params->alfreq_scaling == 4) gsl_matrix_set(alfreqs, i, j, f);
			else gsl_matrix_set(alfreqs, i, j, f-m);
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

void CountData::set_cov(){
	/*
	 * Calculate covariance matrix in blocks on SNPs
	 *  cov[i,j] = mean( cov[i,j]_k ) over all k blocks
	 */
	gsl_matrix_free(cov);
	gsl_matrix_free(cov_var);
	//cout << npop << "\n";
	cov = gsl_matrix_alloc(npop, npop);
	cov_var = gsl_matrix_alloc(npop, npop);
	//cov_var2 = gsl_matrix_alloc(npop, npop);
	//initialize block estimation of covariance matrix
	vector<vector<vector<double> > > cov_block;
	for (int i = 0; i < npop; i++){
		vector<vector<double> > tmp1;
		for(int j = 0; j < npop; j++){
			vector<double> tmp;
			tmp1.push_back(tmp);
		}
		cov_block.push_back(tmp1);
	}
	vector<string> popnames = list_pops();
	cov_samp.clear();
	for (int i = 0; i < popnames.size(); i++){
		map<string, vector<double> > tmp;
		for (int j = 0; j < popnames.size(); j++){
			vector<double> tmp2;
			tmp.insert(make_pair(popnames[j], tmp2));
		}
		cov_samp.insert(make_pair(popnames[i], tmp));
	}
	//if trimming covariances, get the amount to trim
	map<string,double> trim;
	double sumtrim = 0;
	for ( map<string, int>::iterator it = pop2id.begin(); it!= pop2id.end(); it++){
		int id = it->second;
		string pop = it->first;
		double meanhzy = mean_hzy.find(id)->second;
		double mean_n = mean_ninds.find(id)->second;
		double t = meanhzy / (4.0* mean_n);
		sumtrim+= t;
		//cout << pop  << " "<< t << " "<< meanhzy << " "<< mean_n << "\n";
		trim.insert(make_pair(pop, t));
	}
	//calculate the covariance matrix in each block
	cout << "Estimating covariance matrix in "<< nblock << " blocks of size "<< params->window_size <<"\n"; cout.flush();
	for (int k = 0; k < nblock ; k++){
		int index = 0;
		for(int i = 0; i < npop; i++){
			for (int j = i; j < npop; j++){
				double c = 0;
				for (int n = k*params->window_size; n < (k+1)*params->window_size; n++){
					if (isnan(gsl_matrix_get(alfreqs, n, i))) continue;
					double toadd = gsl_matrix_get(alfreqs, n, i) * gsl_matrix_get(alfreqs, n, j);
					c+= toadd;
				}
				double cov = c/ (double) params->window_size;
				//cout << k << " "<< index << " "<< cov << "\n";
				string p1 = id2pop[i];
				string p2 = id2pop[j];
				if (params->sample_size_correct){
					double bias1 = trim[p1];
					double bias2 = trim[p2];
					cov = cov + bias1/ (double) npop + bias2/ (double) npop - sumtrim/( (double) npop* (double) npop);
					if (i ==j) cov = cov - trim[p1];
				}
				cov_samp[p1][p2].push_back(cov);
				if (p1 != p2) cov_samp[p2][p1].push_back(cov);
				//gsl_matrix_set(cov_samp, k, index, cov);
				//cout << k << " "<< index << " "<< gsl_matrix_get(cov_samp, k, index) << "\n";
				cov_block[i][j].push_back(cov);
				index++;
			}
		}
	}
	//ofstream tout("test");
	//calculate the mean, standard error of covariance estimates
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			vector<double> all_covs = cov_block[i][j];
			double sum = 0;
			for (vector<double>::iterator it = all_covs.begin(); it != all_covs.end(); it++) sum+= *it;
			double mean = sum/nblock;
			gsl_matrix_set(cov, i, j, mean);
			gsl_matrix_set(cov, j, i, mean);

			// and standard error
			sum = 0;
			for (vector<double>::iterator it = all_covs.begin(); it != all_covs.end(); it++) sum+= (*it-mean)*(*it-mean);
			//double sd = sqrt(sum/ (double) nblock);
			double c = sqrt(sum) /(double) nblock;
			//cout << i << " "<< j << " "<< sd << " "<< c << "\n";
			gsl_matrix_set(cov_var, i, j, c);
			gsl_matrix_set(cov_var, j, i, c);
		}
	}
	//get the covariance in the estimates of the covariance matrix
	gsl_matrix_set_zero(cov_cov);
	/*for(int i = 0; i < ncomp; i++){
		double meani = 0;
		for (int k = 0; k < nblock; k++) meani+= gsl_matrix_get(cov_samp, k, i);
		meani = meani/ (double) nblock;

		for (int j = i ; j < ncomp; j++){
			double meanj = 0;
			for (int k = 0; k < nblock; k++) meanj+= gsl_matrix_get(cov_samp, k, j);
			meanj = meanj/ (double) nblock;
			double sum = 0;
			for (int k = 0; k < nblock; k++)		sum += (gsl_matrix_get(cov_samp, k, i)- meani )*(gsl_matrix_get(cov_samp, k, j) - meanj);
			sum = sum / (double) nblock;
			//sum = sqrt(sum)/ (double) nblock;
			gsl_matrix_set(cov_cov, i, j, sum );
			gsl_matrix_set(cov_cov, j, i, sum );
		}
	}
	*/
}


void CountData::set_cov_f2(){
	/*
	 * Calculate matrix of f_2 statistics in blocks
	 *
	 */
	gsl_matrix_free(cov);
	gsl_matrix_free(cov_var);
	//cout << npop << "\n";
	cov = gsl_matrix_alloc(npop, npop);
	cov_var = gsl_matrix_alloc(npop, npop);
	//cov_var2 = gsl_matrix_alloc(npop, npop);
	//initialize block estimation of covariance matrix
	vector<vector<vector<double> > > cov_block;
	for (int i = 0; i < npop; i++){
		vector<vector<double> > tmp1;
		for(int j = 0; j < npop; j++){
			vector<double> tmp;
			tmp1.push_back(tmp);
		}
		cov_block.push_back(tmp1);
	}
	vector<string> popnames = list_pops();
	cov_samp.clear();
	for (int i = 0; i < popnames.size(); i++){
		map<string, vector<double> > tmp;
		for (int j = 0; j < popnames.size(); j++){
			vector<double> tmp2;
			tmp.insert(make_pair(popnames[j], tmp2));
		}
		cov_samp.insert(make_pair(popnames[i], tmp));
	}
	//if trimming covariances, get the amount to trim
	map<string,double> trim;
	double sumtrim = 0;
	for ( map<string, int>::iterator it = pop2id.begin(); it!= pop2id.end(); it++){
		int id = it->second;
		string pop = it->first;
		double meanhzy = mean_hzy.find(id)->second;
		double mean_n = mean_ninds.find(id)->second;
		double t = meanhzy / (4.0* mean_n);
		sumtrim+= t;
		//cout << pop  << " "<< t << " "<< meanhzy << " "<< mean_n << "\n";
		trim.insert(make_pair(pop, t));
	}
	//calculate the covariance matrix in each block
	cout << "Estimating f_2 matrix in "<< nblock << " blocks of size "<< params->window_size <<"\n"; cout.flush();
	for (int k = 0; k < nblock ; k++){
		int index = 0;
		for(int i = 0; i < npop; i++){
			for (int j = i; j < npop; j++){
				double c = 0;
				for (int n = k*params->window_size; n < (k+1)*params->window_size; n++){
					if (isnan(gsl_matrix_get(alfreqs, n, i))) continue;
					double toadd = (gsl_matrix_get(alfreqs, n, i) - gsl_matrix_get(alfreqs, n, j));
					toadd = toadd*toadd;
					c+= toadd;
				}
				double cov = c/ (double) params->window_size;
				//cout << k << " "<< index << " "<< cov << "\n";
				string p1 = id2pop[i];
				string p2 = id2pop[j];
				if (params->sample_size_correct){
					double bias1 = trim[p1];
					double bias2 = trim[p2];
					cov = cov - bias1 -bias2;
				}
				cov_samp[p1][p2].push_back(cov);
				if (p1 != p2) cov_samp[p2][p1].push_back(cov);
				//gsl_matrix_set(cov_samp, k, index, cov);
				//cout << k << " "<< index << " "<< gsl_matrix_get(cov_samp, k, index) << "\n";
				cov_block[i][j].push_back(cov);
				index++;
			}
		}
	}
	//ofstream tout("test");
	//calculate the mean, standard error of covariance estimates
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			vector<double> all_covs = cov_block[i][j];
			double sum = 0;
			for (vector<double>::iterator it = all_covs.begin(); it != all_covs.end(); it++) sum+= *it;
			double mean = sum/nblock;
			gsl_matrix_set(cov, i, j, mean);
			gsl_matrix_set(cov, j, i, mean);

			// and standard error
			sum = 0;
			for (vector<double>::iterator it = all_covs.begin(); it != all_covs.end(); it++) sum+= (*it-mean)*(*it-mean);
			//double sd = sqrt(sum/ (double) nblock);
			double c = sqrt(sum) /(double) nblock;
			//cout << i << " "<< j << " "<< sd << " "<< c << "\n";
			gsl_matrix_set(cov_var, i, j, c);
			gsl_matrix_set(cov_var, j, i, c);
		}
	}
	//get the covariance in the estimates of the covariance matrix
	gsl_matrix_set_zero(cov_cov);
}


void CountData::set_cov2(){
	/*
	 * Get SE of covariance matrix without blocks
	 */
	gsl_matrix_free(cov_var2);
	cov_var2 = gsl_matrix_alloc(npop, npop);
	gsl_matrix *tmpcov = gsl_matrix_alloc(npop, npop);
	//calculate the covariance matrix

	for(int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			double c = 0;
			for (int k = 0; k < nsnp; k++){
				if (isnan(gsl_matrix_get(alfreqs, k, i))) continue;
				double toadd = gsl_matrix_get(alfreqs, k, i) * gsl_matrix_get(alfreqs, k, j);
				c+= toadd;
			}
			double cov = c/(double) nsnp;
			//cout << i << " "<< j << " "<< cov << "\n";
			gsl_matrix_set(tmpcov, i, j, cov);
			gsl_matrix_set(tmpcov, j, i, cov);
		}
	}

	//calculate the SE
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			double mean = gsl_matrix_get(tmpcov, i, j);
			double sum = 0;
			for (int k = 0; k < nsnp; k++){
				double s = gsl_matrix_get(alfreqs, k, i) * gsl_matrix_get(alfreqs, k, j);

				sum+= (s-mean)*(s-mean);
			}

			//double sd = sqrt(sum/ (double) nsnp);
			double c = sqrt(sum)/ (double) nsnp;
			//if (i == 0 && j ==2){
			//	cout << "sum "<< sum << " "<< nsnp << " "<< c << "\n";
			//}
			gsl_matrix_set(cov_var2, i, j, c);
			gsl_matrix_set(cov_var2, j, i, c);
		}
	}
	gsl_matrix_free(tmpcov);
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

void CountData::print_cov_cov(string outfile){
	ogzstream out(outfile.c_str());
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			out << id2pop[i]<< "."<< id2pop[j] << " ";
		}
	}
	out <<  "\n";
	int index1 = 0;
	for (int l = 0; l < npop; l++){
		for (int m = l; m< npop; m++){
			out << id2pop[l]<< "."<< id2pop[m] << " ";
			int index =0;
			for (int i = 0; i < npop; i++){
				for (int j = i; j < npop; j++){
					//cout << index1 << " " << index << "\n";
					out << gsl_matrix_get(cov_cov, index1, index)<< " ";
					index++;
				}
			}
			index1++;
			out << "\n";
		}
	}
}



void CountData::print_cov_samp(string outfile){
	ogzstream out(outfile.c_str());
	vector<string> popnames = list_pops();
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			out << popnames[i]<< "."<< popnames[j] << " ";
		}
	}
	out <<  "\n";
	for (int k = 0; k < nblock; k++){
		for (int i = 0; i < npop; i++){
			for (int j = i; j < npop; j++){
				string p1 = popnames[i];
				string p2 = popnames[j];
				out << cov_samp[p1][p2].at(k) << " ";
			}
		}

		out << "\n";
	}

}


void CountData::print_alfreqs(string outfile){
	ogzstream out(outfile.c_str());
	for(map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++)	out << it->first << " ";
	out << "\n";
	for (int i = 0; i < nsnp ; i++){
		for(map<string, int>::iterator it2 = pop2id.begin(); it2 != pop2id.end(); it2++)	 out << " "<< gsl_matrix_get(alfreqs, i, it2->second);
		out << "\n";
	}
}

void CountData::print_fst(string outfile){
	gsl_matrix* fst = gsl_matrix_alloc(npop, npop);
	ogzstream out(outfile.c_str());
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			double total = 0;
			for (int k = 0; k < nsnp; k++){
				double f1 = gsl_matrix_get(alfreqs, k, i);
				double f2 = gsl_matrix_get(alfreqs, k, j);
				double f_hat = (f1+f2)/2;
				double between = 2*f_hat*(1-f_hat);
				double within = f1*(1-f1)+ f2*(1-f2);
				double fst_ijk;
				if (between < 1e-8) fst_ijk = 0;
				else fst_ijk = (between-within)/between;
				//cout << fst_ijk
				total+= fst_ijk;
			}
			total = total / (double) nsnp;
			gsl_matrix_set(fst, i, j, total);
			gsl_matrix_set(fst, j, i, total);
		}
	}
	for(map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++)	out << it->first << " ";
	out << "\n";
	for(map<string, int>::iterator it = pop2id.begin(); it != pop2id.end(); it++){
		out << it->first;
		for(map<string, int>::iterator it2 = pop2id.begin(); it2 != pop2id.end(); it2++)	 out << " "<< gsl_matrix_get(fst, it->second, it2->second);
		out << "\n";
	}
}

void CountData::process_scatter(){
	// project scatter matrix into (npop-1) space
	// get matrix for transformation
	// get determinant and gamma

	scatter_det = 0;
	scatter_gamma = 0;
	int n = nsnp-1;
	size_t pop = npop;
	int s;

	//first do SVD on scatter matrix
	gsl_matrix * A = gsl_matrix_alloc(pop,pop);
	gsl_matrix * VT = gsl_matrix_alloc(pop,pop);
	gsl_vector * S = gsl_vector_alloc(pop);
	gsl_vector * work = gsl_vector_alloc(pop);
	gsl_matrix_memcpy( A, scatter );

	gsl_linalg_SV_decomp(A, VT, S, work);

	// Now copy the first npop=1 eigenvectors to U

	for (int i = 0; i < npop-1; i++){
		for(int j = 0; j < npop; j++){
			gsl_matrix_set(U, i, j, gsl_matrix_get(A, i, j));
		}
	}

	// And transform scatter into m-1 space
	// S' = U S U^T

	gsl_matrix * US = gsl_matrix_alloc(pop-1, pop);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, U, scatter, 0.0, US);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, US, U, 0.0, scatter_prime);

	// Now do LU decomposition and get determinant of scatter_prime
	gsl_matrix_free(A);
	A = gsl_matrix_alloc(npop-1, npop-1);
	gsl_matrix_memcpy( A, scatter_prime );
	gsl_permutation * p = gsl_permutation_alloc(pop-1);
	gsl_linalg_LU_decomp( A, p, &s );
	scatter_det = gsl_linalg_LU_lndet( A );


	//get the log sum of the gammas
	//cout << npop << " "<< n << "\n";
	scatter_gamma = ( (double) (pop-1) * ( (double)  npop-2.0) /4.0) * log (M_PI);
	//cout << scatter_gamma << " sg1\n";
	for (int i = 1; i <= pop-1; i++) scatter_gamma+= gsl_sf_lngamma( (double) n/2.0 + (1.0- (double) i)/2.0);
	//cout << scatter_gamma << " sg2\n";
	cout << "scatter_gamma "<< scatter_gamma << "\n";
	cout << "scatter_det "<< scatter_det << "\n";
	gsl_matrix_free(A);
	gsl_matrix_free(VT);
	gsl_matrix_free(US);
	gsl_vector_free(S);
	gsl_vector_free(work);
	gsl_permutation_free(p);

}

void CountData::process_cov(){
	cov_var = gsl_matrix_alloc(npop, npop);
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			double c = 0;

			double tmp_cov = gsl_matrix_get(cov, i, j);
			for (int k = 0; k < nsnp; k++){
				double toadd = gsl_matrix_get(alfreqs, k, i) * gsl_matrix_get(alfreqs, k, j);
				c += (toadd-tmp_cov) * (toadd-tmp_cov);
			}
			gsl_matrix_set(cov_var, i, j, c/nsnp);
			gsl_matrix_set(cov_var, j, i, c/nsnp);
		}
	}
	//double tmp = 0;
	//int todivide =0;
	//for (int i = 0; i < npop; i++){
	//	for (int j = i; j < npop; j++)	{
	//		//cout << i << " "<< j << " "<< gsl_matrix_get(var_matrix, i, j)<< "\n";
	//		tmp+= gsl_matrix_get(var_matrix, i, j);
	//		todivide++;
	//	}
	//}
	//tmp = tmp/ (double) todivide;
	//cout << "tmp "<< tmp << "\n";
	//cov_var = tmp;
}

double CountData::get_cov(string pop1, string pop2){

	if (pop2id.find(pop1)  == pop2id.end()) {
		cerr << "No population "<< pop1 << "\n";
		exit(1);
	}
	if (pop2id.find(pop2)  == pop2id.end()) {
		cerr << "No population "<< pop2 << "\n";
		exit(1);
	}


	int p1 = pop2id[pop1];
	int p2 = pop2id[pop2];
	double toreturn = gsl_matrix_get(cov, p1, p2);
	return toreturn;
}


double CountData::get_scatter(string pop1, string pop2){

	if (pop2id.find(pop1)  == pop2id.end()) {
		cerr << "No population "<< pop1 << "\n";
		exit(1);
	}
	if (pop2id.find(pop2)  == pop2id.end()) {
		cerr << "No population "<< pop2 << "\n";
		exit(1);
	}


	int p1 = pop2id[pop1];
	int p2 = pop2id[pop2];
	double toreturn = gsl_matrix_get(scatter, p1, p2);
	return toreturn;
}

double CountData::get_cov_var(string pop1, string pop2){
	if (pop2id.find(pop1)  == pop2id.end()) {
		cerr << "No population "<< pop1 << "\n";
		exit(1);
	}
	if (pop2id.find(pop2)  == pop2id.end()) {
		cerr << "No population "<< pop2 << "\n";
		exit(1);
	}


	int p1 = pop2id[pop1];
	int p2 = pop2id[pop2];
	double toreturn = gsl_matrix_get(cov_var, p1, p2);
	return toreturn;
}

void CountData::print_cov(string outfile){
	ogzstream out(outfile.c_str());
	vector<string> pops = list_pops();
	for (int i = 0; i < pops.size(); i++) out << pops.at(i)<< " ";
	out << "\n";
	for (int i = 0; i < pops.size(); i++){
		out << pops.at(i);
		for(int j = 0; j < pops.size(); j++)	 out << " "<< gsl_matrix_get(cov, pop2id[pops[i]], pop2id[pops[j]]);
		out << "\n";
	}

}


void CountData::print_cov_var(string outfile){
	ogzstream out(outfile.c_str());
	vector<string> pops = list_pops();
	for (int i = 0; i < pops.size(); i++) out << pops.at(i)<< " ";
	out << "\n";
	for (int i = 0; i < pops.size(); i++){
		out << pops.at(i);
		for(int j = 0; j < pops.size(); j++)	 out << " "<< gsl_matrix_get(cov_var, pop2id[pops[i]], pop2id[pops[j]]);
		out << "\n";
	}

}


void CountData::print_cov_var2(string outfile){
	ogzstream out(outfile.c_str());
	vector<string> pops = list_pops();
	for (int i = 0; i < pops.size(); i++) out << pops.at(i)<< " ";
	out << "\n";
	for (int i = 0; i < pops.size(); i++){
		out << pops.at(i);
		for(int j = 0; j < pops.size(); j++)	 out << " "<< gsl_matrix_get(cov_var2, pop2id[pops[i]], pop2id[pops[j]]);
		out << "\n";
	}

}


void CountData::set_ne(){
	ne = 0;
	double tmpne = 0;
	int pairs = 0;
	for (int i = 0; i < npop; i++){
		for(int j = i ; j < npop; j++){
			double tmp =gsl_matrix_get(cov_var2, i, j)/gsl_matrix_get(cov_var, i, j);
			tmp = tmp*tmp;
			//if (i == 0&& j == 2){
			//	cout << "tmp "<< tmp <<"\n";
			//}
			tmpne+= tmp;
			pairs++;
		}
	}
	//cout << tmpne << " "<< pairs << "\n";
	tmpne = tmpne/ (double) pairs;
	tmpne = tmpne* nsnp;
	ne = int(tmpne);
}


void CountData::set_ne2(){
	ne2 = 0;
	int p = npop;
	gsl_matrix * A = gsl_matrix_alloc(npop, npop);
	gsl_matrix * VT = gsl_matrix_alloc(p,p);
	gsl_vector * S = gsl_vector_alloc(p);
	gsl_vector * work = gsl_vector_alloc(p);
	gsl_matrix_memcpy( A, scatter );

	gsl_linalg_SV_decomp(A, VT, S, work);

	double s  = 0;
	double s2 = 0;
	//for (int i = 0; i < p; i++){
	//	double eig = gsl_vector_get(S, i);
	//	cout << eig << "\n";
	//}
	for (int i = 0; i < p-1; i++){
		double eig = gsl_vector_get(S, i);
		eig = eig*eig;
		s+= eig;
		s2+= eig*eig;
		//cout << eig << " "<< eig*eig << "\n";
	}
	//cout << s<< " "<< s2 << "\n";
	double num = ( (double) p+1.0)* s*s;
	double denom = ( ((double) p-1.0)* s2 ) - (s*s);
	//cout << "num denom " << num << " "<< denom << "\n";
	double tmp = num/denom;
	ne2 = int(tmp);
	gsl_matrix_free(A);
	gsl_matrix_free(VT);
	gsl_vector_free(S);
	gsl_vector_free(work);

}

int rwishart(gsl_rng *r, const int n, const int dof, const gsl_matrix *scale, gsl_matrix *result){
  /* Wishart distribution random number generator */
  /*
   *    n        gives the dimension of the random matrix
   *    dof      degrees of freedom
   *    scale    scale matrix of dimension n x n
   *    result   output variable with a single random matrix Wishart distribution generation
   */
  int k,l;
  gsl_matrix *work = gsl_matrix_calloc(n,n);

  for(k=0; k<n; k++){
    gsl_matrix_set( work, k, k, sqrt( gsl_ran_chisq( r, (dof-k) ) ) );
    for(l=0; l<k; l++){
      gsl_matrix_set( work, k, l, gsl_ran_ugaussian(r) );
    }
  }
  gsl_matrix_memcpy(result,scale);
  gsl_linalg_cholesky_decomp(result);
  gsl_blas_dtrmm(CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,1.0,result,work);
  gsl_blas_dsyrk(CblasUpper,CblasNoTrans,1.0,work,0.0,result);
  gsl_matrix_free(work);
  return 0;
}

string CountData::get_pop_in_index(int index){
	string toreturn;
	map<string , int>::iterator it = pop2id.begin();
	while (it != pop2id.end()){
		if (it->second == index) return it->first;
		it++;
	}
	if (it == pop2id.end()) {
		cerr << "Trying to get index "<< index << " in CountData, none found\n";
		exit(1);
	}
	return toreturn;
}

/*
void CountData::set_ncomp_ef(){
	ncomp_ef = 0;
	gsl_matrix * S = gsl_matrix_alloc(ncomp, ncomp);
	gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, cov_samp, cov_samp, 0.0, S);

	gsl_matrix * A = gsl_matrix_alloc(ncomp,ncomp);
	gsl_matrix * VT = gsl_matrix_alloc(ncomp,ncomp);
	gsl_vector * Sv = gsl_vector_alloc(ncomp);
	gsl_vector * work = gsl_vector_alloc(ncomp);
	gsl_matrix_memcpy( A, S );
	gsl_linalg_SV_decomp(A, VT, Sv, work);
	//for (int i = 0; i < ncomp ; i++){
	//	cout << gsl_vector_get(Sv, i) << "\n";
	//}
	while (ncomp_ef < ncomp && (gsl_vector_get(Sv, ncomp_ef)  > 1e-10)) ncomp_ef++; //this is the number of eigenvectors to use
	//cout << "Effective number of comparisons "<< ncomp_ef << "\n";
	gsl_matrix_free(S);
	gsl_matrix_free(A);
	gsl_matrix_free(VT);
	gsl_vector_free(Sv);
	gsl_vector_free(work);
}
*/

pair<double, double> CountData::calculate_f4(){
	double mean, se;
	vector<double> f4_block;

	//calculate the covariance matrix in each block
	cout << "Estimating f_4 in "<< nblock << " blocks of size "<< params->window_size <<"\n"; cout.flush();
	for (int k = 0; k < nblock ; k++){
		double c = 0;
		for (int n = k*params->window_size; n < (k+1)*params->window_size; n++){
			if (isnan(gsl_matrix_get(alfreqs, n, 0))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, 1))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, 2))) continue;
			if (isnan(gsl_matrix_get(alfreqs, n, 3))) continue;
			double toadd = (gsl_matrix_get(alfreqs, n, 1) - gsl_matrix_get(alfreqs, n, 0))*(gsl_matrix_get(alfreqs, n, 3) - gsl_matrix_get(alfreqs, n, 2) );
			//cout << gsl_matrix_get(alfreqs, n, 0) << " "<< gsl_matrix_get(alfreqs, n, 1) << " "<< gsl_matrix_get(alfreqs, n, 2)<< " "<< gsl_matrix_get(alfreqs, n, 3)<< " "<< toadd << "\n";
			c+= toadd;
		}
		double cov = c/ (double) params->window_size;
		f4_block.push_back(cov);

	}
	//calculate the mean, standard error of covariance estimates
	double sum = 0;
	for (vector<double>::iterator it = f4_block.begin(); it != f4_block.end(); it++) sum+= *it;
	mean = sum/nblock;

	// and standard error
	sum = 0;
	for (vector<double>::iterator it = f4_block.begin(); it != f4_block.end(); it++) sum+= (*it-mean)*(*it-mean);
	//double sd = sqrt(sum/ (double) nblock);
	se = sqrt(sum) /(double) nblock;

	return make_pair(mean, se);
}

map<string, map<string, map<string, double> > > CountData::calculate_f3s(){
	map<string, map<string, map<string, double> > > toreturn;
	for(int i = 0; i< npop; i++){
		string p1 = id2pop[i];
		map<string, map<string, double> > tmp;
		for(int j = 0; j< npop ; j++){
			string p2 = id2pop[j];

			map<string, double> tmp2;
			for (int k = 0; k < npop; k++){

				string p3 = id2pop[k];
				if (p1 == p2 || p1 == p3 || p2 == p3) continue;
				double sum = 0;
				for (int l = 0; l < nsnp; l++){
					double n = mean_ninds[i];
					double f1 = gsl_matrix_get(alfreqs, l,i);
					double f2 = gsl_matrix_get(alfreqs, l, j);
					double f3 = gsl_matrix_get(alfreqs, l, k);
					double hz1 = f1*(1-f1);
					double s = (f1-f2)*(f1-f3) - hz1/ (2*n);
					sum += s;
				}
				sum = sum/nsnp;
				tmp2.insert(make_pair(p3, sum));
			}
			tmp.insert(make_pair(p2, tmp2));
		}
		toreturn.insert(make_pair(p1, tmp));
	}
	return toreturn;
}

double CountData::calculate_f2(int p1, int p2){
	double toreturn = 0;
	for (int i = 0; i < nsnp; i++){
		double n1 = mean_ninds[p1];
		double n2 = mean_ninds[p2];
		double f1 = gsl_matrix_get(alfreqs, i, p1);
		double f2 = gsl_matrix_get(alfreqs, i, p2);

		double diff = f1-f2;
		double hz1 = f1*(1-f1);
		double hz2 = f2*(1-f2);
		double toadd = diff*diff - hz1/(2*n1) - hz2/ (2*n2);
		toreturn += toadd;
	}
	toreturn = toreturn/ (double) nsnp;
	return toreturn;
}
void CountData::set_cov_jackknife(int which){
	gsl_matrix_set_zero(cov);
	for(int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			string p1 = id2pop[i];
			string p2 = id2pop[j];
			double m = 0;
			for (int k = 0; k < nblock; k++){
				if (k == which) continue;
				m+= cov_samp[p1][p2].at(k);
			}
			m = m/ (double) (nblock-1);
			gsl_matrix_set(cov, i, j, m);
			gsl_matrix_set(cov, j, i, m);
		}
	}

}


void CountData::set_cov_bootstrap(gsl_rng *r){
gsl_matrix_set_zero(cov);
	for(int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			string p1 = id2pop[i];
			string p2 = id2pop[j];
			double m = 0;
			for (int k = 0; k < nblock; k++){
				int rint = gsl_rng_uniform_int(r, nblock);
				m+= cov_samp[p1][p2].at(rint);
			}
			m = m / (double) nblock;
			gsl_matrix_set(cov, i, j, m);
			gsl_matrix_set(cov, j, i, m);
		}
	}
}


void CountData::set_cov_fromsamp(int which){
	gsl_matrix_set_zero(cov);
	for(int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			string p1 = id2pop[i];
			string p2 = id2pop[j];
			double m = cov_samp[p1][p2].at(which);
			gsl_matrix_set(cov, i, j, m);
			gsl_matrix_set(cov, j, i, m);
		}
	}

}
