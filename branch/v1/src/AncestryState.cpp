/*
 * AncestryState.cpp
 *
 *  Created on: Nov 11, 2011
 *      Author: pickrell
 */

#include "AncestryState.h"

AncestryState::AncestryState(){

}
AncestryState::AncestryState(CountData* co, CountData* co2,  PhyloPop_params* p ){
	counts = co;
	counts2 = co2;
	params= p;
	c = 0.0264966;
	sigma = gsl_matrix_alloc(3, 3);
	sigma_cor = gsl_matrix_alloc(3, 3);
	gsl_matrix_set_zero(sigma);
	gsl_matrix_set_zero(sigma_cor);
}

void AncestryState::set_sigma(int which){
	double indhzy = counts->mean_hzy[0];
	double indn = counts->mean_ninds[0];
	double c3 = indhzy/(2*indn);
	gsl_matrix_set(sigma, 0, 0, c+ c3);
	gsl_matrix_set(sigma, 1, 1, c);
	gsl_matrix_set(sigma, 2, 2, c);
	gsl_matrix_set(sigma, 2, 1, 0);
	gsl_matrix_set(sigma, 1, 2, 0);
	if (which  == 0){
		gsl_matrix_set(sigma, 0,1, c);
		gsl_matrix_set(sigma, 1, 0, c);
		gsl_matrix_set(sigma, 0, 2, 0);
		gsl_matrix_set(sigma, 2, 0, 0);
	}
	else if (which == 1){
		gsl_matrix_set(sigma, 0,1, 0);
		gsl_matrix_set(sigma, 1, 0, 0);
		gsl_matrix_set(sigma, 0, 2, c);
		gsl_matrix_set(sigma, 2, 0, c);
	}
	else if (which == 2){
		gsl_matrix_set(sigma, 0,1, c/2);
		gsl_matrix_set(sigma, 1, 0, c/2);
		gsl_matrix_set(sigma, 0, 2, c/2);
		gsl_matrix_set(sigma, 2, 0, c/2);
	}
}



void AncestryState::set_sigmacor_from_sigma(){

	gsl_matrix_set_zero(sigma_cor);
	double c1 = 1.0 / (double) 3.0;
	double c2 = c1*c1;
	double shared = 0;
	for(int i = 0; i < 3; i++){
		for (int j = 0; j < 3; j++){
			shared += gsl_matrix_get(sigma, i, j);
		}
	}

	shared = shared * c2;
	//cout << "here "<< shared << " "<< current_npops << " \n";
	for (int i = 0; i < 3; i++){
		for (int j = i; j < 3; j++){
			//cout << i << " "<< j;
			double vij = gsl_matrix_get(sigma, i, j);
			double sum_i = 0;
			double sum_j = 0;
			for (int k = 0; k < 3; k++){
				sum_i += gsl_matrix_get(sigma, k, i);
				sum_j += gsl_matrix_get(sigma, j, k);
			}
			sum_i = sum_i * c1;
			sum_j = sum_j * c1;
			double w_ij = vij - sum_i - sum_j + shared;
			//cout << " "<< vij << " "<< sum_i << " "<< sum_j << " "<< w_ij << "\n";
			gsl_matrix_set(sigma_cor, i, j, w_ij);
			gsl_matrix_set(sigma_cor, j, i, w_ij);
		}
	}
}



double AncestryState::llik(){
	double toreturn = 0;
	for (int i = 0; i < 3; i++){
		for (int j = i; j < 3; j++){
			double pred = gsl_matrix_get(sigma_cor, i, j);
			double obs = gsl_matrix_get(counts->cov, i, j);
			double se = gsl_matrix_get(counts2->cov_var, i, j);

			double dif = obs-pred;
			double scale = 1;
			//cout << pred << " "<< obs << " "<< se << "\n";
			if (params->smooth_lik) scale = params->smooth_scale;
			double toadd = gsl_ran_gaussian_pdf(dif, se * scale);
			//cout << pred << " "<< obs << " "<< se << " "<< toadd << " "<< log(toadd) << "\n";
			//cout << p1<< " "<< p2 << " "<< toadd << " "<< log(toadd) << "\n";
			toreturn+= log(toadd);

			//}
		}
	}
	//cout <<  "\n";
	//cout << toreturn << "\n";
	//cout << "tmp "<< (double) tmp / (double) total <<"\n";
	return toreturn;
}

void AncestryState::print_ancestry_llik(string outfile){

	ofstream out(outfile.c_str());
	///for (int i = 0; i < 1; i++){
	for (int i = 0; i < counts->nblock; i++){
		counts->set_cov_fromsamp(i);
		out << i << " ";
		for (int j = 0; j < 3; j++){
			set_sigma(j);
			set_sigmacor_from_sigma();
			out << llik() << " ";
			cout << llik() << " ";
		}
		out << "\n";
		cout << "\n";
	}
}
