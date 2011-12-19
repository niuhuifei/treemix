/*
 * SNP.cpp
 *
 *  Created on: Dec 12, 2011
 *      Author: pickrell
 */
#include "SNP.h"


SNP::SNP(vector<double> f, vector<string> n, F2_matrix* fm, map<string, double> * t){
	f2 = fm;
	names = n;
	freqs = f;
	lambda = 1;
	phi = (1+sqrt(5))/2;
	resphi = 2-phi;
	maxit = 20;
	trim = t;
}

double SNP::ss(){
	double toreturn = 0;
	for (int i = 0; i < names.size(); i++){
		for (int j = i+1; j < names.size(); j++){
			double tmp1 = freqs[i] - freqs[j];
			double t1 = trim->find(names[i])->second;
			double t2 = trim->find(names[j])->second;
			tmp1 = (tmp1* tmp1)- t1 - t2;
			double tmp2 = lambda*(f2->get(names[i], names[j]));
			double tmp3 = tmp1-tmp2;
			//cout << names[i] << " "<< names[j] << " "<< freqs[i] << " "<< freqs[j] <<" "<< f2->get(names[i], names[j]) << " "<< tmp1 << " "<< tmp2 << " "<< tmp3 << "\n";
			toreturn+= tmp3*tmp3;
		}
	}
	return toreturn;
}

int SNP::golden_section_lambda(double min, double guess, double max, double tau, int* nit){
	double x;

	if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
	else x = guess - resphi *(guess-min);
	if (fabs(max-min) < tau * (fabs(guess)+fabs(max)) || *nit > maxit) {
		double new_logweight = (min+max)/2;
		lambda = exp(new_logweight);
		return 0;
	}
	*nit = *nit+1;
	lambda = exp(x);

	double f_x = ss();

	lambda = exp(guess);

	double f_guess = ss();


	if (f_x < f_guess){
		if ( (max-guess) > (guess-min) )	return golden_section_lambda(guess, x, max, tau, nit);

		else return golden_section_lambda(min, x, guess, tau, nit);
	}
	else{
		if ( (max - guess) > (guess - min)  ) return golden_section_lambda(min, guess, x, tau, nit);
		else return golden_section_lambda(x, guess, max, tau, nit);
	}
}


void SNP::optimize_lambda(){

	double min, max, guess;
	guess = 1;
	guess = log(guess);
	min = -10;
	max =10;
	int nit = 0;
	golden_section_lambda(min, guess, max, 0.001, &nit);

}
