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

	NNLS_SOLVER nnls_o(4, 10);
	nnls_o.maxIter() = 10;
	double* A = new double[4 * 2];
	//A[0] = 0.0372;
	//A[1] = 0.2869;
	//A[2] = 0.6816;
	//A[3] = 0.7071;
	//A[4] = 0.6233;
	//A[5] = 0.6245;
	//A[6] = 0.6344;
	//A[7] = 0.6170;

	A[0] = 0.0372;
	A[1] = 0.6861;
	A[2] = 0.6233;
	A[3] = 0.6344;
	A[4] = 0.2869;
	A[5] = 0.7071;
	A[6] = 0.6245;
	A[7] = 0.6170;

	//A[0] = 1;
	//A[1] = 2;
	//A[2] = 3;
	//A[3] = 4;
	//A[4] = 2;
	//A[5] = 2;
	//A[6] = 3;
	//A[7] = 4;
	double* b = new double[4];
	double* x = new double[2];
	double* w = new double[2];
	double* zz = new double[4];
	double rNorm;
	b[0] = 0.8587;
	b[1] = 0.1781;
	b[2] = 0.0747;
	b[3] = 0.8405;
	//b[0] = 2;
	//b[1] = 4;
	//b[2] = 6;
	//b[3] = 8;
	int maxIter = 10;
	int * index = new int[6];
	int mode;
	//nnls(A, 2, 4, 2, b, x, &rNorm, w, zz, index,  &mode, maxIter);
	bool converged = nnls_o.solve(A, 2, b, x, rNorm);
	//cout << converged << "\n";
	for (int i = 0; i < 2; i++) cout << x[i] << "\n";
	//cout << rNorm << "\n";
	//a = new double*[HEIGHT];
	//for (int i = 0; i < HEIGHT; ++i)	a[i] = new double[WIDTH];

	  	  // Assign values
	//a[0][0] = 3.6;
	//a[1][2] = 4.0;

	  // De-Allocate memory to prevent memory leak
	// for (int i = 0; i < HEIGHT; ++i) delete [] a[i];
	// delete [] a;
	/*
	PhyloPop_params p;
	p.snpinfo = true;
	CountData counts("testin.gz", &p);
	F2_matrix test("TreeMix.modelcov.gz");
	pair<vector<string>, vector<double> > tmpsnp = counts.get_freqs(3);

	map<string,double> trim;
	double sumtrim = 0;
	for ( map<string, int>::iterator it = counts.pop2id.begin(); it!= counts.pop2id.end(); it++){
		int id = it->second;
		string pop = it->first;
		double meanhzy = counts.mean_hzy.find(id)->second;
		double mean_n = counts.mean_ninds.find(id)->second;
		double t = meanhzy / (4.0* mean_n);
		sumtrim+= t;
		//cout << pop  << " "<< t << " "<< meanhzy << " "<< mean_n << "\n";
		trim.insert(make_pair(pop, t));
	}

	SNP testsnp(tmpsnp.second, tmpsnp.first, &test, &trim);
	testsnp.optimize_lambda();
	cout << "\n";
	cout << testsnp.lambda << " "<< testsnp.ss() << "\n";
*/
	return 0;

}
