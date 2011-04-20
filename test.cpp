/*
 * test.cpp
 *
 *  Created on: Mar 29, 2011
 *      Author: pickrell
 */

//#include "Tree.h"
//#include "BinaryTree.h"
#include "MCMC.h"
#include "CountData.h"

int main(){
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_ranlxs2;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, 100);
	string testnewick = "((YRI:0.8,LWD:0.4):0.1,(CEU:0.4,(JPT:0.2,CHB:0.3):0.2):0.5);";
	ogzstream treeout("treeout.gz");
	ogzstream mout("meanout.gz");
	CountData counts("testin_counts.gz");
	//vector<vector<double> > freqs = counts.get_alfreqs();
	for (int i = 0; i < counts.nsnp ; i++){
		for (int j = 0; j < counts.npop; j++){
			std::cerr << gsl_matrix_get(counts.alfreqs, i, j) << " ";
		}
		std::cerr << "\n";
	}

	counts.set_cov();
	cout <<  "\n";
	for (int i = 0; i < counts.npop ; i++){
		for (int j = 0; j < counts.npop; j++){
			cout << gsl_matrix_get(counts.cov, i, j) << " ";
		}
		cout << "\n";
	}
/*	MCMC_params p;
	State teststate(testnewick, &counts, &p);
	//teststate.read_thetas("testin_thetas", 5, 1000);
	teststate.tree->randomize_tree(r);
	teststate.compute_sigma();
	teststate.init_thetas();
	teststate.init_means();
	teststate.init_liks();
	cout << teststate.countdata->npop << " "<< teststate.countdata->nsnp << " "<< teststate.current_lik<< "\n";
	MCMC testmcmc(&teststate, &p);
	testmcmc.run(r, treeout, mout);

	vector<PhyloPop_Tree::iterator< PhyloPop_Tree::NodeData> > vec = teststate.tree->get_inorder_traversal(counts.npop);
	for (int i = 0; i < vec.size(); i++){
		cout << vec[i]->m_id << " "<< vec[i]->m_len << " "<< vec[i]->m_time<< "\n";
	}
	teststate.tree->set_node_heights(vec);
	cout << "\n";
	for (int i = 0; i < vec.size(); i++){
			cout << vec[i]->m_id << " "<< vec[i]->m_len << " "<< vec[i]->m_time<< "\n";
		}
	cout <<"\n";
	teststate.tree->randomize_tree(r);
	vec = teststate.tree->get_inorder_traversal(counts.npop);
	for (int i = 0; i < vec.size(); i++){
		cout << vec[i]->m_id << " "<< vec[i]->m_len << " "<< vec[i]->m_time<< "\n";
	}
	cout << teststate.tree->get_newick_format() << "\n";
*/
/*
	cout << "\n";
	teststate.compute_sigma();
	for (int i = 0 ; i < 10000; i++){
	teststate.update(r);
	tmp << gsl_vector_get(teststate.means, 0) << "\n";
	}
	cout << "\n";
	vec = teststate.tree->get_inorder_traversal(counts.npop);
	for (int i = 0; i < vec.size(); i++){
			cout << vec[i]->m_id << " "<< vec[i]->m_len << " "<< vec[i]->m_time<< "\n";
		}

	cout << "llik: "<< teststate.llik() << "\n";
*/
};