/*
 * State.cpp
 *
 *  Created on: Mar 31, 2011
 *      Author: pickrell
 */

#include "State.h"

State::State(string newick, CountData* counts){
	countdata = counts;
	tree = new PhyloPop_Tree::BinaryTree<PhyloPop_Tree::NodeData>(newick, countdata->pop2id);
	thetas = gsl_matrix_alloc(counts->nsnp, counts->npop);
	means = gsl_matrix_alloc(counts->nsnp, counts->npop);
	sigma = gsl_matrix_alloc(counts->npop, counts->npop);
	gsl_matrix_set_zero(thetas);
	gsl_matrix_set_zero(means);
	gsl_matrix_set_zero(sigma);
}

State::State(const State& oldstate){
	tree = oldstate.tree->copy();
	countdata = oldstate.countdata;
	thetas = gsl_matrix_alloc(oldstate.countdata->nsnp, oldstate.countdata->npop);
	gsl_matrix_memcpy(thetas, oldstate.thetas);
	means = gsl_matrix_alloc(oldstate.countdata->nsnp, oldstate.countdata->npop);
	gsl_matrix_memcpy(means, oldstate.means);
	sigma = gsl_matrix_alloc(oldstate.countdata->npop, oldstate.countdata->npop);
	gsl_matrix_memcpy(sigma, oldstate.sigma);
}

void State::compute_sigma(){
	map<int, PhyloPop_Tree::iterator<PhyloPop_Tree::NodeData> > id2node = tree->get_tips(tree->getRoot());
	for( int i = 0; i < countdata->npop; i++){
		for (int j = i; j < countdata->npop; j++){
			cout <<i << " "<< j  <<" "; cout.flush();
			if (i == j){
				double dist = tree->get_dist_to_root(id2node[i]);
				cout << dist << "\n"; cout.flush();
				gsl_matrix_set(sigma, i, j, dist);
				cout << "here \n"; cout.flush();
			}
			else{
				PhyloPop_Tree::iterator<PhyloPop_Tree::NodeData> lca = tree->get_LCA(tree->getRoot(), id2node[i], id2node[j]);
				double dist = tree->get_dist_to_root(lca);
				cout << dist << "\n";
				gsl_matrix_set(sigma, i, j, dist);
				gsl_matrix_set(sigma, j, i, dist);
				cout << "here \n"; cout.flush();
			}
		}
	}
}

void State::print_sigma(){
	for(int i = 0; i < countdata->npop; i++){
		for (int j = 0; j < countdata->npop; j++){
			cout << gsl_matrix_get(sigma, i, j) << " ";
		}
		cout << "\n";
	}


}

double State::llik(){
	double toreturn = 0;
	for (int i = 0; i < countdata->nsnp; i++){
		size_t r = i;
		gsl_vector* m = gsl_vector_alloc(countdata->npop);
		gsl_vector* theta_snp = gsl_vector_alloc(countdata->npop);
		gsl_matrix_get_row(m, means, r);
		gsl_matrix_get_row(theta_snp, thetas, r);
		double tmp = dmvnorm(countdata->npop, theta_snp, m, sigma);
		cout << tmp << "\n";
		toreturn+= log(tmp);
	}
	return toreturn;
}
