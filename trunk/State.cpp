/*
 * State.cpp
 *
 *  Created on: Mar 31, 2011
 *      Author: pickrell
 */

#include "State.h"

State::State(string newick, CountData* counts, MCMC_params* p){
	countdata = counts;
	params =p;
	tree = new PhyloPop_Tree::BinaryTree<PhyloPop_Tree::NodeData>(newick, countdata->pop2id);
	thetas = gsl_matrix_alloc(counts->nsnp, counts->npop);
	means = gsl_vector_alloc(counts->nsnp);
	sigma = gsl_matrix_alloc(counts->npop, counts->npop);
	gsl_matrix_set_zero(thetas);
	gsl_vector_set_zero(means);
	gsl_matrix_set_zero(sigma);
	//traversal = tree.get_inorder_traversal(countdata->npop);
}

State::State(const State& oldstate){
	tree = oldstate.tree->copy();
	countdata = oldstate.countdata;
	params = oldstate.params;
	//traversal = oldstate.traversal;
	thetas = gsl_matrix_alloc(oldstate.countdata->nsnp, oldstate.countdata->npop);
	gsl_matrix_memcpy(thetas, oldstate.thetas);
	means = gsl_vector_alloc(oldstate.countdata->nsnp);
	gsl_vector_memcpy(means, oldstate.means);
	sigma = gsl_matrix_alloc(oldstate.countdata->npop, oldstate.countdata->npop);
	gsl_matrix_memcpy(sigma, oldstate.sigma);
}

void State::compute_sigma(){
	map<int, PhyloPop_Tree::iterator<PhyloPop_Tree::NodeData> > id2node = tree->get_tips(tree->getRoot());
	for( int i = 0; i < countdata->npop; i++){
		for (int j = i; j < countdata->npop; j++){
			//cout <<i << " "<< j  <<" "; cout.flush();
			if (i == j){
				double dist = tree->get_dist_to_root(id2node[i]);
				//cout << dist << "\n"; cout.flush();
				gsl_matrix_set(sigma, i, j, dist);
				//cout << "here \n"; cout.flush();
			}
			else{
				PhyloPop_Tree::iterator<PhyloPop_Tree::NodeData> lca = tree->get_LCA(tree->getRoot(), id2node[i], id2node[j]);
				double dist = tree->get_dist_to_root(lca);
				//cout << dist << "\n";
				gsl_matrix_set(sigma, i, j, dist);
				gsl_matrix_set(sigma, j, i, dist);
				//cout << "here \n"; cout.flush();
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
		double tmp = llik_snp(i);
		//cout << tmp << "\n";
		toreturn+= tmp;
	}
	return toreturn;
}


double State::llik_snp(int i){

	/*
	 * the log-likelihood at a SNP is MVN(theta | m, sigma), where m is the ancestral allele frequency (shared by all the populations)
	 *  and sigma is the covariance matrix determined by the shape of the tree
	 */
	double toreturn = 0;
	size_t r = i;
	gsl_vector* m = gsl_vector_alloc(countdata->npop);
	gsl_vector* theta_snp = gsl_vector_alloc(countdata->npop);
	gsl_vector_set_all(m, gsl_vector_get(means, r));
	gsl_matrix_get_row(theta_snp, thetas, r);
	/*
	  for(int i = 0; i < countdata->npop; i++){
		cout << gsl_vector_get(m, i) << " "<< gsl_vector_get(theta_snp, i) <<"\n";
	}

	cout << "\n";
	for (int j = 0; j < countdata->npop; j++){
		for (int k = 0; k < countdata->npop; k++){
			cout << gsl_matrix_get(sigma, j, k)<< " ";
		}
		cout << "\n";
	}
	*/
	toreturn = log(dmvnorm(countdata->npop, theta_snp, m, sigma));
	return toreturn;
}

void State::update(gsl_rng* r){
	update_means(r);
	update_tree(r);
}

void State::update_tree(gsl_rng* r){
	//cout << "updating tree \n"; cout.flush();
	double oldlik = llik();
	PhyloPop_Tree::Tree<PhyloPop_Tree::NodeData>* oldtree = tree->copy();

	//propose new tree
	vector<PhyloPop_Tree::iterator<PhyloPop_Tree::NodeData> > trav = propose_tree(r);

	//check to make sure the maximum distance to the root is not greater than B
	double maxdist = 0;
	for(int i = 0; i < trav.size(); i++){
		if (trav[i]->m_time > maxdist) maxdist = trav[i]->m_time;
	}
	//cout << "max dist "<< maxdist << "\n";
	//if it's ok, do a metropolis update
	if (maxdist < params->B){
		double newlik = llik();
		double ratio = exp(newlik-oldlik);
		cout << oldlik  << " "<< newlik << " "<< ratio << "\n";
		if (ratio < 1){
			double acc = gsl_rng_uniform(r);
			if (acc > ratio){
				delete tree;
				tree = oldtree;
				compute_sigma();
			}
		}
	}
	else{
		delete tree;
		tree = oldtree;
		compute_sigma();
	}
	//cout << "updated \n"; cout.flush();

}

vector<PhyloPop_Tree::iterator<PhyloPop_Tree::NodeData> > State::propose_tree(gsl_rng* r){
	//1. flip nodes
	//cout << "flipping\n"; cout.flush();
	tree->flip_sons(tree->getRoot(), r);
	//2. get the traversal
	//cout << "get traversal 1\n"; cout.flush();
	vector<PhyloPop_Tree::iterator<PhyloPop_Tree::NodeData> > trav = tree->get_inorder_traversal(countdata->npop);
	//3. fiddle with read depths
	//cout << "perturb\n"; cout.flush();
	tree->perturb_node_heights(trav, params->epsilon, r);
	//4. rebuild the tree
	//cout << "rebild\n"; cout.flush();
	tree->build_tree(trav);
	//5. update the branch lengths
	//cout << "update branch lengths\n"; cout.flush();
	tree->update_branch_lengths(tree->getRoot());
	//6. get the new traversal
	//cout << "traverse 2\n"; cout.flush();
	trav = tree->get_inorder_traversal(countdata->npop);
	//7. reset the heights
	//cout << "reset heights\n"; cout.flush();
	tree->set_node_heights(trav);
	//8. recompute the correlation matrix
	//cout << "compute sigma\n"; cout.flush();
	compute_sigma();
	//cout << "done \n"; cout.flush();
	return trav;
}

void State::update_means(gsl_rng* r){
	//cout << "nsnps: "<< countdata->nsnp << "\n";
	for(int i = 0; i < countdata->nsnp; i++){
		update_mean(r, i);
	}
}

void State::update_mean(gsl_rng* r, int i){
	/*
	 * Metropolis update of ancestral allele frequency parameters. Prior is Beta(lambda, lambda) (symmetric)
	 *  proposal is oldparameter+ N(0, s^2), constrained to fall within [0,1]
	 */
	double oldm = gsl_vector_get(means, i);
	double newm = oldm + gsl_ran_gaussian(r, params->s2);
	if (newm > 0 && newm < 1){
		double oldpost = llik_snp(i) * gsl_ran_beta_pdf(oldm, params->lambda, params->lambda);
		gsl_vector_set(means, i, newm);
		double newpost = llik_snp(i) *gsl_ran_beta_pdf(newm, params->lambda, params->lambda);
		double ratio = exp(newpost-oldpost);
		if (ratio < 1){
			double acc = gsl_rng_uniform(r);
			if (acc > ratio) gsl_vector_set(means, i, oldm);
		}
		//cout << oldm << " "<< newm << " "<< oldpost << " "<< newpost <<  " "<< ratio<<"\n";
	}
}
