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
	snp_liks = gsl_vector_alloc(counts->nsnp);
	winv = gsl_matrix_alloc(counts->npop, counts->npop);
	ax =0;
	gsl_matrix_set_zero(thetas);
	gsl_vector_set_zero(means);
	gsl_matrix_set_zero(sigma);
	gsl_matrix_set_zero(winv);
	gsl_vector_set_zero(snp_liks);
	//traversal = tree.get_inorder_traversal(countdata->npop);
}


State::State(const State& oldstate){
	tree = oldstate.tree->copy();
	countdata = oldstate.countdata;
	params = oldstate.params;
	winv = gsl_matrix_alloc(oldstate.countdata->npop, oldstate.countdata->npop);
	gsl_matrix_memcpy(winv, oldstate.winv);
	ax =oldstate.ax;
	//traversal = oldstate.traversal;
	thetas = gsl_matrix_alloc(oldstate.countdata->nsnp, oldstate.countdata->npop);
	gsl_matrix_memcpy(thetas, oldstate.thetas);
	means = gsl_vector_alloc(oldstate.countdata->nsnp);
	gsl_vector_memcpy(means, oldstate.means);
	sigma = gsl_matrix_alloc(oldstate.countdata->npop, oldstate.countdata->npop);
	gsl_matrix_memcpy(sigma, oldstate.sigma);
	snp_liks = gsl_vector_alloc(oldstate.countdata->nsnp);
	gsl_vector_memcpy(snp_liks, oldstate.snp_liks);
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
	set_sigma_inv();

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
	toreturn = log(dens_mvnorm(theta_snp, m));
	gsl_vector_free(m);
	gsl_vector_free(theta_snp);
	return toreturn;
}

void State::update(gsl_rng* r){
	update_means(r);
	update_tree(r);
}

void State::update_tree(gsl_rng* r){
	double oldlik = current_lik;

	//flip the nodes
	tree->flip_sons(tree->getRoot(), r);

	//copy the old tree
	PhyloPop_Tree::Tree<PhyloPop_Tree::NodeData>* oldtree = tree->copy();

	//propose new tree
	vector<PhyloPop_Tree::iterator<PhyloPop_Tree::NodeData> > trav = propose_tree(r);

	//check to make sure the maximum distance to the root is not greater than B
	double maxdist = 0;
	for(int i = 0; i < trav.size(); i++){
		if (trav[i]->m_time > maxdist) maxdist = trav[i]->m_time;
	}

	//if it's ok, do a metropolis update
	if (maxdist < params->B){
		double newlik = llik();
		double ratio = exp(newlik-oldlik);
		if (ratio < 1){
			double acc = gsl_rng_uniform(r);
			if (acc > ratio){
				delete tree;
				tree = oldtree;
				compute_sigma();
			}
			else{
				delete oldtree;
				current_lik = newlik;
			}
		}
		else{
			delete oldtree;
			current_lik = newlik;
		}
	}

	//if not, switch back to the old tre
	else{
		delete tree;
		tree = oldtree;
		compute_sigma();
	}
}

void State::init_liks(){
	double total;
	for(int i = 0; i < countdata->nsnp; i++){
		double tmp = llik_snp(i);
		gsl_vector_set(snp_liks, i, tmp);
		total+= tmp;
	}
	current_lik = total;
}

vector<PhyloPop_Tree::iterator<PhyloPop_Tree::NodeData> > State::propose_tree(gsl_rng* r){
	//1. flip nodes
	//cout << "flipping\n"; cout.flush();
	//tree->flip_sons(tree->getRoot(), r);
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
	double total= 0;
	for(int i = 0; i < countdata->nsnp; i++){
		total+= update_mean(r, i);
	}
	current_lik = total;
}

double State::update_mean(gsl_rng* r, int i){
	/*
	 * Metropolis update of ancestral allele frequency parameters. Prior is Beta(lambda, lambda) (symmetric)
	 *  proposal is oldparameter+ N(0, s^2), constrained to fall within [0,1]
	 */
	double oldm = gsl_vector_get(means, i);
	double newm = oldm + gsl_ran_gaussian(r, params->s2);
	double toreturn = gsl_vector_get(snp_liks, i);
	if (newm > 0 && newm < 1){
		double oldpost = toreturn * gsl_ran_beta_pdf(oldm, params->lambda, params->lambda);
		gsl_vector_set(means, i, newm);
		double newlik = llik_snp(i);
		double newpost = newlik *gsl_ran_beta_pdf(newm, params->lambda, params->lambda);
		double ratio = exp(newpost-oldpost);
		if (ratio < 1){
			double acc = gsl_rng_uniform(r);
			if (acc > ratio) gsl_vector_set(means, i, oldm);
			else {
				gsl_vector_set(snp_liks, i, newlik);
				toreturn = newlik;
			}
		}
		else{
			gsl_vector_set(snp_liks, i, newlik);
			toreturn = newlik;
		}
	}
	return toreturn;
}



void State::print_state(ogzstream& treefile, ogzstream& mfile){
	string t = tree->get_newick_format();
	//cout <<t << "\n";
	treefile << t << "\n";
	//cout << "printed tree\n"; cout.flush();
	for (int i = 0; i < countdata->nsnp; i++){
		mfile << gsl_vector_get(means, i) << " ";
	}
	mfile << "\n";
	//cout << llik() << "\n";

}

void State::read_thetas(string thetafile, int npop, int nsnp){
	countdata->npop = npop;
	countdata->nsnp = nsnp;
	gsl_vector_free(means);
	gsl_matrix_free(sigma);
	gsl_matrix_free(thetas);
	gsl_vector_free(snp_liks);
	thetas = gsl_matrix_alloc(countdata->nsnp, countdata->npop);
	means = gsl_vector_alloc(countdata->nsnp);
	snp_liks = gsl_vector_alloc(countdata->nsnp);
	sigma = gsl_matrix_alloc(countdata->npop, countdata->npop);
	gsl_matrix_set_zero(thetas);
	gsl_vector_set_zero(means);
	gsl_vector_set_zero(snp_liks);
	gsl_matrix_set_zero(sigma);

	ifstream in(thetafile.c_str());
	struct stat stFileInfo;
	int intStat;
	string st, buf;
    intStat = stat(thetafile.c_str(), &stFileInfo);
    if (intStat !=0){
    	std::cerr<< "ERROR: cannot open file " << in << "\n";
    	exit(1);
	}
    int i = 0;
    while(getline(in, st)){
			vector<string> line;
            string buf;
            stringstream ss(st);
            while (ss>> buf){
                    line.push_back(buf);
            }
            for(int j = 0; j < line.size(); j++)	gsl_matrix_set(thetas, i, j, atof(line[j].c_str()));
            i++;
    }
}

void State::set_sigma_inv(){
	ax = 0;
	size_t pop = countdata->npop;
	int s;
	gsl_matrix * work = gsl_matrix_alloc(pop,pop);
	gsl_matrix * newinv = gsl_matrix_alloc(pop,pop);
	gsl_matrix_memcpy( work, sigma );

	gsl_permutation * p = gsl_permutation_alloc(countdata->npop);
	gsl_linalg_LU_decomp( work, p, &s );
	gsl_linalg_LU_invert( work, p, newinv );
	ax = gsl_linalg_LU_det( work, s );
	gsl_matrix_free( work );
	gsl_permutation_free( p );
	gsl_matrix_memcpy( winv, newinv);
	gsl_matrix_free(newinv);
}

double State::dens_mvnorm(const gsl_vector* x, const gsl_vector* mean){
	// x are the thetas, mean are the ancestral allele frequencies. assume we've already done the inversion of sigma
	double ay;
	gsl_vector *ym, *xm;

	xm = gsl_vector_alloc(countdata->npop);
	gsl_vector_memcpy( xm, x);
	gsl_vector_sub( xm, mean );
	ym = gsl_vector_alloc(countdata->npop);


	gsl_blas_dsymv(CblasUpper,1.0,winv,xm,0.0,ym);
	gsl_blas_ddot( xm, ym, &ay);
	gsl_vector_free(xm);
	gsl_vector_free(ym);
	ay = exp(-0.5*ay)/sqrt( pow((2*M_PI),countdata->npop)*ax );
	return ay;
}
