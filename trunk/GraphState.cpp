/*
 * GraphState.cpp
 *
 *  Created on: Apr 25, 2011
 *      Author: pickrell
 */


#include "GraphState.h"


GraphState::GraphState(string newick, CountData* counts, MCMC_params* p){
	countdata = counts;
	params =p;
	tree = new PopGraph(newick);
	tree_bk = new PopGraph(newick);
	sigma = gsl_matrix_alloc(counts->npop, counts->npop);
	gsl_matrix_set_zero(sigma);
	cout << tree->get_newick_format() <<"\n";
	//cout << "here\n"; cout.flush();
	compute_sigma();
	//cout <<  "here2\n"; cout.flush();

}



void GraphState::print_sigma(){
	for(int i = 0; i < countdata->npop; i++){
		for (int j = 0; j < countdata->npop; j++){
			cout << gsl_matrix_get(sigma, i, j) << " ";
		}
		cout << "\n";
	}
}

void GraphState::init(){
	current_lik = llik();
}


void GraphState::compute_sigma(){
	map<string, Graph::vertex_descriptor > popname2tip = tree->get_tips(tree->root);
	for( map<string, int>::iterator it1 = countdata->pop2id.begin(); it1 != countdata->pop2id.end(); it1++){
		for (map<string, int>::iterator it2 = countdata->pop2id.begin(); it2 != countdata->pop2id.end(); it2++){
			string p1 = it1->first;
			string p2 = it2->first;
			int i = it1->second;
			int j = it2->second;
			if (i == j){
				double dist = tree->get_dist_to_root(popname2tip[p1]);
				gsl_matrix_set(sigma, i, j, dist);
			}
			else{
				Graph::vertex_descriptor lca = tree->get_LCA(tree->root, popname2tip[p1], popname2tip[p2]);
				double dist = tree->get_dist_to_root(lca);
				gsl_matrix_set(sigma, i, j, dist);
				gsl_matrix_set(sigma, j, i, dist);
			}
		}
	}
}

void GraphState::read_sigma(string infile){
	gsl_matrix_free(sigma);
	vector<vector<double> > tmpcov;

	ifstream in(infile.c_str());
    vector<string> line;
    struct stat stFileInfo;
    int intStat;
    string st, buf;
    intStat = stat(infile.c_str(), &stFileInfo);
    if (intStat !=0){
            std::cerr<< "ERROR: cannot open file " << infile << "\n";
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
    int npop = countdata->npop;
    sigma = gsl_matrix_alloc(npop, npop);
    for (int i = 0; i < npop; i++){
    	for (int j = 0; j< npop; j++) gsl_matrix_set(sigma, i, j, tmpcov[i][j]);
    }

}
double GraphState::llik(){
	return dens_wishart();
}


void GraphState::print_state(ogzstream& treefile){
	string t = tree->get_newick_format();
	treefile << t << "\n";
}

void GraphState::local_update_tree(gsl_rng* r){
	local_update_tree_topology(r);
	local_update_tree_branches(r);
}
void GraphState::local_update_tree_topology(gsl_rng* r){
	//cout << "local update\n"; cout.flush();
	double oldlik = current_lik;

	//copy the old tree, get the traversal
	tree_bk->copy(tree);
	vector<Graph::vertex_descriptor> trav = tree->get_inorder_traversal(countdata->npop);
	//int nnode = 2*countdata->npop  - 1;

	//pick a random internal node (odd positions in inorder traversal)
	double ran = gsl_rng_uniform(r) * ((double) countdata->npop-1);
	int ranindex = int(2*ceil(ran) -1);

	//do the local update
	tree->local_update(trav[ranindex], r);
	trav = tree->get_inorder_traversal(countdata->npop);
	compute_sigma();

	//check to make sure the tree isn't too large
	double maxdist = 0;
	for(int i = 0; i < trav.size(); i++){
		if ( tree->get_height(trav[i]) > maxdist) maxdist = tree->get_height(trav[i]);
	}

	//if it's ok, do a metropolis update
	if (maxdist < params->B){
		//cout << "local update\n"; cout.flush();
		double newlik = llik();
		double ratio = exp(newlik-oldlik);
		//cout << newlik << " "<< oldlik << "\n";
		if (ratio < 1){
			double acc = gsl_rng_uniform(r);
			//cout << oldlik << " "<< newlik << " "<< ratio << "\n";
			if (acc > ratio){
				PopGraph * tmptree = tree;
				tree = tree_bk;
				tree_bk = tmptree;
				//tree->print();
				compute_sigma();
			}
			else{
				current_lik = newlik;
			}
		}
		else{
			current_lik = newlik;
		}
	}

	//if not, switch back to the old tree
	else{
		PopGraph * tmptree = tree;
		tree = tree_bk;
		tree_bk = tmptree;
		compute_sigma();
	}

}


void GraphState::local_update_tree_branches(gsl_rng* r){
	double oldlik = current_lik;

	//copy the old tree, get the traversal
	tree_bk->copy(tree);
	vector<Graph::vertex_descriptor> trav = tree->get_inorder_traversal(countdata->npop);
	//int nnode = 2*countdata->npop  - 1;

	//pick a random internal node (odd positions in inorder traversal)
	double ran = gsl_rng_uniform(r) * ((double) countdata->npop-1);
	int ranindex = int(2*ceil(ran) -1);

	//do the local update
	tree->local_update_branches(trav[ranindex], r, params->epsilon);
	trav = tree->get_inorder_traversal(countdata->npop);
	compute_sigma();

	//check to make sure the tree isn't too large
	double maxdist = 0;
	for(int i = 0; i < trav.size(); i++){
		if ( tree->get_height(trav[i]) > maxdist) maxdist = tree->get_height(trav[i]);
	}

	//if it's ok, do a metropolis update
	if (maxdist < params->B){
		double newlik = llik();
		double ratio = exp(newlik-oldlik);
		if (ratio < 1){
			double acc = gsl_rng_uniform(r);
			//cout << oldlik << " "<< newlik << " "<< ratio << "\n";
			if (acc > ratio){
				PopGraph * tmptree = tree;
				tree = tree_bk;
				tree_bk = tmptree;
				//tree->print();
				compute_sigma();
			}
			else{
				current_lik = newlik;
			}
		}
		else{
			current_lik = newlik;
		}
	}

	//if not, switch back to the old tree
	else{
		PopGraph * tmptree = tree;
		tree = tree_bk;
		tree_bk = tmptree;
		compute_sigma();
	}

}

void GraphState::update_tree(gsl_rng* r){
	double oldlik = current_lik;

	//flip the nodes
	tree->flip_sons(tree->root, r);

	//copy the old tree

	tree_bk->copy(tree);
	//cout <<"nope!\n";

	//propose new tree
	vector<Graph::vertex_descriptor> trav = propose_tree(r);
	//check to make sure the maximum distance to the root is not greater than B
	double maxdist = 0;
	for(int i = 0; i < trav.size(); i++){
		if ( tree->get_height(trav[i]) > maxdist) maxdist = tree->get_height(trav[i]);
	}

	//tree->print();
	//if it's ok, do a metropolis update
	if (maxdist < params->B){
		double newlik = llik();
		double ratio = exp(newlik-oldlik);
		if (ratio < 1){
			double acc = gsl_rng_uniform(r);
			if (acc > ratio){
				PopGraph * tmptree = tree;
				tree = tree_bk;
				tree_bk = tmptree;
				//tree->print();
				compute_sigma();
			}
			else{
				current_lik = newlik;
			}
		}
		else{
			current_lik = newlik;
		}
	}

	//if not, switch back to the old tree
	else{
		PopGraph * tmptree = tree;
		tree = tree_bk;
		tree_bk = tmptree;
		compute_sigma();
	}
}

vector<Graph::vertex_descriptor > GraphState::propose_tree(gsl_rng* r){

	//2. get the traversal
	//cout << "get traversal 1\n"; cout.flush();
	vector<Graph::vertex_descriptor> trav = tree->get_inorder_traversal(countdata->npop);
	//3. fiddle with read depths
	//cout << "perturb\n"; cout.flush();
	tree->perturb_node_heights(trav, params->epsilon, r);
	//4. rebuild the tree
	//cout << "rebild\n"; cout.flush();
	tree->build_tree(trav);
	//5. update the branch lengths
	//cout << "update branch lengths\n"; cout.flush();
	tree->update_branch_lengths(tree->root);
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


double GraphState::dens_wishart(){
	// density of the wishart distribution with covariance matrix sigma, n = number of snps-1, p = number of populations
	// 		scatter matrix has been stored in countdata->scatter, the ln(determinant) is in countdata->scatter_det
	// 		and the ln of the relevant multiariate gamma is in scatter_gamma
	//
	// density is ( [n-p-1]/2 * scatter_det - [1/2] trace [sigma^-1 * scatter] - [np/2] ln(2) - [n/2] ln(det(sigma)) - scatter_gamma

	double toreturn = 0;
	int s;
	double ld;
	int p = countdata->npop;
	int n = countdata->nsnp -1;

	//copy the covariance matrix over, declare inverse
	gsl_matrix * work = gsl_matrix_alloc(p,p);
	gsl_matrix_memcpy( work, sigma );
	gsl_permutation * perm = gsl_permutation_alloc(p);
	gsl_matrix * inv  = gsl_matrix_alloc(p, p);
	gsl_matrix * ViU = gsl_matrix_alloc(p, p);

	//do LU decomposition
	gsl_linalg_LU_decomp( work, perm, &s );

	//invert sigma
	gsl_linalg_LU_invert( work, perm, inv );

	//get log of determinant
	ld = gsl_linalg_LU_lndet( work );

	//multiply inverse of cov by scatter, get trace
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, inv, countdata->scatter, 0.0, ViU);
	double trace = 0;
	for (int i = 0; i < countdata->npop; i++) trace+= gsl_matrix_get(ViU, i, i);

	toreturn+= ( (double) n- (double) p-1.0 )/2.0 * countdata->scatter_det - trace/2.0;
	toreturn+= -( (double) n* (double) p/2.0) * log(2.0);
	toreturn += -((double) n/2.0)*ld;
	toreturn+= -countdata->scatter_gamma;

	gsl_matrix_free(work);
	gsl_matrix_free(inv);
	gsl_matrix_free(ViU);
	gsl_permutation_free(perm);

	return toreturn;
}


void GraphState::init_tree(gsl_rng* r){
	for (int i = 0; i < 50; i++){
		tree->flip_sons(tree->root, r);
		tree->randomize_tree(r);
	}

}

void GraphState::set_branches_ls(){
	vector<Graph::vertex_descriptor> i_nodes = tree->get_inorder_traversal_noroot(countdata->npop);

	//initialize the workspace
	int n = countdata->npop * countdata->npop;
	int p = 2*countdata->npop -2;
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, p);
	gsl_matrix * X = gsl_matrix_alloc(n, p);
	gsl_vector * y  = gsl_vector_alloc(n);
	gsl_vector * c = gsl_vector_alloc(p);
	gsl_matrix * cov = gsl_matrix_alloc(p, p);
	double chisq;

	//set up the workspace
	map<string, Graph::vertex_descriptor > popname2tip = tree->get_tips(tree->root);
	int index = 0;
	for( map<string, int>::iterator it1 = countdata->pop2id.begin(); it1 != countdata->pop2id.end(); it1++){
		for (map<string, int>::iterator it2 = countdata->pop2id.begin(); it2 != countdata->pop2id.end(); it2++){
			string p1 = it1->first;
			string p2 = it2->first;
			//cout << p1 << " "<< p2 << " here\n"; cout.flush();
			int i = it1->second;
			int j = it2->second;
			double empirical_cov = gsl_matrix_get(countdata->cov, i, j);
			set<Graph::vertex_descriptor> path;
			if (i == j) path = tree->get_path_to_root(popname2tip[p1]);

			else{
				Graph::vertex_descriptor lca = tree->get_LCA(tree->root, popname2tip[p1], popname2tip[p2]);
				path = tree->get_path_to_root(lca);
			}
			//cout << "here2\n"; cout.flush();
			gsl_vector_set(y, index, empirical_cov);
			for (int i2 =0 ; i2 < i_nodes.size(); i2++){
				Graph::vertex_descriptor tmp = i_nodes[i2];
				if  (path.find(tmp) == path.end() ) gsl_matrix_set(X, index, i2, 0);
				else  gsl_matrix_set(X, index, i2, 1);
			}
			index++;
		}
	}

	// fit the least squares estimates
	gsl_multifit_linear(X, y, c, cov, &chisq, work);

	//and put in the solutions
	for( int i = 0; i < i_nodes.size(); i++){
		Graph::vertex_descriptor v = i_nodes[i];
		graph_traits<Graph>::in_edge_iterator in_i = in_edges(v, tree->g).first;
		double l = gsl_vector_get(c, i);
		//if (l < 0) l = 1E-8;
		tree->g[*in_i].len = l;
	}

	//free memory
	gsl_multifit_linear_free(work);
	gsl_matrix_free(X);
	gsl_vector_free(y);
	gsl_vector_free(c);
	gsl_matrix_free(cov);

}

double GraphState::llik_normal(){
	double v = countdata->cov_var;
	double toreturn = 0;
	int npop = countdata->npop;
	for (int i = 0; i < npop; i++){
		for (int j = i; j < npop; j++){
			double pred = gsl_matrix_get(sigma, i, j);
			double obs = gsl_matrix_get(countdata->cov, i, j);
			double dif = obs-pred;
			double toadd = gsl_ran_gaussian_pdf(dif, sqrt(v));
			//cout << pred << " "<< obs << " "<< toadd << "\n";
			toreturn+= log(toadd);
		}
	}
	return toreturn;
}
