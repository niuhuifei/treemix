/*
 * GraphState2.cpp
 *
 *  Created on: Jun 28, 2011
 *      Author: pickrell
 */

#include "GraphState2.h"


GraphState2::GraphState2(CountData* counts, PhyloPop_params* pa){
	params = pa;
	countdata = counts;
	allpopnames = counts->list_pops();
	vector<string> startpops;
	startpops.push_back(allpopnames[0]); startpops.push_back(allpopnames[1]); startpops.push_back(allpopnames[2]);

	tree = new PopGraph(startpops);
	tree_bk = new PopGraph(startpops);

	current_npops = 3;
	sigma = gsl_matrix_alloc(current_npops, current_npops);
	sigma_cor = gsl_matrix_alloc(current_npops, current_npops);
	scatter = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(sigma);

	//tree->print();
	//cout << "llik: "<< llik_normal() << "\n";
	//cout << tree->get_newick_format() << "\n";
	set_branches_ls();
	//cout << tree->get_newick_format() << "\n";
	//cout << "here2\n"; cout.flush();
	compute_sigma();
	//cout << "here\n"; cout.flush();
	process_scatter();
	//cout << "not here\n"; cout.flush();
	current_llik = llik();

	//vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal(3);

	local_hillclimb(3);
	cout << "Starting from:\n";
	cout << tree->get_newick_format() << "\n";
}


void GraphState2::compute_sigma(){
	map<string, Graph::vertex_descriptor > popname2tip = tree->get_tips(tree->root);

	for (int i = 0; i < current_npops; i++){
		for (int j = i; j < current_npops; j++){

			string p1 = allpopnames[i];
			string p2 = allpopnames[j];

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

void GraphState2::print_sigma(){
	for(int i = 0; i < current_npops; i++){
		cout << allpopnames[i] << " ";
		for (int j = 0; j < current_npops; j++){
			cout << gsl_matrix_get(sigma, i, j) << " ";
		}
		cout << "\n";
	}
	cout << "\n";
	for(int i = 0; i < current_npops; i++){
		cout << allpopnames[i] << " ";
		for (int j = 0; j < current_npops; j++){
			cout << gsl_matrix_get(sigma_cor, i, j) << " ";
		}
		cout << "\n";
	}
}

void GraphState2::set_graph(string newick){
	tree->set_graph(newick);
	compute_sigma();
}

void GraphState2::set_branches_ls_old(){

	vector<Graph::vertex_descriptor> i_nodes = tree->get_inorder_traversal_noroot(current_npops);


	//initialize the workspace
	int n = current_npops * current_npops; //n is the total number of entries in the covariance matrix
	int p = 2*current_npops -2; // p is the number of branches lengths to be estimated
	int total = countdata->npop; //total is the total number of populations (for the bias correction)
	double inv_total = 1.0/ (double) total;
	double inv_total2 = 1.0/ ( (double) total * (double) total);

	//cout << "inv "<< inv_total << " "<< inv_total2 << "\n";
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, p);
	gsl_matrix * X = gsl_matrix_alloc(n, p); //X holds the weights on each branch. In the simplest model, this is 1s and 0s corresponding to whether the branch
	                                         // is in the path to the root. With the bias correction, will contain terms with 1/n and 1/n^2, where n is the total number of populations
	gsl_vector * y  = gsl_vector_alloc(n);  // y contains the entries of the empirical covariance matrix
	gsl_vector * c = gsl_vector_alloc(p);   // c will be estimated, contains the entries of the fitted covariance matrix
	gsl_matrix * cov = gsl_matrix_alloc(p, p);
	double chisq;
	gsl_matrix_set_zero(X);
	map<string, Graph::vertex_descriptor> popname2tip = tree->get_tips(tree->root);
	//for (map<string, Graph::vertex_descriptor>::iterator it = popname2tip.begin(); it != popname2tip.end(); it++){
	//	cout << it->first << " "<< tree->g[it->second].index << "\n";
	//}


	//set up the workspace


	// Set up the matrices from the tree
	int index = 0;
	for( int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];
			//cout << p1 << " "<< p2 << "\n";cout.flush();
			double empirical_cov = countdata->get_cov(p1, p2);
			//cout << p1 << " "<< p2 << " "<< empirical_cov << "\n";cout.flush();
			set<Graph::vertex_descriptor> path;

			if (i == j) path = tree->get_path_to_root(popname2tip[p1]);

			else{
				//cout << p1 << " "<< p2 << " "<< tree->g[popname2tip[p2]].index << "\n";
				Graph::vertex_descriptor lca = tree->get_LCA(tree->root, popname2tip[p1], popname2tip[p2]);
				path = tree->get_path_to_root(lca);

			}

			gsl_vector_set(y, index, empirical_cov);
			for (int i2 =0 ; i2 < i_nodes.size(); i2++){
				Graph::vertex_descriptor tmp = i_nodes[i2];
				if  (path.find(tmp) != path.end() ) {
					gsl_matrix_set(X, index, i2, gsl_matrix_get(X, index, i2)+1);
					if (params->bias_correct){
						for (int z = 0; z < n ; z++){
							//if (z == 0) cout << i2 << "\n";
							gsl_matrix_set(X, z, i2, gsl_matrix_get(X, z, i2) + inv_total2);
						}
					}
				}
			}
			index++;
		}
	}

	//for(int i = 0; i < n ; i ++){
	//	cout << gsl_vector_get(y, i) << " ";
	//	for (int j = 0; j < p; j++){
	//		cout << gsl_matrix_get(X, i, j) << " ";
	//	}
	//	cout << "\n\n";
	//}
	//cout << "\n";
	// Now add the bias terms
	index = 0;
	if (params->bias_correct){
		for( int i = 0; i < current_npops; i++){
			for (int j = 0; j < current_npops; j++){
				string p1 = allpopnames[i];
				string p2 = allpopnames[j];
				set<Graph::vertex_descriptor> path_ik;
				set<Graph::vertex_descriptor> path_jk;
				for (int k = 0; k < current_npops; k++){
					string p3 = allpopnames[k];
					if (i == k) path_ik = tree->get_path_to_root(popname2tip[p1]);
					else{
						Graph::vertex_descriptor lca = tree->get_LCA(tree->root, popname2tip[p1], popname2tip[p3]);
						path_ik = tree->get_path_to_root(lca);
					}
					if (j == k) path_jk = tree->get_path_to_root(popname2tip[p2]);
					else{
					//cout << p1 << " "<< p2 << " "<< tree->g[popname2tip[p2]].index << "\n";
						Graph::vertex_descriptor lca = tree->get_LCA(tree->root, popname2tip[p2], popname2tip[p3]);
						path_jk = tree->get_path_to_root(lca);
					}
					for (int i2 =0 ; i2 < i_nodes.size(); i2++){
						Graph::vertex_descriptor tmp = i_nodes[i2];
						if (i == j){
							if  (path_ik.find(tmp) != path_ik.end() ) {


								gsl_matrix_set(X, index, i2, gsl_matrix_get(X, index, i2) - 2*inv_total);
							//if (index == 0) cout << gsl_matrix_get(X, index, i2) << "\n";
							}
						}
						else{
							if  (path_ik.find(tmp) != path_ik.end() ) {
								gsl_matrix_set(X, index, i2, gsl_matrix_get(X, index, i2) - inv_total);
								//if (index == 1) cout << i2 << "\n";
							}
							if  (path_jk.find(tmp) != path_jk.end() ) {
								gsl_matrix_set(X, index, i2, gsl_matrix_get(X, index, i2) - inv_total);
								//if (index == 1) cout << i2 << "\n";
							}
						}

					}
				}
				index++;
			}
		}
	}

	//for(int i = 0; i < n ; i ++){
	//		cout << gsl_vector_get(y, i) << " ";
	//		for (int j = 0; j < p; j++){
	//			cout << gsl_matrix_get(X, i, j) << " ";
	//		}
	//		cout << "\n\n";
	//	}


	// fit the least squares estimates
	gsl_multifit_linear(X, y, c, cov, &chisq, work);

	//and put in the solutions in the graph
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


void GraphState2::set_branches_ls(){

	// edit so only one parameter for branch length including root

	vector<Graph::vertex_descriptor> i_nodes = tree->get_inorder_traversal_noroot(current_npops);
	set<Graph::vertex_descriptor> root_adj = tree->get_root_adj();
	vector<Graph::vertex_descriptor> i_nodes2;
	for (int i = 0; i < i_nodes.size(); i++){
		if (root_adj.find(i_nodes[i]) == root_adj.end()) i_nodes2.push_back( i_nodes[i] );
	}
	int joint_index = i_nodes2.size();

	//initialize the workspace
	int n = current_npops * current_npops; //n is the total number of entries in the covariance matrix
	int p = 2*current_npops -3; // p is the number of branches lengths to be estimated
	int total = countdata->npop; //total is the total number of populations (for the bias correction)
	double inv_total = 1.0/ (double) total;
	double inv_total2 = 1.0/ ( (double) total * (double) total);

	//cout << "inv "<< inv_total << " "<< inv_total2 << "\n";
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, p);
	gsl_matrix * X = gsl_matrix_alloc(n, p); //X holds the weights on each branch. In the simplest model, this is 1s and 0s corresponding to whether the branch
	                                         // is in the path to the root. With the bias correction, will contain terms with 1/n and 1/n^2, where n is the total number of populations
	gsl_vector * y  = gsl_vector_alloc(n);  // y contains the entries of the empirical covariance matrix
	gsl_vector * c = gsl_vector_alloc(p);   // c will be estimated, contains the entries of the fitted covariance matrix
	gsl_matrix * cov = gsl_matrix_alloc(p, p);
	double chisq;
	gsl_matrix_set_zero(X);
	map<string, Graph::vertex_descriptor> popname2tip = tree->get_tips(tree->root);
	//for (map<string, Graph::vertex_descriptor>::iterator it = popname2tip.begin(); it != popname2tip.end(); it++){
	//	cout << it->first << " "<< tree->g[it->second].index << "\n";
	//}


	//set up the workspace


	// Set up the matrices from the tree
	int index = 0;
	for( int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];
			//cout << p1 << " "<< p2 << "\n";cout.flush();
			double empirical_cov = countdata->get_cov(p1, p2);
			//cout << p1 << " "<< p2 << " "<< empirical_cov << "\n";cout.flush();
			set<Graph::vertex_descriptor> path;

			if (i == j) path = tree->get_path_to_root(popname2tip[p1]);

			else{
				//cout << p1 << " "<< p2 << " "<< tree->g[popname2tip[p2]].index << "\n";
				Graph::vertex_descriptor lca = tree->get_LCA(tree->root, popname2tip[p1], popname2tip[p2]);
				path = tree->get_path_to_root(lca);

			}

			gsl_vector_set(y, index, empirical_cov);

			for (int i2 =0 ; i2 < i_nodes2.size(); i2++){
				Graph::vertex_descriptor tmp = i_nodes2[i2];
				if  (path.find(tmp) != path.end() ) {
					gsl_matrix_set(X, index, i2, gsl_matrix_get(X, index, i2)+1);
					if (params->bias_correct){
						for (int z = 0; z < n ; z++){
							//if (z == 0) cout << i2 << "\n";
							gsl_matrix_set(X, z, i2, gsl_matrix_get(X, z, i2) + inv_total2);
						}
					}
				}
			}
			for (set<Graph::vertex_descriptor>::iterator it = root_adj.begin(); it != root_adj.end() ;it++){
				if (path.find(*it) != path.end() ){
					gsl_matrix_set(X, index, joint_index, gsl_matrix_get(X, index, joint_index)+0.5);
					if (params->bias_correct) {
						for (int z = 0; z < n ; z++){
							gsl_matrix_set(X, z, joint_index, gsl_matrix_get(X, z, joint_index) + inv_total2/2);
						}
					}
				}
			}
			index++;
		}
	}

	//for(int i = 0; i < n ; i ++){
	//	cout << gsl_vector_get(y, i) << " ";
	//	for (int j = 0; j < p; j++){
	//		cout << gsl_matrix_get(X, i, j) << " ";
	//	}
	//	cout << "\n\n";
	//}
	//cout << "\n";
	// Now add the bias terms
	index = 0;
	if (params->bias_correct){
		for( int i = 0; i < current_npops; i++){
			for (int j = 0; j < current_npops; j++){
				string p1 = allpopnames[i];
				string p2 = allpopnames[j];
				set<Graph::vertex_descriptor> path_ik;
				set<Graph::vertex_descriptor> path_jk;
				for (int k = 0; k < current_npops; k++){
					string p3 = allpopnames[k];
					if (i == k) path_ik = tree->get_path_to_root(popname2tip[p1]);
					else{
						Graph::vertex_descriptor lca = tree->get_LCA(tree->root, popname2tip[p1], popname2tip[p3]);
						path_ik = tree->get_path_to_root(lca);
					}
					if (j == k) path_jk = tree->get_path_to_root(popname2tip[p2]);
					else{
					//cout << p1 << " "<< p2 << " "<< tree->g[popname2tip[p2]].index << "\n";
						Graph::vertex_descriptor lca = tree->get_LCA(tree->root, popname2tip[p2], popname2tip[p3]);
						path_jk = tree->get_path_to_root(lca);
					}
					for (int i2 =0 ; i2 < i_nodes2.size(); i2++){
						Graph::vertex_descriptor tmp = i_nodes2[i2];
						if (i == j){
							if  (path_ik.find(tmp) != path_ik.end() ) {


								gsl_matrix_set(X, index, i2, gsl_matrix_get(X, index, i2) - 2*inv_total);
							//if (index == 0) cout << gsl_matrix_get(X, index, i2) << "\n";
							}
						}
						else{
							if  (path_ik.find(tmp) != path_ik.end() ) {
								gsl_matrix_set(X, index, i2, gsl_matrix_get(X, index, i2) - inv_total);
								//if (index == 1) cout << i2 << "\n";
							}
							if  (path_jk.find(tmp) != path_jk.end() ) {
								gsl_matrix_set(X, index, i2, gsl_matrix_get(X, index, i2) - inv_total);
								//if (index == 1) cout << i2 << "\n";
							}
						}

					}
					for (set<Graph::vertex_descriptor>::iterator it = root_adj.begin(); it != root_adj.end(); it++){
						if (i == j){
							if  (path_ik.find(*it) != path_ik.end() ) {


								gsl_matrix_set(X, index, joint_index, gsl_matrix_get(X, index, joint_index) - inv_total);
							//if (index == 0) cout << gsl_matrix_get(X, index, i2) << "\n";
							}
						}
						else{
							if  (path_ik.find(*it) != path_ik.end() ) {
								gsl_matrix_set(X, index, joint_index, gsl_matrix_get(X, index, joint_index) - inv_total/2);
								//if (index == 1) cout << i2 << "\n";
							}
							if  (path_jk.find(*it) != path_jk.end() ) {
								gsl_matrix_set(X, index, joint_index, gsl_matrix_get(X, index, joint_index) - inv_total/2);
								//if (index == 1) cout << i2 << "\n";
							}
						}
					}
				}
				index++;
			}

		}
	}

	//for(int i = 0; i < n ; i ++){
	//		cout << gsl_vector_get(y, i) << " ";
	//		for (int j = 0; j < p; j++){
	//			cout << gsl_matrix_get(X, i, j) << " ";
	//		}
	//		cout << "\n\n";
	//	}


	// fit the least squares estimates
	gsl_multifit_linear(X, y, c, cov, &chisq, work);

	//and put in the solutions in the graph
	for( int i = 0; i < i_nodes2.size(); i++){
		Graph::vertex_descriptor v = i_nodes2[i];
		graph_traits<Graph>::in_edge_iterator in_i = in_edges(v, tree->g).first;
		double l = gsl_vector_get(c, i);
		//if (l < 0) l = 1E-8;
		tree->g[*in_i].len = l;
	}
	for(set<Graph::vertex_descriptor>::iterator it = root_adj.begin(); it != root_adj.end(); it++){
		Graph::vertex_descriptor v = *it;
		graph_traits<Graph>::in_edge_iterator in_i = in_edges(v, tree->g).first;
		double l = gsl_vector_get(c, joint_index);
		//if (l < 0) l = 1E-8;
		tree->g[*in_i].len = l/2;
	}

	//and the corrected covariance matrix
	index = 0;
	for (int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			double pred = 0;
			for (int k = 0; k < p; k++) {
				pred += gsl_matrix_get(X, index, k) * gsl_vector_get(c, k);
			}

			//cout << index << " "<< i << " "<< j << " "<< pred << "\n";
			gsl_matrix_set(sigma_cor, i, j, pred);
			//cout << "not here\n";
			index++;
		}

	}


	//free memory
	gsl_multifit_linear_free(work);
	gsl_matrix_free(X);
	gsl_vector_free(y);
	gsl_vector_free(c);
	gsl_matrix_free(cov);

}


double GraphState2::llik_normal(){
	double toreturn = 0;
	for (int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];
			double pred = gsl_matrix_get(sigma_cor, i, j);
			double obs = countdata->get_cov(p1, p2);
			double var = countdata->get_cov_var(p1, p2);
			double dif = obs-pred;
			double toadd = gsl_ran_gaussian_pdf(dif, sqrt(var));
			//cout << i << " "<< j << " "<< obs << " "<< pred << " "<< var << "\n";
			toreturn+= log(toadd);
		}
	}
	return toreturn;
}

int GraphState2::local_hillclimb(int inorder_index){
	// if there was a rearrangement, return 1. otw 0.
	//

	vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal(current_npops);
	double llik1, llik2;
	//cout << tree->get_newick_format() << " "<< llik() << " "<< current_llik << " l0\n";
	tree_bk->copy(tree);


	tree->local_rearrange(inorder[inorder_index], 1);
	//cout << tree->get_newick_format() << "\n";
	set_branches_ls();

	//compute_sigma();
	llik1 =  llik();
	//cout << tree->get_newick_format() << " "<< llik1 << " l1\n";
	tree->copy(tree_bk);

	inorder = tree->get_inorder_traversal(current_npops);
	tree->local_rearrange(inorder[inorder_index], 2);

	set_branches_ls();
	//cout << tree->get_newick_format()  << "\n";
	//compute_sigma();
	llik2 =  llik();
	//cout << tree->get_newick_format() << " "<< llik2 << " l2\n";
	//cout << current_llik << " "<< llik1 << " "<< llik2 << "\n";
	tree->copy(tree_bk);
	inorder = tree->get_inorder_traversal(current_npops);
	if (current_llik < llik1 || current_llik < llik2){
		if (llik1 > llik2){
			tree->local_rearrange(inorder[inorder_index], 1);
			set_branches_ls();
			//compute_sigma();
			current_llik = llik1;
			return 1;
		}
		else{
			tree->local_rearrange(inorder[inorder_index], 2);
			//cout << tree->get_newick_format() << " "<< llik2 << " "<< llik() << "\n";
			set_branches_ls();
			//compute_sigma();
			current_llik = llik2;
			return 1;
		}
	}
	return 0;
}

int GraphState2::many_local_hillclimb(){
	int leninorder = 2*current_npops -1;
	int toreturn = 0;
	for (int i = 1; i < leninorder; i+=2){
		toreturn += local_hillclimb(i);
	}
	return toreturn;
}

void GraphState2::iterate_hillclimb(){
	int changes = many_local_hillclimb();
	while (changes > 0) changes = many_local_hillclimb();
}

void GraphState2::add_pop(){

	gsl_matrix_free(sigma);
	sigma = gsl_matrix_alloc(current_npops+1, current_npops+1);
	gsl_matrix_set_zero(sigma);
	gsl_matrix_free(sigma_cor);
	sigma_cor = gsl_matrix_alloc(current_npops+1, current_npops+1);
	gsl_matrix_set_zero(sigma_cor);

	current_npops++;
	process_scatter();
	current_npops--;

	string toadd = allpopnames[current_npops];
	cout << "Adding "<< toadd << "\n";
	vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal(current_npops);
	int max_index;
	double max_llik;
	tree_bk->copy(tree);

	for (int i = 0; i < inorder.size(); i++){

		Graph::vertex_descriptor tmp = tree->add_tip(inorder[i], toadd);
		//tree->print();
		current_npops++;
		//cout << "here\n"; cout.flush();
		set_branches_ls();

		//compute_sigma();


		double llk = llik();
		//cout << "testing to add "<< tree->get_newick_format()<< " "<< llik() << " "<< max_llik<<"\n";
		//cout << i << " "<< llk << "\n"; cout.flush();
		if (i == 0){
			max_index = i;
			max_llik = llk;
		}
		else if (llk > max_llik){
			max_index = i;
			max_llik = llk;
		}
		//cout << "removing "<< tree->g[tmp].index;
		tree->copy(tree_bk);
		//tree->remove_tip(tmp);
		//tree->print(); cout << "\n";
		current_npops--;
		inorder = tree->get_inorder_traversal(current_npops);
	}
	Graph::vertex_descriptor tmp = tree->add_tip(inorder[max_index], toadd);
	current_npops++;
	set_branches_ls();
	//compute_sigma();
	//cout << "added "<< tree->get_newick_format()<< " "<< llik() << " "<< max_llik <<"\n";
	current_llik = max_llik;
	//cout << "will process scatter\n"; cout.flush();

}

double GraphState2::llik(){
	return llik_normal();
	//return llik_wishart();
}

double GraphState2::llik_wishart(){
	// density of the wishart distribution with covariance matrix sigma, n = number of snps-1, p = number of populations
	// 		scatter matrix has been stored in countdata->scatter, the ln(determinant) is in countdata->scatter_det
	// 		and the ln of the relevant multiariate gamma is in scatter_gamma
	//
	// density is ( [n-p-1]/2 * scatter_det - [1/2] trace [sigma^-1 * scatter] - [np/2] ln(2) - [n/2] ln(det(sigma)) - scatter_gamma

	double toreturn = 0;
	int s;
	double ld;
	int p = current_npops;
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
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, inv, scatter, 0.0, ViU);
	double trace = 0;
	for (int i = 0; i < p ; i++) trace+= gsl_matrix_get(ViU, i, i);

	toreturn+= ( (double) n- (double) p-1.0 )/2.0 * scatter_det - trace/2.0;
	toreturn+= -( (double) n* (double) p/2.0) * log(2.0);
	toreturn += -((double) n/2.0)*ld;
	toreturn+= -scatter_gamma;

	gsl_matrix_free(work);
	gsl_matrix_free(inv);
	gsl_matrix_free(ViU);
	gsl_permutation_free(perm);

	return toreturn;
}


void GraphState2::process_scatter(){
	scatter_det = 0;
	scatter_gamma = 0;

	int n = countdata->nsnp-1;
	size_t pop = current_npops;
	gsl_matrix_free(scatter);
	scatter = gsl_matrix_alloc(pop, pop);
	for (int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			//cout << i <<  " "<< j << "\n";
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];
			gsl_matrix_set(scatter, i, j, countdata->get_scatter(p1, p2));
		}
	}

	int s;

	gsl_matrix * work = gsl_matrix_alloc(pop,pop);
	gsl_permutation * p = gsl_permutation_alloc(pop);



	//do LU decomposition and get determinant
	gsl_linalg_LU_decomp( work, p, &s );
	scatter_det = gsl_linalg_LU_lndet( work );
	//get the log sum of the gammas
	scatter_gamma = ( (double) current_npops * ( (double)  current_npops-1.0) /4.0) * log (M_PI);
	for (int i = 1; i <= current_npops; i++) scatter_gamma+= gsl_sf_lngamma( (double) n/2.0 + (1.0- (double) i)/2.0);
	//cout << "scatter_gamma "<< scatter_gamma << "\n";
}


void GraphState2::print_sigma_cor(string outfile){
	ogzstream out(outfile.c_str());
	for (int i = 0; i < current_npops; i++) out << allpopnames.at(i)<< " ";
	out << "\n";
	for (int i = 0; i < current_npops; i++){
		out << allpopnames.at(i);
		for(int j = 0; j < current_npops; j++)	 out << " "<< gsl_matrix_get(sigma_cor, i, j);
		out << "\n";
	}

}

void GraphState2::add_mig(){
}

string GraphState2::get_trimmed_newick(){
	map<string, double> trim;
	string toreturn;

	for ( map<string, int>::iterator it = countdata->pop2id.begin(); it!= countdata->pop2id.end(); it++){
		int id = it->second;
		string pop = it->first;
		double meanhzy = countdata->mean_hzy.find(id)->second;
		double mean_n = countdata->mean_ninds.find(id)->second;
		double t = meanhzy / (4.0* mean_n);
		//cout << pop  << " "<< t << "\n";
		trim.insert(make_pair(pop, t));
	}
	toreturn = tree->get_newick_format(&trim);
	return toreturn;
}
