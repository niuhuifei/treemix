/*
 * GraphState2.cpp
 *
 *  Created on: Jun 28, 2011
 *      Author: pickrell
 */

#include "GraphState2.h"


GraphState2::GraphState2(CountData* counts){
	countdata = counts;
	allpopnames = counts->list_pops();
	vector<string> startpops;
	startpops.push_back(allpopnames[0]); startpops.push_back(allpopnames[1]); startpops.push_back(allpopnames[2]);

	tree = new PopGraph(startpops);
	tree_bk = new PopGraph(startpops);

	current_npops = 3;
	sigma = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(sigma);

	tree->print();
	//cout << "llik: "<< llik_normal() << "\n";
	set_branches_ls();
	compute_sigma();
	current_llik = llik_normal();
	vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal(3);

	local_hillclimb(3);
}


void GraphState2::compute_sigma(){
	map<string, Graph::vertex_descriptor > popname2tip = tree->get_tips(tree->root);
	//tree->print();
	for (int i = 0; i < current_npops; i++){
		for (int j = i; j < current_npops; j++){
			//cout << i << " "<< j << "\n";
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];
			//cout << p1 << " "<< p2 << "\n";
			if (i == j){
				double dist = tree->get_dist_to_root(popname2tip[p1]);
				gsl_matrix_set(sigma, i, j, dist);
			}
			else{

				Graph::vertex_descriptor lca = tree->get_LCA(tree->root, popname2tip[p1], popname2tip[p2]);
				//cout << p1 << " "<< p2 << " "<< tree->g[lca].index << "\n";
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
}



void GraphState2::set_branches_ls(){
	vector<Graph::vertex_descriptor> i_nodes = tree->get_inorder_traversal_noroot(current_npops);


	//initialize the workspace
	int n = current_npops * current_npops; //n is the total number of entries in the covariance matrix
	int p = 2*current_npops -2; // p is the number of branches lengths to be estimated
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, p);
	gsl_matrix * X = gsl_matrix_alloc(n, p);
	gsl_vector * y  = gsl_vector_alloc(n);
	gsl_vector * c = gsl_vector_alloc(p);
	gsl_matrix * cov = gsl_matrix_alloc(p, p);
	double chisq;

	map<string, Graph::vertex_descriptor> popname2tip = tree->get_tips(tree->root);
	//set up the workspace

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



double GraphState2::llik_normal(){
	double toreturn = 0;
	for (int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];
			double pred = gsl_matrix_get(sigma, i, j);
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

void GraphState2::local_hillclimb(int inorder_index){
	vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal(current_npops);
	double llik1, llik2;
	tree_bk->copy(tree);

	tree->local_rearrange(inorder[inorder_index], 1);
	set_branches_ls();
	compute_sigma();
	llik1 =  llik_normal();
	tree->copy(tree_bk);

	inorder = tree->get_inorder_traversal(current_npops);
	tree->local_rearrange(inorder[inorder_index], 2);
	set_branches_ls();
	compute_sigma();
	llik2 =  llik_normal();

	//cout << current_llik << " "<< llik1 << " "<< llik2 << "\n";
	tree->copy(tree_bk);
	inorder = tree->get_inorder_traversal(current_npops);
	if (current_llik < llik1 || current_llik < llik2){
		if (llik1 > llik2){
			tree->local_rearrange(inorder[inorder_index], 1);
			set_branches_ls();
			compute_sigma();
			current_llik = llik1;
		}
		else{
			tree->local_rearrange(inorder[inorder_index], 2);
			set_branches_ls();
			compute_sigma();
			current_llik = llik2;
		}
	}
}

void GraphState2::add_pop(){

	gsl_matrix_free(sigma);
	sigma = gsl_matrix_alloc(current_npops+1, current_npops+1);
	gsl_matrix_set_zero(sigma);
	string toadd = allpopnames[current_npops];
	cout << "Adding "<< toadd << "\n";
	vector<Graph::vertex_descriptor> inorder_noroot = tree->get_inorder_traversal_noroot(current_npops);
	int max_index;
	double max_llik;

	for (int i = 0; i < inorder_noroot.size(); i++){

		Graph::vertex_descriptor tmp = tree->add_tip(inorder_noroot[i], toadd);

		current_npops++;
		//tree->print();
		set_branches_ls();
		compute_sigma();
		double llk = llik_normal();

		cout << i << " "<< llk << "\n"; cout.flush();
		if (i == 0){
			max_index = i;
			max_llik = llk;
		}
		else if (llk > max_llik){
			max_index = i;
			max_llik = llk;
		}
		tree->remove_tip(tmp);
		current_npops--;
		inorder_noroot = tree->get_inorder_traversal_noroot(current_npops);
	}
	Graph::vertex_descriptor tmp = tree->add_tip(inorder_noroot[max_index], toadd);
	current_npops++;
	set_branches_ls();
	compute_sigma();
	current_llik = max_llik;
}
