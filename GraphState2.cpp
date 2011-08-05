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
	tree_bk2 =new PopGraph(startpops);
	tree_bk3 = new PopGraph(startpops);

	current_npops = 3;
	sigma = gsl_matrix_alloc(current_npops, current_npops);
	sigma_cor = gsl_matrix_alloc(current_npops, current_npops);
	scatter = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(sigma);

	set_branches_ls();
	process_scatter();
	current_llik = llik();

	//vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal(3);

	local_hillclimb(3);
	cout << "Starting from:\n";
	cout << tree->get_newick_format() << "\n";
	phi = (1+sqrt(5))/2;
	resphi = 2-phi;
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
}


void GraphState2::set_graph_from_file(string infile){
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
    getline(in, st);
    current_npops = allpopnames.size();
	gsl_matrix_free(sigma);
	sigma = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(sigma);
	gsl_matrix_free(sigma_cor);
	sigma_cor = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(sigma_cor);

    cout << "Reading tree topology from file:\n";
    cout << st << "\n";
	tree->set_graph(st);

	map<string, Graph::vertex_descriptor> tips = tree->get_tips(tree->root);
	//cerr << "ERROR: Input Newick string with "<< tips.size() << " populations. Input file has "<< current_npops  <<"\n";
	//for (map<string, Graph::vertex_descriptor>::iterator it = tips.begin(); it != tips.end(); it++){
	//	cout << it->first << " "<< tree->g[it->second].name << "\n";
	//}
	if (tips.size()  != current_npops){
		cerr << "ERROR: Input Newick string with "<< tips.size() << " populations. Input file has "<< current_npops  <<"\n";
		exit(1);
	}
	for ( vector<string>::iterator it = allpopnames.begin(); it != allpopnames.end(); it++){
		if (tips.find(*it) == tips.end() ){
			cerr << "ERROR: No population "<< *it << " in Newick string\n";
			exit(1);
		}
	}
	//tree->print();
	set_branches_ls();
	cout << "ln(lk): "<< llik() << "\n";
}

map<Graph::vertex_descriptor, int> GraphState2::get_v2index(){

	map<Graph::vertex_descriptor, int> vertex2index;

	vector<Graph::vertex_descriptor> i_nodes = tree->get_inorder_traversal(current_npops); //get descriptors for all the nodes
	set<Graph::vertex_descriptor> root_adj = tree->get_root_adj(); //get the ones next to the root
	vector<Graph::vertex_descriptor> i_nodes2;  //remove the ones next to the root

	int index = 0;

	for (int i = 0; i < i_nodes.size(); i++){
		if (root_adj.find(i_nodes[i]) == root_adj.end() && !tree->g[ i_nodes[i] ].is_root) {
			i_nodes2.push_back( i_nodes[i] );
			vertex2index.insert(make_pair(i_nodes[i], index));
			index++;
		}
	}

	int joint_index = i_nodes2.size(); //the index of the parameter for the sum of the two branch lengths next to the root

	for(set<Graph::vertex_descriptor>::iterator it = root_adj.begin(); it != root_adj.end(); it++) vertex2index[*it] = joint_index;
	index++;
	//get all the migration nodes
	vector<Graph::vertex_descriptor> mig_nodes;
	for(Graph::vertex_iterator it = vertices(tree->g).first; it != vertices(tree->g).second; it++){
		if ( tree->g[*it].is_mig ) {
			mig_nodes.push_back(*it);
			vertex2index[*it] = index;
			index++;
		}
	}
	return vertex2index;

}
void GraphState2::set_branches_ls(){

	/* one parameter for each non-migration node (minus 1 for the root and 1 for the unidentifiable branch length next to the root), one for each migration nodde
	 *
	 *
    */

	map<Graph::vertex_descriptor, int> vertex2index;
	map<Graph::vertex_descriptor, float> vertex2frac;
	vector<Graph::vertex_descriptor> i_nodes = tree->get_inorder_traversal(current_npops); //get descriptors for all the nodes
	set<Graph::vertex_descriptor> root_adj = tree->get_root_adj(); //get the ones next to the root
	vector<Graph::vertex_descriptor> i_nodes2;  //remove the ones next to the root

	int index = 0;

	for (int i = 0; i < i_nodes.size(); i++){
		if (root_adj.find(i_nodes[i]) == root_adj.end() && !tree->g[ i_nodes[i] ].is_root) {
			i_nodes2.push_back( i_nodes[i] );
			vertex2index.insert(make_pair(i_nodes[i], index));
			vertex2frac.insert(make_pair(i_nodes[i], 1));
			index++;
		}
	}

	int joint_index = i_nodes2.size(); //the index of the parameter for the sum of the two branch lengths next to the root

	for(set<Graph::vertex_descriptor>::iterator it = root_adj.begin(); it != root_adj.end(); it++) vertex2index[*it] = joint_index;
	index++;
		//get all the migration nodes
	vector<Graph::vertex_descriptor> mig_nodes;
	for(Graph::vertex_iterator it = vertices(tree->g).first; it != vertices(tree->g).second; it++){
		if ( tree->g[*it].is_mig ) {
			mig_nodes.push_back(*it);
			vertex2index[*it] = index;
			index++;
		}
	}

	//initialize the workspace
	int n = current_npops * current_npops; //n is the total number of entries in the covariance matrix
	int p = 2*current_npops -3 + mig_nodes.size(); // p is the number of branches lengths to be estimated
	int total = countdata->npop; //total is the total number of populations (for the bias correction)
	double inv_total = 1.0/ (double) total;
	double inv_total2 = 1.0/ ( (double) total * (double) total);

	//set up the workspace
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, p);
	gsl_matrix * X = gsl_matrix_alloc(n, p); //X holds the weights on each branch. In the simplest model, this is 1s and 0s corresponding to whether the branch (weighted by the migration weights)
	                                         //   is in the path to the root. With the bias correction, will contain terms with 1/n and 1/n^2, where n is the total number of populations

	gsl_vector * y  = gsl_vector_alloc(n);  // y contains the entries of the empirical covariance matrix
	gsl_vector * c = gsl_vector_alloc(p);   // c will be estimated, contains the entries of the fitted covariance matrix
	gsl_matrix * cov = gsl_matrix_alloc(p, p);
	double chisq;
	gsl_matrix_set_zero(X);

	//get all paths to the root for all tips
	map<string, Graph::vertex_descriptor> popname2tip = tree->get_tips(tree->root);
	map<string, set<pair<double, set<Graph::vertex_descriptor> > > > name2paths;
	for( map<string, Graph::vertex_descriptor>::iterator it = popname2tip.begin(); it != popname2tip.end(); it++){
		set<pair<double, set<Graph::vertex_descriptor> > > tmpset = tree->get_paths_to_root(it->second);
		name2paths.insert(make_pair(it->first, tmpset));
	}


	// Set up the matrices from the tree
	index = 0;
	for( int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];

			double empirical_cov = countdata->get_cov(p1, p2);
			gsl_vector_set(y, index, empirical_cov);

			set<pair<double, set<Graph::vertex_descriptor> > > paths_1 = name2paths[p1];
			set<pair<double, set<Graph::vertex_descriptor> > > paths_2 = name2paths[p2];


			for(set<pair<double, set<Graph::vertex_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
				for(set<pair<double, set<Graph::vertex_descriptor> > >::iterator it2 = paths_2.begin(); it2 != paths_2.end(); it2++){
					for( set<Graph::vertex_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
						if (tree->g[*it3].is_root) continue;
						if (it2->second.find(*it3) != it2->second.end()){
							int addindex = vertex2index[*it3];
							double add = it1->first*it2->first;
							double inv2add = it1->first*it2->first*inv_total2;
							if ( root_adj.find(*it3) != root_adj.end()) {
									add = add/2;
									inv2add = inv2add/2;
							}
							gsl_matrix_set(X, index, addindex, gsl_matrix_get(X, index, addindex)+add);
							for (int z = 0; z < n ; z++){
								gsl_matrix_set(X, z, addindex, gsl_matrix_get(X, z, addindex) + inv2add);
							}
						}
					}
				}
			}
			index++;
		}
	}
/*
	for(int i = 0; i < n ; i ++){
		cout << gsl_vector_get(y, i) << " ";
		for (int j = 0; j < p; j++){
			cout << gsl_matrix_get(X, i, j) << " ";
		}
		cout << "\n";
	}
	cout << "\n";
*/
	// Now add the bias terms
	index = 0;
	for( int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];
			set<pair<double, set<Graph::vertex_descriptor> > > paths_1 = name2paths[p1];
			set<pair<double, set<Graph::vertex_descriptor> > > paths_2 = name2paths[p2];

			for (int k = 0; k < current_npops; k++){
				string p3 = allpopnames[k];
				set<pair<double, set<Graph::vertex_descriptor> > > paths_3 = name2paths[p3];
				for(set<pair<double, set<Graph::vertex_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
					for(set<pair<double, set<Graph::vertex_descriptor> > >::iterator it2 = paths_3.begin(); it2 != paths_3.end(); it2++){
						for( set<Graph::vertex_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
							if (tree->g[*it3].is_root) continue;
							if (it2->second.find(*it3) != it2->second.end()){
								int addindex = vertex2index[*it3];
								double weight = it1->first*it2->first;
								double invadd = weight*inv_total;
								if ( root_adj.find(*it3) != root_adj.end()) invadd = invadd/2;
								gsl_matrix_set(X, index, addindex, gsl_matrix_get(X, index, addindex)-invadd);
							}
						}
					}
				}

				for(set<pair<double, set<Graph::vertex_descriptor> > >::iterator it1 = paths_2.begin(); it1 != paths_2.end(); it1++){
					for(set<pair<double, set<Graph::vertex_descriptor> > >::iterator it2 = paths_3.begin(); it2 != paths_3.end(); it2++){
						for( set<Graph::vertex_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
							if (tree->g[*it3].is_root) continue;
							if (it2->second.find(*it3) != it2->second.end()){
								int addindex = vertex2index[*it3];
								double weight = it1->first*it2->first;
								double invadd = weight*inv_total;
								if ( root_adj.find(*it3) != root_adj.end()) invadd = invadd/2;
								gsl_matrix_set(X, index, addindex, gsl_matrix_get(X, index, addindex)-invadd);
							}
						}
					}
				}
			}
			index++;
		}

	}

/*
		for(int i = 0; i < n ; i ++){
			cout << gsl_vector_get(y, i) << " ";
			for (int j = 0; j < p; j++){
				cout << gsl_matrix_get(X, i, j) << " ";
			}
			cout << "\n";
		}
		cout << "\n\n";
*/




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

void GraphState2::set_branch_coefs(gsl_matrix* X, gsl_vector* y, map<Graph::edge_descriptor, int>* edge2index, map<Graph::edge_descriptor, double>* edge2frac){

	//
	// y = Xc
	//
	// y is the observed matrix, X contains the contribution of each branch length to each entry in y
	//
	//

	gsl_matrix_set_zero(X);
	// get all the paths to the root from each tip

	map<string, Graph::vertex_descriptor> popname2tip = tree->get_tips(tree->root);
	map<string, set<pair<double, set<Graph::edge_descriptor> > > > name2paths;
	for( map<string, Graph::vertex_descriptor>::iterator it = popname2tip.begin(); it != popname2tip.end(); it++){
		set<pair<double, set<Graph::edge_descriptor> > > tmpset = tree->get_paths_to_root_edge(it->second);
		name2paths.insert(make_pair(it->first, tmpset));
	}

	double inv_total = 1.0/ (double) countdata->npop;
	double inv_total2 = 1.0/ ( (double) countdata->npop * (double) countdata->npop);

	// Set up the matrices from the tree
	int index = 0;
	for( int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];

			double empirical_cov = countdata->get_cov(p1, p2);
			gsl_vector_set(y, index, empirical_cov);

			set<pair<double, set<Graph::edge_descriptor> > > paths_1 = name2paths[p1];
			set<pair<double, set<Graph::edge_descriptor> > > paths_2 = name2paths[p2];


			for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
				for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it2 = paths_2.begin(); it2 != paths_2.end(); it2++){
					for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
						if (it2->second.find(*it3) != it2->second.end()){
							int addindex = edge2index->find(*it3)->second;
							double frac = edge2frac->find(*it3)->second;
							double add = it1->first * it2->first *frac;
							double inv2add = add*inv_total2;
							gsl_matrix_set(X, index, addindex, gsl_matrix_get(X, index, addindex)+add);
							for (int z = 0; z < current_npops*current_npops ; z++){
								gsl_matrix_set(X, z, addindex, gsl_matrix_get(X, z, addindex) + inv2add);
							}
						}
					}
				}
			}
			index++;
		}
	}
/*
	for(int i = 0; i < n ; i ++){
		cout << gsl_vector_get(y, i) << " ";
		for (int j = 0; j < p; j++){
			cout << gsl_matrix_get(X, i, j) << " ";
		}
		cout << "\n";
	}
	cout << "\n";
*/
	// Now add the bias terms
	index = 0;
	for( int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];
			set<pair<double, set<Graph::edge_descriptor> > > paths_1 = name2paths[p1];
			set<pair<double, set<Graph::edge_descriptor> > > paths_2 = name2paths[p2];

			for (int k = 0; k < current_npops; k++){
				string p3 = allpopnames[k];
				set<pair<double, set<Graph::edge_descriptor> > > paths_3 = name2paths[p3];
				for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_1.begin(); it1 != paths_1.end(); it1++){
					for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it2 = paths_3.begin(); it2 != paths_3.end(); it2++){
						for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
							if (it2->second.find(*it3) != it2->second.end()){
								int addindex = edge2index->find(*it3)->second;
								double frac = edge2frac->find(*it3)->second;
								double weight = it1->first*it2->first*frac;
								double invadd = weight*inv_total;
								gsl_matrix_set(X, index, addindex, gsl_matrix_get(X, index, addindex)-invadd);
							}
						}
					}
				}

				for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it1 = paths_2.begin(); it1 != paths_2.end(); it1++){
					for(set<pair<double, set<Graph::edge_descriptor> > >::iterator it2 = paths_3.begin(); it2 != paths_3.end(); it2++){
						for( set<Graph::edge_descriptor>::iterator it3 = it1->second.begin(); it3 != it1->second.end(); it3++){
							if (it2->second.find(*it3) != it2->second.end()){
								int addindex = edge2index->find(*it3)->second;
								double frac = edge2frac->find(*it3)->second;
								double weight = it1->first*it2->first *frac;
								double invadd = weight*inv_total;
								gsl_matrix_set(X, index, addindex, gsl_matrix_get(X, index, addindex)-invadd);
							}
						}
					}
				}
			}
			index++;
		}

	}



}

void GraphState2::optimize_weights(){
	/*
	 *  Go through each migration edge, optimize (restricting to be in range [0,1]) by normal optimization using logistic function
	 *
	 */

	// get list of migration edges
	vector<Graph::edge_descriptor> mig_edges = tree->get_mig_edges();
	double start_llik = current_llik;
	bool done = false;
	int nit = 0;
	while(!done){
		for (vector<Graph::edge_descriptor>::iterator it = mig_edges.begin(); it != mig_edges.end(); it++){
			double min, max, guess;
			guess = tree->g[*it].weight;
			guess = exp(guess) / (1+exp(guess));
			min = params->minweight;
			max = params->maxweight;
			golden_section_weight( *it, min, guess, max, params->tau);
		}
		if (current_llik < start_llik+1E-8) done = true;
		else start_llik = current_llik;
		nit++;
	}

}

int GraphState2::golden_section_weight(Graph::edge_descriptor e, double min, double guess, double max, double tau){
	double x;


	if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
	else x = guess - resphi *(guess-min);
	if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
		double new_logweight = (min+max)/2;
		double neww = 1/ (1+exp(-new_logweight));
		tree->g[e].weight = neww;
		set_branches_ls_wmig();
		current_llik = llik();
		return 0;
	}
	double w = 1/(1+exp(-x));
	tree->g[e].weight = w;
	set_branches_ls_wmig();
	double f_x = -llik();

	w = 1/(1+exp(-guess));
	tree->g[e].weight = w;
	set_branches_ls_wmig();
	double f_guess = -llik();

	if (f_x < f_guess){
		if ( (max-guess) > (guess-min) ) return golden_section_weight(e, guess, x, max, tau);
		else return golden_section_weight(e, min, x, guess, tau);
	}
	else{
		if ( (max - guess) > (guess - min)  ) return golden_section_weight(e, min, guess, x, tau);
		else return golden_section_weight(e, x, guess, max, tau);
	}
}


void GraphState2::set_branches_ls_wmig(){

	/* one parameter for each non-migration node (minus 1 for the root and 1 for the unidentifiable branch length next to the root), one for each migration node

	   Complication when doing this: paths to root in terms of edges (migration coming into nodes makes nodes not possible).
	   Many edge lengths are not identifiable, so have a single parameter which is their sum. Need to figure out which edge goes with which parameter, how to weight them.
    */
	map<Graph::vertex_descriptor, int> vertex2index;
	vector<Graph::vertex_descriptor> i_nodes = tree->get_inorder_traversal(current_npops); //get descriptors for all the nodes
	set<Graph::vertex_descriptor> root_adj = tree->get_root_adj(); //get the ones next to the root
	vector<Graph::vertex_descriptor> i_nodes2;  //remove the ones next to the root

	int index = 0;

	for (int i = 0; i < i_nodes.size(); i++){
		if (root_adj.find(i_nodes[i]) == root_adj.end() && !tree->g[ i_nodes[i] ].is_root) {
			i_nodes2.push_back( i_nodes[i] );
			vertex2index.insert(make_pair(i_nodes[i], index));
			index++;
		}
	}

	int joint_index = i_nodes2.size(); //the index of the parameter for the sum of the two branch lengths next to the root

	for(set<Graph::vertex_descriptor>::iterator it = root_adj.begin(); it != root_adj.end(); it++) vertex2index[*it] = joint_index;
	index++;

	vector<Graph::vertex_descriptor> mig_nodes;
	for(Graph::vertex_iterator it = vertices(tree->g).first; it != vertices(tree->g).second; it++){
		if ( tree->g[*it].is_mig ) {
			mig_nodes.push_back(*it);
			vertex2index[*it] = index;
			index++;
		}
	}

	// now get the edge to index and edge to fraction maps
	map<Graph::edge_descriptor, int> edge2index;
	map<Graph::edge_descriptor, double> edge2frac;
	for( Graph::edge_iterator it = edges(tree->g).first; it != edges(tree->g).second; it++){
		Graph::vertex_descriptor index_vertex, t;
		double f = 1.0;
		t = target(*it, tree->g);
		if (tree->g[*it].is_mig) index_vertex = source(*it, tree->g);
		else if ( tree->g[t].is_mig) {
			index_vertex = tree->get_child_node_mig(t);
			f = tree->g[*it].len / tree->get_parent_node(index_vertex).second;
		}
		else index_vertex = t;
		if (root_adj.find(index_vertex) != root_adj.end()) f = f/2;
		int i;
		if (vertex2index.find(index_vertex) == vertex2index.end()){
			cerr << "Error in least squares estimation: vertex "<< tree->g[index_vertex].index << " not found in the list of vertices\n";
			exit(1);
		}
		else i = vertex2index[index_vertex];
		edge2index.insert(make_pair(*it, i));
		edge2frac.insert(make_pair(*it, f));
	}

	//initialize the workspace
	int n = current_npops * current_npops; //n is the total number of entries in the covariance matrix
	int p = 2*current_npops -3 + mig_nodes.size(); // p is the number of branches lengths to be estimated
	int total = countdata->npop; //total is the total number of populations (for the bias correction)
	double inv_total = 1.0/ (double) total;
	double inv_total2 = 1.0/ ( (double) total * (double) total);

	//set up the workspace
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, p);
	gsl_matrix * X = gsl_matrix_alloc(n, p); //X holds the weights on each branch. In the simplest model, this is 1s and 0s corresponding to whether the branch (weighted by the migration weights)
	                                         //   is in the path to the root. With the bias correction, will contain terms with 1/n and 1/n^2, where n is the total number of populations

	gsl_vector * y  = gsl_vector_alloc(n);  // y contains the entries of the empirical covariance matrix
	gsl_vector * c = gsl_vector_alloc(p);   // c will be estimated, contains the entries of the fitted covariance matrix
	gsl_matrix * cov = gsl_matrix_alloc(p, p);
	double chisq;
	gsl_matrix_set_zero(X);

	set_branch_coefs(X, y, &edge2index, &edge2frac);

	// fit the least squares estimates
	gsl_multifit_linear(X, y, c, cov, &chisq, work);

	//and put in the solutions in the graph
	for( Graph::edge_iterator it = edges(tree->g).first; it != edges(tree->g).second; it++){
		int i = edge2index[*it];
		double frac = edge2frac[*it];
		double l = gsl_vector_get(c, i);
		tree->g[*it].len = l*frac;
	}

	//and the corrected covariance matrix
	index = 0;
	for (int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			double pred = 0;
			for (int k = 0; k < p; k++) {
				pred += gsl_matrix_get(X, index, k) * gsl_vector_get(c, k);
			}
			gsl_matrix_set(sigma_cor, i, j, pred);
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
		for (int j = i; j < current_npops; j++){
			string p1 = allpopnames[i];
			string p2 = allpopnames[j];
			double pred = gsl_matrix_get(sigma_cor, i, j);
			double obs = countdata->get_cov(p1, p2);
			double var = countdata->get_cov_var(p1, p2);
			double dif = obs-pred;
			double toadd = gsl_ran_gaussian_pdf(dif, sqrt(var));
			//cout << i << " "<< j << " "<< obs << " "<< pred << " "<< var << " "<< toadd << "\n";
			toreturn+= log(toadd);
		}
	}
	//cout << tree->get_newick_format() << "\n";
	//cout << toreturn << "\n";
	return toreturn;
}

int GraphState2::local_hillclimb(int inorder_index){
	// if there was a rearrangement, return 1. otw 0.
	//

	vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal(current_npops);
	double llik1, llik2, llik3;

	tree_bk->copy(tree);


	llik1 = current_llik;

	tree->local_rearrange(inorder[inorder_index], 1);
	set_branches_ls();
	llik2 =  llik();
	tree_bk2->copy(tree);


	tree->copy(tree_bk);
	inorder = tree->get_inorder_traversal(current_npops);
	tree->local_rearrange(inorder[inorder_index], 2);

	set_branches_ls();
	llik3 =  llik();
	tree_bk3->copy(tree);
	tree->copy(tree_bk);
	if (llik1 < llik2 || llik1 < llik3){
		if (llik2 > llik3){
			tree->copy(tree_bk2);
			current_llik = llik2;
			return 1;
		}
		else{
			tree->copy(tree_bk3);
			current_llik = llik3;
			return 1;
		}
	}
	return 0;
}

int GraphState2::global_hillclimb(int inorder_index){
	// take the node in inorder_index, try putting it above all the other nodes
	// return 1 if there is an improvement, 0 otw
	double max = current_llik;
	int maxindex = inorder_index;
	tree_bk ->copy(tree);
	int toreturn = 0;
	for (int i = 0; i < 2*current_npops-1; i++){
		vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal(current_npops);
		Graph::vertex_descriptor v1 = inorder[inorder_index];
		Graph::vertex_descriptor v2 = inorder[i];

		//make sure its a reasonable move (can't attach a clade within itself)
		if (i == inorder_index) continue;
		if ( tree->get_parent_node( v1 ).first == tree->get_parent_node( v2 ).first) continue;
		if ( tree->get_parent_node(v1).first == v2) continue;

		set<Graph::vertex_descriptor> path = tree->get_path_to_root( v2 );
		if ( path.find(v1) != path.end() ) continue;

		tree->global_rearrange(v1, v2);
		set_branches_ls();
		double lk =  llik();
		if ( lk > max ){
			max = lk;
			maxindex = i;
			toreturn = 1;
		}
		tree->copy(tree_bk);
	}
	if (toreturn ==1){
		vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal(current_npops);
		tree->global_rearrange( inorder[ inorder_index], inorder[maxindex] );
		set_branches_ls();
		double lk =  llik();
		current_llik = lk;
		cout << tree->get_newick_format() <<"\n";
		cout << "ln(lk): "<< lk <<"\n"; cout.flush();
	}
	return toreturn;

}
int GraphState2::many_local_hillclimb(){
	int leninorder = 2*current_npops -1;
	int toreturn = 0;
	for (int i = 1; i < leninorder; i+=2){
		toreturn += local_hillclimb(i);
	}
	return toreturn;
}

int GraphState2::many_global_hillclimb(){
	int leninorder = 2*current_npops -1;
	int toreturn = 0;
	for (int i = 0; i < leninorder; i++){
		toreturn += global_hillclimb(i);
	}
	return toreturn;
}

void GraphState2::iterate_hillclimb(){
	int changes = many_local_hillclimb();
	//cout << "Hill climbing "<< changes << " changes\n";
	while (changes > 0) {
		changes = many_local_hillclimb();
		//cout << "Hill climbing "<< changes << " changes\n";
	}
}


void GraphState2::iterate_global_hillclimb(){
	int changes = many_global_hillclimb();
	while (changes > 0) changes = many_global_hillclimb();

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
		current_npops++;
		set_branches_ls();
		double llk = llik();
		if (i == 0){
			max_index = i;
			max_llik = llk;
		}
		else if (llk > max_llik){
			max_index = i;
			max_llik = llk;
		}
		tree->copy(tree_bk);
		current_npops--;
		inorder = tree->get_inorder_traversal(current_npops);
	}
	Graph::vertex_descriptor tmp = tree->add_tip(inorder[max_index], toadd);
	current_npops++;
	set_branches_ls();
	current_llik = max_llik;
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
	vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal(current_npops);
	int maxst;
	int maxsp;
	int size = inorder.size();
	tree_bk->copy(tree);
	double maxllk = current_llik;

	for (int i = 0; i < size; i++){
		for (int j = 0; j < size; j++){

			inorder = tree->get_inorder_traversal(current_npops);
			if (tree->is_legal_migration(inorder[i], inorder[j])){
				Graph::edge_descriptor e = tree->add_mig_edge( inorder[i], inorder[j]);
				optimize_weights();
				//Graph::vertex_descriptor v = target(e, tree->g);
				//Graph::in_edge_iterator it = in_edges(v, tree->g).first;
				//bool neg = false;
				//while (it != in_edges(v, tree->g).second){
				//	if (tree->g[*it].len < 0)	neg = true;
				//	it++;
				//}
				//if (neg) {
				//	tree->copy(tree_bk);
				//	continue;
				// }
				if (current_llik > maxllk){
					maxst = i;
					maxsp = j;
					maxllk = current_llik;
				}

			}
			tree->copy(tree_bk);
		}
	}
	cout << "here2\n"; cout.flush();
	cout << maxst <<  " "<< maxsp << " here\n";
	inorder = tree->get_inorder_traversal(current_npops);
	cout << tree->g[ inorder[maxst]].index <<  " "<< tree->g[ inorder[maxsp]].index<< " here2\n";
	tree->add_mig_edge(inorder[maxst], inorder[maxsp]);
	optimize_weights();
}

bool GraphState2::add_mig_targeted(){
	// find the largest residual, try migration events in the vicinity
	// return true if an event is added, false otw

	//1. find the largest residual
	bool toreturn =false;
	string pop1, pop2;
	double max = 0;
	for (int i = 0; i < current_npops; i++){
		for (int j = i; j < current_npops; j++){
			if( i == j) continue;
			double cov = countdata->get_cov( allpopnames[i], allpopnames[j] );
			double fitted = gsl_matrix_get(sigma_cor, i, j);
			double diff = cov-fitted;
			if (diff > max){
				pop1 = allpopnames[i];
				pop2 = allpopnames[j];
			}
		}
	}
	cout << "Targeting migration to "<< pop1 <<" and "<< pop2 << "\n"; cout.flush();
	//2. Get the paths to the root for both
	map<string, Graph::vertex_descriptor> p2node = tree->get_tips(tree->root);
	set< pair< double, set<Graph::vertex_descriptor> > > p1 = tree->get_paths_to_root(p2node[pop1]);
	set< pair< double, set<Graph::vertex_descriptor> > > p2 = tree->get_paths_to_root(p2node[pop2]);
	set<Graph::vertex_descriptor> p1_s, p2_s;
	for (set< pair< double, set<Graph::vertex_descriptor> > >::iterator it = p1.begin(); it != p1.end(); it++){
		for (set<Graph::vertex_descriptor>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
			if (!tree->g[*it2].is_mig) p1_s.insert(*it2);
		}
	}
	for (set< pair< double, set<Graph::vertex_descriptor> > >::iterator it = p2.begin(); it != p2.end(); it++){
		for (set<Graph::vertex_descriptor>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) {
			if (!tree->g[*it2].is_mig) p2_s.insert(*it2);
		}
	}

	//3. try migration to all the pairwise combinations
	double max_llik = current_llik;
	tree_bk->copy(tree);
	Graph::vertex_descriptor max_1, max_2;
	for (set<Graph::vertex_descriptor>::iterator it = p1_s.begin(); it != p1_s.end(); it++){
		for (set<Graph::vertex_descriptor>::iterator it2 = p2_s.begin(); it2 != p2_s.end(); it2++){
			if ( tree->is_legal_migration(*it, *it2)){
				Graph::edge_descriptor e = tree->add_mig_edge( *it, *it2);
				optimize_weights();
				if (current_llik > max_llik){
					max_1 = *it;
					max_2 = *it2;
					max_llik = current_llik;
				}
				//remove_vertex( source(e, tree->g), tree->g);
				//remove_edge(e, tree->g);
			}
			if ( tree->is_legal_migration(*it2, *it)){
				Graph::edge_descriptor e = tree->add_mig_edge( *it2, *it);
				optimize_weights();
				if (current_llik > max_llik){
					max_1 = *it;
					max_2 = *it2;
					max_llik = current_llik;
				}
				//remove_vertex( source(e, tree->g), tree->g);
				//remove_edge(e, tree->g);
			}
		}
	}
	if (toreturn == true)	{
		tree->add_mig_edge(max_1, max_2);
		optimize_weights();
	}
	return toreturn;

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
