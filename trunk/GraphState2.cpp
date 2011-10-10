/*
 * GraphState2.cpp
 *
 *  Created on: Jun 28, 2011
 *      Author: pickrell
 */

#include "GraphState2.h"

GraphState2::GraphState2(){
}

GraphState2::GraphState2(CountData* counts, PhyloPop_params* pa){
	params = pa;
	countdata = counts;
	allpopnames = counts->list_pops();
	unsigned int seed = unsigned( time(NULL));
	//seed = 1317332588;
	cout << "SEED: "<< seed << "\n";
	srand ( seed );
	random_shuffle(allpopnames.begin(), allpopnames.end() );
	for (int i = 0; i <  allpopnames.size(); i++)	popname2index.insert(make_pair( allpopnames[i], i));
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
	scatter_det = counts->scatter_det;
	scatter_gamma = counts->scatter_gamma;
	gsl_matrix_set_zero(sigma);

	set_branches_ls();
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

void GraphState2::set_graph(string vfile, string efile){
	igzstream vin(vfile.c_str());
	tree->g.clear();
    vector<string> line;
    struct stat stFileInfo;
    int intStat;
    current_npops = 0;
    tree->popnames.clear();

    intStat = stat(vfile.c_str(), &stFileInfo);
    if (intStat !=0){
            std::cerr<< "ERROR: cannot open file " << vfile << "\n";
            exit(1);
    }
    string st;
    while(getline(vin, st)){
    	string buf;
		stringstream ss(st);
		line.clear();
		while (ss>> buf){
			line.push_back(buf);
		}
		int index = atoi(line[0].c_str());
		if (index >= tree->indexcounter) tree->indexcounter = index+1;
		string name = line[1];
		string root = line[2];
		string mig = line[3];
		string tip = line[4];
		Graph::vertex_descriptor v= add_vertex(tree->g);
		tree->g[v].index = index;
		tree->g[v].name = name;
		tree->g[v].mig_frac = 1;
		if (root == "ROOT") {
			tree->g[v].is_root = true;
			tree->root = v;
		}
		else tree->g[v].is_root = false;
		if (mig == "MIG") {
			tree->g[v].is_mig = true;
			tree->g[v].mig_frac = 0.5;
		}
		else tree->g[v].is_mig = false;
		if (tip == "TIP") {
			tree->g[v].is_tip = true;
			current_npops++;
			tree->popnames.push_back(name);
		}
		else tree->g[v].is_tip = false;
		tree->g[v].rev = false;
     }

	igzstream ein(efile.c_str());
	map<int, Graph::vertex_descriptor> i2v = tree->index2vertex();
    intStat = stat(efile.c_str(), &stFileInfo);
    if (intStat !=0){
            std::cerr<< "ERROR: cannot open file " << efile << "\n";
            exit(1);
    }
    while(getline(ein, st)){
    	string buf;
		stringstream ss(st);
		line.clear();
		while (ss>> buf){
			line.push_back(buf);
		}
		int index1 = atoi(line[0].c_str());
		int index2 = atoi(line[1].c_str());
		float len = atof(line[2].c_str());
		float w = atof(line[3].c_str());
		string mig= line[4];
		Graph::edge_descriptor e = add_edge( i2v[index1], i2v[index2], tree->g).first;
		tree->g[e].len = len;
		tree->g[e].weight = w;
		if (mig == "MIG") {
			//cout << "here\n";
			tree->g[e].is_mig = true;
			//float mig_frac = atof(line[5].c_str());
			//tree->set_mig_frac(e, mig_frac);
			//cout << "not here\n";
		}
		else tree->g[e].is_mig = false;
    }
    ein.close();
    igzstream ein2(efile.c_str());
    while(getline(ein2, st)){
     	string buf;
 		stringstream ss(st);
 		line.clear();
 		while (ss>> buf){
 			line.push_back(buf);
 		}
 		int index1 = atoi(line[0].c_str());
 		int index2 = atoi(line[1].c_str());
 		string mig= line[4];
 		Graph::edge_descriptor e = edge( i2v[index1], i2v[index2], tree->g).first;
 		if (mig == "MIG") {
 			//cout << "here\n";
 			//tree->g[e].is_mig = true;
 			float mig_frac = atof(line[5].c_str());
 			tree->set_mig_frac(e, mig_frac);
 			//cout << "not here\n";
 		}
     }
	gsl_matrix_free(sigma);
	sigma = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(sigma);
	gsl_matrix_free(sigma_cor);
	sigma_cor = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(sigma_cor);

    set_branches_ls_wmig();

    current_llik = llik();
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
    cout << st << "\n"; cout.flush();
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
	current_llik = llik();
	cout << "ln(lk): "<< current_llik << "\n";
}



void GraphState2::set_graph_from_string(string newick){
	set_graph(newick);
	map<string, Graph::vertex_descriptor> tips = tree->get_tips(tree->root);
	allpopnames.clear();
	for (map<string, Graph::vertex_descriptor>::iterator it = tips.begin(); it != tips.end(); it++){
		allpopnames.push_back(it->first);
		if (countdata->pop2id.find(it->first) == countdata->pop2id.end()){
			cerr << "ERROR: cannot find population "<< it->first<< "\n";
			exit(1);
		}

	}
    current_npops = allpopnames.size();
	gsl_matrix_free(sigma);
	sigma = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(sigma);
	gsl_matrix_free(sigma_cor);
	sigma_cor = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_set_zero(sigma_cor);

	//tree->print();
	set_branches_ls();
	current_llik = llik();
	cout << "ln(lk): "<< current_llik << "\n";
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
						if ( tree->g[*it3].is_mig) continue;
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
							if (tree->g[*it3].is_mig)continue;
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
								if (tree->g[*it3].is_mig)continue;
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
			cout << tree->g[*it].weight << "\n";
			}
		if (current_llik < start_llik+0.1) done = true;
		else start_llik = current_llik;
		nit++;
	}

}


void GraphState2::optimize_fracs(){
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
			guess = tree->g[source(*it, tree->g)].mig_frac;
			guess = exp(guess) / (1+exp(guess));
			min = params->minweight;
			max = params->maxweight;
			golden_section_frac( *it, min, guess, max, params->tau);
			cout << tree->g[source(*it, tree->g)].mig_frac << " frac\n";
			}
		if (current_llik < start_llik+0.1) done = true;
		else start_llik = current_llik;
		nit++;
	}

}

void GraphState2::quick_optimize_weight(Graph::edge_descriptor e){
	tree->g[e].weight = 0.1;
	set_branches_ls_wmig();
	double best = 0.1;
	double start_llik = llik();
	for (double w = 0.3; w < 1; w+=0.2){
		tree->g[e].weight = w;
		set_branches_ls_wmig();
		double test_llik = llik();
		if (test_llik > start_llik){
			start_llik  = test_llik;
			best = w;
		}
	}
	tree->g[e].weight = best;
	set_branches_ls_wmig();
	current_llik = llik();
}



void GraphState2::optimize_weight(Graph::edge_descriptor e){

	if (params->quick){
		quick_optimize_weight(e);
	}
	else{
		double start_llik = current_llik;
		bool done = false;
		int nit = 0;
		while(!done && nit < params->maxit){
			//cout << nit << "\n"; cout.flush();
			double min, max, guess;
			guess = tree->g[e].weight;
			guess = exp(guess) / (1+exp(guess));
			min = params->minweight;
			max = params->maxweight;
			golden_section_weight(e, min, guess, max, params->tau);
			if (current_llik < start_llik+1E-8) done = true;
			else start_llik = current_llik;
			//cout << guess << " "<< start_llik << " "<< current_llik << "\n";
			nit++;
		}
	}

}



void GraphState2::optimize_frac(Graph::edge_descriptor e){


		double start_llik = current_llik;
		bool done = false;
		int nit = 0;
		while(!done && nit < params->maxit){
			//cout << nit << "\n"; cout.flush();
			double min, max, guess;
			guess = tree->g[source(e, tree->g)].mig_frac;
			guess = exp(guess) / (1+exp(guess));
			min = params->minweight;
			max = params->maxweight;
			golden_section_frac(e, min, guess, max, params->tau);
			if (current_llik < start_llik+0.01) done = true;
			else start_llik = current_llik;
			//cout << guess << " "<< start_llik << " "<< current_llik << "\n";
			nit++;
		}


}
int GraphState2::golden_section_weight(Graph::edge_descriptor e, double min, double guess, double max, double tau){
	double x;

	//cout << guess << "\n"; cout.flush();
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
	//cout << "here\n"; cout.flush();
	set_branches_ls_wmig();
	//cout << "not here\n"; cout.flush();
	double f_x = -llik();

	w = 1/(1+exp(-guess));
	tree->g[e].weight = w;
	//cout << "here2\n"; cout.flush();
	set_branches_ls_wmig();
	//cout << "not here2\n"; cout.flush();
	double f_guess = -llik();

	if (f_x < f_guess){
		if ( (max-guess) > (guess-min) )	return golden_section_weight(e, guess, x, max, tau);

		else return golden_section_weight(e, min, x, guess, tau);
	}
	else{
		if ( (max - guess) > (guess - min)  ) return golden_section_weight(e, min, guess, x, tau);
		else return golden_section_weight(e, x, guess, max, tau);
	}
}


int GraphState2::golden_section_frac(Graph::edge_descriptor e, double min, double guess, double max, double tau){
	double x;

	//cout << guess << "\n"; cout.flush();
	if ( (max - guess) > (guess - min)) x = guess + resphi *( max - guess);
	else x = guess - resphi *(guess-min);
	if (fabs(max-min) < tau * (fabs(guess)+fabs(max))) {
		double new_logweight = (min+max)/2;
		double neww = 1/ (1+exp(-new_logweight));
		tree->set_mig_frac(e, neww);
		set_branches_ls_wmig();
		current_llik = llik();
		return 0;
	}
	double w = 1/(1+exp(-x));
	tree->set_mig_frac(e, w);
	//cout << "here\n"; cout.flush();
	set_branches_ls_wmig();
	//cout << "not here\n"; cout.flush();
	double f_x = -llik();

	w = 1/(1+exp(-guess));
	tree->set_mig_frac(e, w);
	//tree->g[e].weight = w;
	//cout << "here2\n"; cout.flush();
	set_branches_ls_wmig();
	//cout << "not here2\n"; cout.flush();
	double f_guess = -llik();

	if (f_x < f_guess){
		if ( (max-guess) > (guess-min) )	return golden_section_frac(e, guess, x, max, tau);

		else return golden_section_frac(e, min, x, guess, tau);
	}
	else{
		if ( (max - guess) > (guess - min)  ) return golden_section_weight(e, min, guess, x, tau);
		else return golden_section_frac(e, x, guess, max, tau);
	}
}

void GraphState2::set_branches_ls_wmig(){

	/* one parameter for each non-migration node (minus 1 for the root and 1 for the unidentifiable branch length next to the root)

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
	//cout << "here1\n"; cout.flush();
	vector<Graph::vertex_descriptor> mig_nodes;
	for(Graph::vertex_iterator it = vertices(tree->g).first; it != vertices(tree->g).second; it++){
		if ( tree->g[*it].is_mig ) {
			mig_nodes.push_back(*it);
			vertex2index[*it] = index;
			index++;
		}
	}
	//cout << "here2\n"; cout.flush();
	// now get the edge to index and edge to fraction maps
	map<Graph::edge_descriptor, int> edge2index;
	map<Graph::edge_descriptor, double> edge2frac;
	for( Graph::edge_iterator it = edges(tree->g).first; it != edges(tree->g).second; it++){
		Graph::vertex_descriptor index_vertex, t;
		double f = 1.0;
		double f2 = 1.0;
		t = target(*it, tree->g);
		Graph::vertex_descriptor t2 = source(*it, tree->g);
		//cout << tree->g[t].index << "\n";
		if (tree->g[*it].is_mig) continue;
		else if ( tree->g[t].is_mig) {
			index_vertex = tree->get_child_node_mig(t);
		}
		else index_vertex = t;

		if (tree->g[t].is_mig || tree->g[t2].is_mig){
			if (tree->g[t].is_mig && tree->g[t2].is_mig)	{
				//cout << "here?\n"; cout.flush();
				f = tree->g[t].mig_frac - tree->g[t2].mig_frac;
			}
			else if (tree->g[t].is_mig) {
				//cout << "here2 "<< tree->g[t].mig_frac << "\n";
				f = tree->g[t].mig_frac;
			}
			else if (tree->g[t2].is_mig) {
				//cout << "here3 "<< tree->g[t2].mig_frac << "\n";
				f = 1-tree->g[t2].mig_frac;
			}
		}

		//f = tree->g[*it].len / tree->get_parent_node(index_vertex).second;

		//cout << tree->g[t2].index << " "<< tree->g[t].index << " "<< f<< " "<< tree->g[*it].len << "\n";
		//cout << tree->g[source(*it, tree->g)].index << " "<< tree->g[t].index << " "<< f << " "<< f2 << "\n";
		//cout << tree->g[t].index << "\n";
		if (root_adj.find(index_vertex) != root_adj.end()) f = f/2;
		int i;
		if (vertex2index.find(index_vertex) == vertex2index.end()){
			cerr << "Error in least squares estimation: vertex "<< tree->g[index_vertex].index << " not found in the list of vertices\n";
			exit(1);
		}
		else i = vertex2index[index_vertex];
		//cout << tree->g[source(*it, tree->g)].index<< " "<< tree->g[target(*it, tree->g)].index << " "<<  tree->g[index_vertex].index << " "<< f<< "\n";
		edge2index.insert(make_pair(*it, i));
		edge2frac.insert(make_pair(*it, f));
	}
	//cout <<  "\n";
	//cout << "here3\n"; cout.flush();
	//initialize the workspace
	int n = current_npops * current_npops; //n is the total number of entries in the covariance matrix
	int p = 2*current_npops -3; // p is the number of branches lengths to be estimated
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
	//cout << "here\n";
	set_branch_coefs(X, y, &edge2index, &edge2frac);
	//cout << "here2\n";
	// fit the least squares estimates
	gsl_multifit_linear(X, y, c, cov, &chisq, work);
	//and put in the solutions in the graph
	for( Graph::edge_iterator it = edges(tree->g).first; it != edges(tree->g).second; it++){
		if (tree->g[*it].is_mig) continue;
		int i = edge2index[*it];
		double frac = edge2frac[*it];
		double l = gsl_vector_get(c, i);
		tree->g[*it].len = l*frac;
	}
	//cout << "here3\n";
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

	//cout << "here4\n";
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
			double se = countdata->get_cov_var(p1, p2);
			double dif = obs-pred;
			double toadd = gsl_ran_gaussian_pdf(dif, se);
			if (toadd < 1e-300) 	toadd = 1e-300;

			//cout << pred << " "<< obs << " "<< se << " "<< toadd << "\n";
			toreturn+= log(toadd);
		}
	}
	//cout << "tmp "<< (double) tmp / (double) total <<"\n";
	return toreturn;
}

int GraphState2::local_hillclimb(int inorder_index){
	// if there was a rearrangement, return 1. otw 0.
	//

	vector<Graph::vertex_descriptor> inorder = tree->get_inorder_traversal(current_npops);
	if ( tree->g[ inorder[inorder_index]].is_root) return local_hillclimb_root();

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

int GraphState2::local_hillclimb_root(){
	double best_llik = current_llik;
	tree_bk->copy(tree);
	tree_bk2->copy(tree);
	int toreturn = 0;
	for (int i = 1; i <= 4; i++){
		tree->move_root(i);
		set_branches_ls();
		double tmplik = llik();
		if (tmplik > best_llik){
			toreturn = 1;
			tree_bk->copy(tree);
			best_llik = tmplik;
		}
		tree->copy(tree_bk2);
	}
	if (toreturn == 1){
		tree->copy(tree_bk);
		current_llik = best_llik;
	}
	return toreturn;
}

int GraphState2::global_hillclimb(int inorder_index){
	// take the node in inorder_index, try putting it above all the other nodes
	// return 1 if there is an improvement, 0 otw
	double max = current_llik;
	gsl_matrix* tmpfitted = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted, sigma_cor);

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
		if (!try_mig(v1, v2, tmpfitted)) continue;
		cout << tree->get_newick_format(v1)<< " "<< tree->get_newick_format(v2)<< "\n";
		tree->global_rearrange(v1, v2);

		set_branches_ls();
		double lk =  llik();
		cout << lk << " "<<  max << "\n";
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
	gsl_matrix_free(tmpfitted);
	return toreturn;

}

int GraphState2::iterate_local_hillclimb_wmig(int index){
	int moving = local_hillclimb_wmig(index);
	if (moving  == 0 ) return 0;
	while (moving > 0) {
		//cout << "moving vertex "<< index << " "<< moving << "\n";
		//cout << llik() << "\n";
		//cout << tree->get_newick_format()<<"\n";
		moving = local_hillclimb_wmig(index);
	}
	return 1;
}

void GraphState2::iterate_mig_hillclimb_and_optimweight(pair<int, int> indices){

	int moving1 = iterate_local_hillclimb_wmig(indices.first);
	//cout << "here!\n";
	//tree->print();
	//cout << "\n";
	int moving2 = iterate_local_hillclimb_wmig(indices.second);
	int moving = moving1+moving2;
	while (moving > 0){
		map<int, Graph::vertex_descriptor> vindex = tree->index2vertex();

		Graph::vertex_descriptor v = vindex[indices.second];
		set<Graph::edge_descriptor> inm = tree->get_in_mig_edges(v);
		//cout << "here\n"; cout.flush();
		for (set<Graph::edge_descriptor>::iterator it = inm.begin(); it != inm.end(); it++) optimize_weight(*it);
		//cout << "not here\n"; cout.flush();
		moving1 = iterate_local_hillclimb_wmig(indices.first);
		moving2 = iterate_local_hillclimb_wmig(indices.second);
		moving = moving1+moving2;
		//cout <<"or here\n"; cout.flush();

	}
}


int GraphState2::local_hillclimb_wmig(int index){
	double max = current_llik;
	double lik_bk = current_llik;
	tree_bk->copy(tree);
	tree_bk2->copy(tree);
	gsl_matrix *tmpfitted = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted, sigma_cor);
	int toreturn = 0;
	for (int i = 1; i <=4; i++){
		map<int, Graph::vertex_descriptor> index2v = tree->index2vertex();
		Graph::vertex_descriptor v = index2v[index];
		//if (i == 4 && index == 751) {
		//	tree->print("before_rearrange");
		//}
		bool rearr = tree->local_rearrange_wmig(v, i);
		if ( has_loop() ) {
			//cout << "has loop!\n";
			tree->copy(tree_bk);
			continue;
		}

		//tree->print();
		if (rearr){
			set_branches_ls_wmig();
			double lk = llik();
			//if (i == 4 && index == 751) tree->print("after_rearrange");
			//cout << lk << " "<< max << " "<< index << " "<< i << "\n";
			if (lk > max){
				max = lk;
				//cout << i << "\n";
				toreturn = 1;
				tree_bk2->copy(tree);
				gsl_matrix_memcpy( tmpfitted, sigma_cor);
			}
			tree->copy(tree_bk);
		}
	}
	if (toreturn == 1){
		tree->copy(tree_bk2);
		current_llik = max;
		//cout << lik_bk << " "<< max << "\n";
		//cout << "moving from"<< tree_bk->get_newick_format() <<"\n";
		//cout << "to "<< tree->get_newick_format()<<"\n";
		gsl_matrix_memcpy( sigma_cor, tmpfitted);
	}
	else{
		tree->copy(tree_bk);
		current_llik = lik_bk;
		gsl_matrix_memcpy( sigma_cor, tmpfitted);
	}
	gsl_matrix_free(tmpfitted);
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
	set_branches_ls();
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
	gsl_matrix_free(scatter);
	scatter = gsl_matrix_alloc(current_npops, current_npops);
	for (int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			//cout << allpopnames[i] << " "<< allpopnames[j] << " "<< countdata->get_scatter(allpopnames[i], allpopnames[j]) << "\n";
			gsl_matrix_set(scatter, i, j, countdata->get_scatter( allpopnames[i], allpopnames[j] ));
		}
	}

/*
	for (int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			cout << gsl_matrix_get(scatter, i, j) << " ";
		}
		cout << "\n";
	}
*/
	double toreturn = 0;
	int s;
	double ld;
	int p = current_npops;
	int n = countdata->nsnp -1;

	//set scatter gamma
	scatter_gamma = ( (double) (p-1) * ( (double)  (p-1)-1.0) /4.0) * log (M_PI);
	for (int i = 1; i <= p-1; i++) scatter_gamma+= gsl_sf_lngamma( (double) n/2.0 + (1.0- (double) i)/2.0);

	gsl_matrix * U = gsl_matrix_alloc(p-1, p);
	gsl_matrix *scatter_prime = gsl_matrix_alloc(p-1, p-1);
	gsl_matrix * W_prime = gsl_matrix_alloc(p-1, p-1);
	gsl_matrix * W_inv = gsl_matrix_alloc(p-1, p-1);

	// 1. Take SVD of scatter
	gsl_matrix * A = gsl_matrix_alloc(p,p);
	gsl_matrix * VT = gsl_matrix_alloc(p,p);
	gsl_vector * S = gsl_vector_alloc(p);
	gsl_vector * work = gsl_vector_alloc(p);
	gsl_matrix_memcpy( A, scatter );

	gsl_linalg_SV_decomp(A, VT, S, work);


	// Now copy the first npop=1 eigenvectors to U

	for (int i = 0; i < p-1; i++){
		for(int j = 0; j < p; j++){
			gsl_matrix_set(U, i, j, gsl_matrix_get(A, i, j));
		}
	}
/*
	cout << "\n";
	for (int i = 0; i < current_npops-1; i++){
		for (int j = 0; j < current_npops; j++){
			cout << gsl_matrix_get(U, i, j) << " ";
		}
		cout << "\n";
	}

	cout << "\n";
*/
	// 2. transform scatter into m-1 space
	// S' = U S U^T

	gsl_matrix * US = gsl_matrix_alloc(p-1, p);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, U, scatter, 0.0, US);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, US, U, 0.0, scatter_prime);
/*
	for (int i = 0; i < current_npops-1; i++){
			for (int j = 0; j < current_npops-1; j++){
				cout << gsl_matrix_get(scatter_prime, i, j) << " ";
			}
			cout << "\n";
		}
	cout <<"\n";

	*/
	// 3. Same thing on predicted covariance matrix
	gsl_matrix * UW = gsl_matrix_alloc(p-1, p);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, U, sigma_cor, 0.0, UW);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, UW, U, 0.0, W_prime);


	// 4. Do LU decomposition and get determinant of scatter_prime
	gsl_matrix_free(A);
	A = gsl_matrix_alloc(p-1, p-1);
	gsl_matrix_memcpy( A, scatter_prime );
	gsl_permutation * perm = gsl_permutation_alloc(p-1);
	gsl_linalg_LU_decomp( A, perm, &s );
	scatter_det = gsl_linalg_LU_lndet( A );

	// 5. Do LU decomposition, get inverse of W_prime
	gsl_matrix_free(A);
	gsl_permutation_free(perm);
	A = gsl_matrix_alloc(p-1, p-1);
	gsl_matrix_memcpy( A, W_prime );
	perm = gsl_permutation_alloc(p-1);
	gsl_linalg_LU_decomp( A, perm, &s );
	gsl_linalg_LU_invert( A, perm, W_inv );
	double w_det = gsl_linalg_LU_lndet( A );
/*
	cout << "scatter_prime_det "<< scatter_det << " w_det "<< w_det <<"\n";
	for (int i = 0; i < current_npops; i++){
		for (int j = 0; j < current_npops; j++){
			cout << gsl_matrix_get(sigma_cor, i, j) << " ";
		}
		cout << "\n";
	}

	cout<<"\n";

	for (int i = 0; i < current_npops-1; i++){
			for (int j = 0; j < current_npops-1; j++){
				cout << gsl_matrix_get(W_prime, i, j) << " ";
			}
			cout << "\n";
		}
	cout <<"\n";

	for (int i = 0; i < current_npops-1; i++){
			for (int j = 0; j < current_npops-1; j++){
				cout << gsl_matrix_get(W_inv, i, j) << " ";
			}
			cout << "\n";
		}
	cout <<"\n";

	*/
	//multiply inverse of cov by scatter, get trace
	gsl_matrix * ViU = gsl_matrix_alloc(p-1, p-1);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, W_inv, scatter_prime, 0.0, ViU);

	double trace = 0;
	for (int i = 0; i < p-1 ; i++) trace+= gsl_matrix_get(ViU, i, i);

	toreturn+= ( (double) n- (double) (p-1)-1.0 )/2.0 * scatter_det - trace/2.0;
	toreturn+= -( (double) n* (double) (p-1)/2.0) * log(2.0);
	toreturn += -((double) n/2.0)*w_det;
	toreturn+= -scatter_gamma;

	gsl_matrix_free( U );
	gsl_matrix_free(scatter_prime);
	gsl_matrix_free(W_prime);
	gsl_matrix_free(W_inv);
	gsl_matrix_free(A);
	gsl_matrix_free(VT);
	gsl_vector_free(S);
	gsl_vector_free(work);
	gsl_matrix_free(ViU);
	gsl_matrix_free( US );
	gsl_matrix_free( UW );
	gsl_permutation_free(perm);

	//cout << "llik: "<< toreturn <<"\n";
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

Graph::edge_descriptor GraphState2::add_mig(int index1, int index2){
	Graph::edge_descriptor e;
	map<int, Graph::vertex_descriptor> i2v = tree->index2vertex();
	if ( tree->is_legal_migration( i2v[index1], i2v[index2])){
		e = tree->add_mig_edge(i2v[index1], i2v[index2]);
		//quick_optimize_weight(e);
		optimize_weight(e);
		//quick_optimize_weight(e);
		Graph::vertex_descriptor v = source(e, tree->g);
		cout << tree->g[e].weight <<"\n";
		//for (float f = 0.1; f < 1; f+=0.2){
		//	tree->g[v].mig_frac = f;
			//optimize_weight(e);
		//	set_branches_ls_wmig();
		//	cout << f<< " "<< llik() << "\n";
		//}

	}
	else{
		cerr << "ERROR: not a legal migration between index " << index1 << " and "<< index2 << "\n";
		exit(1);
	}
	current_llik = llik();
	return(e);
}

void GraphState2::rearrange(int index1, int index2){
	map<int, Graph::vertex_descriptor> i2v = tree->index2vertex();
	tree->global_rearrange(i2v[index1], i2v[index2]);
	set_branches_ls_wmig();
	current_llik = llik();
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

pair<string, string> GraphState2::get_max_resid(){
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
				max = diff;
			}
		}
	}
	return make_pair(pop1, pop2);
}


pair<bool, pair<int, int> > GraphState2::add_mig_targeted(){
	// find the largest residual, try migration events in the vicinity
	// return true if an event is added, false otw

	// tmpfitted will have a backup of the current covariance matrix
	gsl_matrix *tmpfitted = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted, sigma_cor);

	//tmpfitted2 will hold the covariance matrix at the tree with the max likelihood
	gsl_matrix *tmpfitted2 = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted2, sigma_cor);

	double llik_bk = current_llik;
	//1. find the largest residual
	pair<bool, pair<int, int> > toreturn;
	toreturn.first = false;
	string pop1, pop2;
	double max = 0;
	int besti = 0;
	int bestj = 0;
	for (int i = 0; i < current_npops; i++){
		for (int j = i; j < current_npops; j++){
			if( i == j) continue;
			double cov = countdata->get_cov( allpopnames[i], allpopnames[j] );
			double fitted = gsl_matrix_get(sigma_cor, i, j);
			double diff = cov-fitted;
			//cout << allpopnames[i] << " "<< allpopnames[j] << " "<< diff << "\n";
			if (diff > max){
				pop1 = allpopnames[i];
				pop2 = allpopnames[j];
				besti = i;
				bestj = j;
				max = diff;
			}
		}
	}

	// try all pairwise combinations within 20% of the max residual
	double max_test = max*0.8;
	set<pair<int, int> > tested; //hold the pairs of vertices that have been tested
	cout << "Targeting migration to vicinity of "<< pop1 <<" and "<< pop2 << "\n"; cout.flush();
	set<pair<string, string> > pops2test;
	for (int i = 0; i < current_npops; i++){
		double cov1 = countdata->get_cov( allpopnames[besti], allpopnames[i]);
		double fitted1 = gsl_matrix_get(sigma_cor, besti, i);
		double diff = cov1-fitted1;
		if (diff >= max_test) pops2test.insert(make_pair( allpopnames[besti], allpopnames[i]));

		double cov2 = countdata->get_cov( allpopnames[bestj], allpopnames[i]);
		double fitted2 = gsl_matrix_get(sigma_cor, bestj, i);
		double diff2 = cov2-fitted2;
		if (diff2 >= max_test) pops2test.insert(make_pair( allpopnames[bestj], allpopnames[i]));
	}

	//2. Get the paths to the root for both
	map<string, Graph::vertex_descriptor> p2node = tree->get_tips(tree->root);


	//3. try migration to all the pairwise combinations

	double max_llik = current_llik;

	//tree_bk holds a backup of the current tree; bk_2 the best tree
	tree_bk->copy(tree);
	tree_bk2->copy(tree);
	pair<int, int> best_edge;
	set<Graph::vertex_descriptor> root_adj = tree->get_root_adj();
	for (set<pair<string, string> >::iterator pit = pops2test.begin(); pit != pops2test.end(); pit++){
		cout << "Trying "<< pit->first<< " "<< pit->second << "\n"; cout.flush();
		set<Graph::vertex_descriptor> p1_s = tree->get_path_to_root(p2node[pit->first]);
		set<Graph::vertex_descriptor> p2_s = tree->get_path_to_root(p2node[pit->second]);
		for (set<Graph::vertex_descriptor>::iterator it = p1_s.begin(); it != p1_s.end(); it++){
			for (set<Graph::vertex_descriptor>::iterator it2 = p2_s.begin(); it2 != p2_s.end(); it2++){
				if (!try_mig(*it, *it2, tmpfitted)) continue;

				cout << tree->get_newick_format(*it)<< " "<< tree->get_newick_format(*it2)<< "\n";
				pair<int, int> totest = make_pair( tree->g[*it].index, tree->g[*it2].index);
				if ( tested.find(totest) == tested.end() && tree->is_legal_migration(*it, *it2)){

					Graph::edge_descriptor e = tree->add_mig_edge( *it, *it2);
					optimize_weight(e);
					cout << "1->2 "<< tree->g[e].weight << " "<< current_llik << " "<< max_llik << "\n";

					if (current_llik > max_llik){
						tree_bk2->copy(tree);
						gsl_matrix_memcpy( tmpfitted2, sigma_cor);
						max_llik = current_llik;
						best_edge.first = tree->g[*it].index;
						best_edge.second = tree->g[*it2].index;
						toreturn.first = true;
					}

					tree->remove_mig_edge(e);


					// test to see if one of these is around the root
					if (root_adj.find(*it) != root_adj.end()){
						set<Graph::vertex_descriptor>::iterator rit = root_adj.begin();
						if (tree->g[*rit].index == tree->g[*it].index) rit++;
							if (tree->is_legal_migration(*rit, *it2)){
								Graph::edge_descriptor e = tree->add_mig_edge( *rit, *it2);
								optimize_weight(e);
								cout << "1 (root) ->2 "<< tree->g[e].weight << " "<< current_llik << " "<< max_llik << "\n";
								if (current_llik > max_llik){
									tree_bk2->copy(tree);
									gsl_matrix_memcpy( tmpfitted2, sigma_cor);
									max_llik = current_llik;
									best_edge.first = tree->g[*rit].index;
									best_edge.second = tree->g[*it2].index;
									toreturn.first = true;
								}

								tree->remove_mig_edge(e);
							}
					}
					if (root_adj.find(*it2) != root_adj.end()){
						set<Graph::vertex_descriptor>::iterator rit = root_adj.begin();
						if (tree->g[*rit].index == tree->g[*it2].index) rit++;
							if (tree->is_legal_migration(*it, *rit)){
								Graph::edge_descriptor e = tree->add_mig_edge( *it, *rit);
								optimize_weight(e);
								cout << "1->2 (root) "<< tree->g[e].weight << " "<< current_llik << " "<< max_llik << "\n";
								if (current_llik > max_llik){
									tree_bk2->copy(tree);
									gsl_matrix_memcpy( tmpfitted2, sigma_cor);
									max_llik = current_llik;
									best_edge.first = tree->g[*it].index;
									best_edge.second = tree->g[*rit].index;
									toreturn.first = true;
								}

								tree->remove_mig_edge(e);
							}
					}
					tested.insert(make_pair( tree->g[*it].index, tree->g[*it2].index));
				}

				totest = make_pair( tree->g[*it2].index, tree->g[*it].index);
				if ( tested.find(totest) == tested.end() && tree->is_legal_migration(*it2, *it)){

					Graph::edge_descriptor e = tree->add_mig_edge( *it2, *it);
					optimize_weight(e);

					cout << "2->1 "<< tree->g[e].weight << " "<< current_llik << " "<< max_llik << "\n";
					if (current_llik > max_llik){
						tree_bk2->copy(tree);
						gsl_matrix_memcpy( tmpfitted2, sigma_cor);
						max_llik = current_llik;
						best_edge.first = tree->g[*it2].index;
						best_edge.second = tree->g[*it].index;
						toreturn.first = true;
					}

					tree->remove_mig_edge(e);

					if (root_adj.find(*it) != root_adj.end()){

						set<Graph::vertex_descriptor>::iterator rit = root_adj.begin();
						if (tree->g[*rit].index == tree->g[*it].index) rit++;
							if (tree->is_legal_migration(*it2, *rit)){

								Graph::edge_descriptor e = tree->add_mig_edge( *it2, *rit);

								optimize_weight(e);
								cout << "2->1 (root) "<< tree->g[e].weight << " "<< current_llik << " "<< max_llik << "\n";
								if (current_llik > max_llik){
									tree_bk2->copy(tree);
									gsl_matrix_memcpy( tmpfitted2, sigma_cor);
									max_llik = current_llik;
									best_edge.first = tree->g[*it2].index;
									best_edge.second = tree->g[*rit].index;
									toreturn.first = true;
								}

								tree->remove_mig_edge(e);
							}

					}
					if (root_adj.find(*it2) != root_adj.end()){

						set<Graph::vertex_descriptor>::iterator rit = root_adj.begin();
						if (tree->g[*rit].index == tree->g[*it2].index) rit++;
							if (tree->is_legal_migration(*rit, *it)){

								Graph::edge_descriptor e = tree->add_mig_edge( *rit, *it);

								optimize_weight(e);
								cout << "2 (root) ->1 "<< tree->g[e].weight << " "<< current_llik << " "<< max_llik << "\n";
								if (current_llik > max_llik){
									tree_bk2->copy(tree);
									gsl_matrix_memcpy( tmpfitted2, sigma_cor);
									max_llik = current_llik;
									best_edge.first = tree->g[*rit].index;
									best_edge.second = tree->g[*it].index;
									toreturn.first = true;
								}


								tree->remove_mig_edge(e);
							}

					}

					tested.insert(make_pair( tree->g[*it2].index, tree->g[*it].index));
				}

				Graph::vertex_descriptor p1 = tree->get_parent_node(*it).first;
				Graph::vertex_descriptor p2 = tree->get_parent_node(*it2).first;

				totest = make_pair( tree->g[p1].index, tree->g[*it2].index);
				if (tested.find(totest) == tested.end() && tree->is_legal_migration(p1, *it2)){
					//cout << "trying "<< tree->get_newick_format(p1)<< " "<< tree->get_newick_format(*it2)<< " rev\n";
					Graph::edge_descriptor e = tree->add_mig_edge( p1, *it2);

					optimize_weight(e);
					//cout << "after "<< tree->get_newick_format(*it)<< " "<< tree->get_newick_format(*it2)<< "\n";
					cout << "p1->2 "<< tree->g[e].weight << " "<< current_llik << " "<< max_llik << "\n";
					if (current_llik > max_llik){
						tree_bk2->copy(tree);
						gsl_matrix_memcpy( tmpfitted2, sigma_cor);
						max_llik = current_llik;
						best_edge.first = tree->g[p1].index;
						best_edge.second = tree->g[*it2].index;
						toreturn.first = true;
					}

					tree->remove_mig_edge(e);
					tested.insert(make_pair( tree->g[p1].index, tree->g[*it2].index));
				}

				totest = make_pair( tree->g[p2].index, tree->g[*it].index);
				if ( tested.find(totest) == tested.end() && tree->is_legal_migration(p2, *it)){
					//cout << "trying "<< tree->get_newick_format(p2)<< " "<< tree->get_newick_format(*it)<< " rev\n";
					Graph::edge_descriptor e = tree->add_mig_edge( p2, *it);

					optimize_weight(e);
					//cout << "after "<< tree->get_newick_format(*it)<< " "<< tree->get_newick_format(*it2)<< "\n";
					cout << "p2->1 "<< tree->g[e].weight << " "<< current_llik << " "<< max_llik << "\n";
					if (current_llik > max_llik){
						tree_bk2->copy(tree);
						gsl_matrix_memcpy( tmpfitted2, sigma_cor);
						max_llik = current_llik;
						best_edge.first = tree->g[p2].index;
						best_edge.second = tree->g[*it].index;
						toreturn.first = true;
					}

					tree->remove_mig_edge(e);
					tested.insert(make_pair( tree->g[p2].index, tree->g[*it].index));
				}
			}
		}
	}
	//cout << "here1!\n"; cout.flush();
	if (toreturn.first == true)	{
		//cout << "doing something\n"; cout.flush();
		tree->copy(tree_bk2);
		gsl_matrix_memcpy( sigma_cor, tmpfitted2);
		current_llik = max_llik;
		toreturn.second.first = best_edge.first;
		toreturn.second.second = best_edge.second;

	}
	else {
		gsl_matrix_memcpy( sigma_cor, tmpfitted);
		current_llik = llik_bk;
		tree->copy(tree_bk);
	}

	//cout << "here!\n"; cout.flush();
	////tree->print();
	gsl_matrix_free(tmpfitted);
	gsl_matrix_free(tmpfitted2);
	return toreturn;

}


pair< pair<bool, bool>, pair<double, pair<int, int> > > GraphState2::add_mig_targeted(string p1, string p2){
	// find the largest residual, try migration events in the vicinity
	// return true if an event is added, false otw

	//2. Get the paths to the root for both
	pair< pair<bool, bool>, pair<double, pair<int, int> > > toreturn;
	toreturn.first.first = false;
	toreturn.first.second = false;
	map<string, Graph::vertex_descriptor> p2node = tree->get_tips(tree->root);

	vector<Graph::vertex_descriptor> p1_s = tree->get_path_to_root_vec(p2node[p1]);
	vector<Graph::vertex_descriptor> p2_s = tree->get_path_to_root_vec(p2node[p2]);

	//3. try migration to all the pairwise combinations
	double max_llik = current_llik;
	tree_bk->copy(tree);
	double migw = 0;
	pair<int, int> best_edge;
	for(int i = 0; i < p1_s.size(); i++){
		for (int j = 0; j < p2_s.size(); j++){

			if ( tree->is_legal_migration(p1_s[i], p2_s[j])){
				Graph::edge_descriptor e = tree->add_mig_edge( p1_s[i] , p2_s[j]);

				optimize_weights();

				if (current_llik > max_llik){
					tree_bk->copy(tree);
					max_llik = current_llik;
					best_edge.first = i;
					best_edge.second = j;
					toreturn.first.first = true;
					migw = tree->g[e].weight;
				}

				tree->remove_mig_edge(e);
			}
			if ( tree->is_legal_migration(p2_s[j], p1_s[i])){
				Graph::edge_descriptor e = tree->add_mig_edge( p2_s[j], p1_s[i]);

				optimize_weights();
				if (current_llik > max_llik){
					tree_bk->copy(tree);
					max_llik = current_llik;
					best_edge.first = i;
					best_edge.second = j;
					toreturn.first.first = true;
					toreturn.first.second = true;
					migw = tree->g[e].weight;
				}
				tree->remove_mig_edge(e);
			}
		}
	}
	if (toreturn.first.first == true)	{

		tree->copy(tree_bk);
		set_branches_ls_wmig();
		toreturn.second.first = migw;
		toreturn.second.second = best_edge;

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
		//cout << pop  << " "<< t << " "<< meanhzy << " "<< mean_n << "\n";
		trim.insert(make_pair(pop, t));
	}
	toreturn = tree->get_newick_format(&trim);
	return toreturn;
}


void GraphState2::print_trimmed(string stem){
	map<string, double> trim;
	string toreturn;

	for ( map<string, int>::iterator it = countdata->pop2id.begin(); it!= countdata->pop2id.end(); it++){
		int id = it->second;
		string pop = it->first;
		double meanhzy = countdata->mean_hzy.find(id)->second;
		double mean_n = countdata->mean_ninds.find(id)->second;
		double t = meanhzy / (4.0* mean_n);
		//cout << pop  << " "<< t << " "<< meanhzy << " "<< mean_n << "\n";
		trim.insert(make_pair(pop, t));
	}
	tree->print(stem, &trim);
}

Graph::vertex_descriptor GraphState2::get_neighborhood(Graph::vertex_descriptor v){
	Graph::vertex_descriptor toreturn = v;
	int i = 0;
	while (i < params->m_neigh && tree->g[ tree->get_parent_node(toreturn).first ].is_root == false){
		toreturn = tree->get_parent_node(toreturn).first;
		i++;
	}
	return toreturn;
}


bool GraphState2::try_mig(Graph::vertex_descriptor v1, Graph::vertex_descriptor v2, gsl_matrix * fitted){
	map<string, Graph::vertex_descriptor> tips1 = tree->get_tips(v1);
	map<string, Graph::vertex_descriptor> tips2 = tree->get_tips(v2);
	for (map<string, Graph::vertex_descriptor>::iterator it1 = tips1.begin(); it1 != tips1.end(); it1++){
		string p1 = it1->first;
		int index1 = popname2index[p1];
		for (map<string, Graph::vertex_descriptor>::iterator it2 = tips2.begin(); it2 != tips2.end(); it2++){
			string p2 = it2->first;
			int index2 = popname2index[p2];
			if (p1==p2) continue;
			double resid = countdata->get_cov(p1, p2) -  gsl_matrix_get(fitted, index1, index2);
			if ( resid < 0) return false;
		}
	}
	return true;
}

void GraphState2::place_root(string r){
	tree->place_root(r);
	cout << "Set root above "<< r << "\n"; cout.flush();
	set_branches_ls();
	current_llik = llik();
	cout << tree->get_newick_format()<< "\n"; cout.flush();
	cout << "ln(lk) = "<< current_llik << "\n";
}

void GraphState2::flip_mig(){
	vector<Graph::edge_descriptor> m = tree->get_mig_edges();
	//tree->print("test2");
	int i = 0;
	while (i < m.size()){
		//cout << i << " in flip\n";
		Graph::vertex_descriptor v = target(m[i], tree->g);
		//cout << tree->g[v].index << " "<< tree->g[source(*it, tree->g)].index << "\n";
		set<Graph::edge_descriptor> inm = tree->get_in_mig_edges(v);
		double w = 0;
		double max = 0;
		Graph::edge_descriptor e;
		for (set<Graph::edge_descriptor>::iterator it3 = inm.begin(); it3!= inm.end(); it3++ ){
			w += tree->g[*it3].weight;
			if (tree->g[*it3].weight > max){
				max = tree->g[*it3].weight;
				e = *it3;
			}
		}
		if (w > 0.6 && tree->g[e].is_mig){
			Graph::vertex_descriptor t = target(e, tree->g);
			Graph::vertex_descriptor s = source(e, tree->g);
			//cout << "flipping "<< tree->g[t].index << " "<< tree->g[s].index << "\n"; cout.flush();
			if ( tree->g[tree->get_parent_node(t).first].is_root) {
				i++;
				continue;
			}
			if ( tree->g[ tree->get_parent_node(t).first].index ==  tree->g[ tree->get_parent_node(s).first].index ) {
				i++;
				continue;
			}
			double totalw = 0;
			set<Graph::edge_descriptor> inm = tree->get_in_mig_edges(t);
			for( set<Graph::edge_descriptor>::iterator it2 = inm.begin(); it2 != inm.end(); it2++) totalw += tree->g[*it2].weight;
			double left = 1-totalw;
			Graph::edge_descriptor inc;
			pair<Graph::in_edge_iterator, Graph::in_edge_iterator> ine = in_edges(t, tree->g);
			while (ine.first != ine.second){
				if (tree->g[*ine.first].is_mig == false) inc = *ine.first;
					ine.first++;
			}

			Graph::vertex_descriptor newm = source(inc, tree->g);
			if (tree->g[newm].is_mig) {
				i++;
				continue;
			}
			set<Graph::edge_descriptor> inm2 = tree->get_in_mig_edges(newm);

			if (inm2.size() > 0 ) {
				i++;
				continue;
			}
			tree->g[newm].is_mig = true;
			tree->g[newm].mig_frac = 0.5;

			tree->g[s].is_mig = false;
			tree->g[s].mig_frac = 0;

			tree->g[inc].is_mig = true;
			tree->g[inc].weight = left;
			tree->g[inc].len = 0;
			//cout << "here?\n"; cout.flush();
			tree->set_mig_frac(inc, 0.5);


			tree->g[e].is_mig = false;
			tree->g[e].len = 1;
			set_branches_ls_wmig();
			optimize_frac(inc);
			optimize_weight(inc);
			current_llik = llik();
			//iterate_movemig( tree->g[newm].index);
			//iterate_local_hillclimb_wmig(tree->g[t].index);
			m = tree->get_mig_edges();
			i = 0;
			continue;
		}
		i++;
	}

}

void GraphState2::trim_mig(){

	vector<Graph::edge_descriptor> m = tree->get_mig_edges();
	int i = 0;
	while (i < m.size()){
			//cout << i <<" in trim\n";
			double w = tree->g[ m[i] ].weight;
			Graph::vertex_descriptor v = target( m[i], tree->g);
			if (w < params->min_migw){
				tree->remove_mig_edge(m[i]);
				set_branches_ls_wmig();
				m = tree->get_mig_edges();
				i = 0;
				continue;
			}
			i++;
	}
}



bool GraphState2::has_loop(){
	bool toreturn = false;
	graph_traits<Graph>::edge_iterator ei, ei_end;
	for (tie(ei, ei_end) = edges(tree->g); ei != ei_end; ++ei){
		if (!tree->g[*ei].is_mig ) continue;
		int si = tree->g[source(*ei, tree->g)].index;
		int ti = tree->g[target(*ei, tree->g)].index;
		tree_bk3->copy(tree);
		Graph::edge_descriptor e;
		graph_traits<Graph>::edge_iterator ei2, ei_end2;
		Graph::vertex_descriptor t1;
		Graph::vertex_descriptor t2;
		for (tie(ei2, ei_end2) = edges(tree_bk3->g); ei2 != ei_end2; ++ei2){
			if (tree_bk3->g[source(*ei2, tree_bk3->g)].index == si && tree_bk3->g[target(*ei2, tree_bk3->g)].index == ti) {
				e = *ei2;
				t1 = source(*ei2, tree_bk3->g);
				t1 = tree_bk3->get_child_node_mig(t1);
				t2 = target(*ei2, tree_bk3->g);
			}
		}
		tree_bk3->remove_mig_edge(e);
		if (!tree_bk3->is_legal_migration(t1, t2)) return true;
	}
	return toreturn;
}

void GraphState2::iterate_movemig(int index){

	pair<bool, int> moving = movemig(index);
	while( moving.first){
		moving = movemig(moving.second);
	}

}

pair<bool, int> GraphState2::movemig( int index ){

	pair<bool, int> toreturn = make_pair(false, 0);
	double max = current_llik;
	double lik_bk = current_llik;
	tree_bk->copy(tree);
	tree_bk2->copy(tree);
	gsl_matrix *tmpfitted = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted, sigma_cor);
	gsl_matrix *tmpfitted2 = gsl_matrix_alloc(current_npops, current_npops);
	gsl_matrix_memcpy( tmpfitted2, sigma_cor);
	tree->print("test");
	for (int i = 0; i < 3 ; i++){
		cout << "moving "<< i <<"\n";
		map<int, Graph::vertex_descriptor> vindex = tree->index2vertex();
		Graph::vertex_descriptor s = vindex[index];
		Graph::vertex_descriptor p = tree->get_parent_node(s).first;
		Graph::vertex_descriptor ch = tree->get_child_node_mig(s);
		Graph::edge_descriptor e = tree->get_out_mig_edge(s);
		Graph::vertex_descriptor t = target(e, tree->g);
		pair<Graph::vertex_descriptor, Graph::vertex_descriptor> ch2 = tree->get_child_nodes(ch);

		if ( i == 0){
			if (tree->g[p].is_root){
				pair<Graph::vertex_descriptor, Graph::vertex_descriptor> ch3 = tree->get_child_nodes(p);
				if ( tree->g[ch3.first].index == tree->g[ch].index )	 p = ch3.second;
				else p = ch3.first;
			}

			tree->remove_mig_edge(e);
			if (tree->is_legal_migration(p, t)){
				Graph::edge_descriptor e2 = tree->add_mig_edge(p, t);
				optimize_weight(e2);
				tree->print("test0");
				cout << i << " "<< current_llik << " "<< max << "\n";
				if ( current_llik > max){
					max = current_llik;
					tree_bk2->copy(tree);
					gsl_matrix_memcpy( tmpfitted2, sigma_cor);
					toreturn.first = true;
					toreturn.second = tree->g[source(e2, tree->g)].index;

				}
			}
			tree->copy(tree_bk);
		}
		else if (i == 1){
			if (tree->g[ch].is_tip) continue;
			tree->remove_mig_edge(e);
			if (tree->is_legal_migration(ch2.first, t)){
				Graph::edge_descriptor e2 = tree->add_mig_edge(ch2.first, t);
				optimize_weight(e2);
				tree->print("test1");
				cout << i << " "<< current_llik << " "<< max << "\n";
				if ( current_llik > max){
					max = current_llik;
					tree_bk2->copy(tree);
					gsl_matrix_memcpy( tmpfitted2, sigma_cor);
					toreturn.first = true;
					toreturn.second = tree->g[source(e2, tree->g)].index;

				}
			}
			tree->copy(tree_bk);
		}
		else if (i == 2){
			if (tree->g[ch].is_tip) continue;
			tree->remove_mig_edge(e);
			if (tree->is_legal_migration(ch2.second, t)){
				Graph::edge_descriptor e2 = tree->add_mig_edge(ch2.second, t);
				optimize_weight(e2);
				tree->print("test2");
				cout << i << " "<< current_llik << " "<< max << "\n";
				if ( current_llik > max){
					max = current_llik;
					tree_bk2->copy(tree);
					gsl_matrix_memcpy( tmpfitted2, sigma_cor);
					toreturn.first = true;
					toreturn.second = tree->g[source(e2, tree->g)].index;

				}
			}
			tree->copy(tree_bk);
		}


	}
	if (toreturn.first){
		tree->copy(tree_bk2);
		gsl_matrix_memcpy( sigma_cor, tmpfitted2);
		current_llik = max;
	}
	else{
		gsl_matrix_memcpy( sigma_cor, tmpfitted);
		current_llik = lik_bk;
		tree->copy(tree_bk);
	}
	gsl_matrix_free(tmpfitted);
	gsl_matrix_free(tmpfitted2);
	return toreturn;
}



