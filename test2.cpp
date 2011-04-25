/*
 * test2.cpp
 *
 *  Created on: Apr 7, 2011
 *      Author: pickrell
 */

#include "PopGraph.h"

int
main (void)
{

	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_ranlxs2;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, 100);


	string testnewick = "((YRI:0.8,LWD:0.4):0.1,(JPT:0.2,CHB:0.3):0.5);";
	PopGraph tmpgraph(testnewick);

	IndexMap index = get(&Node::index, tmpgraph.g);


	cout << "vertices(g) = ";
	pair<vertex_iter, vertex_iter> vp;
	for (vp = vertices(tmpgraph.g); vp.first != vp.second; ++vp.first)
		std::cout << index[*vp.first] <<  ":"<< tmpgraph.g[*vp.first].name << ":"<< tmpgraph.get_dist_to_root(*vp.first)<< " ";
	cout << "\n";

    std::cout << "edges(g) = ";
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(tmpgraph.g); ei != ei_end; ++ei)
        std::cout << "(" << index[source(*ei, tmpgraph.g)]
                  << "," << index[target(*ei, tmpgraph.g)] << ","<< tmpgraph.g[*ei].len <<  ") ";
    std::cout << std::endl;

    vector<Graph::vertex_descriptor> tmpinorder = tmpgraph.get_inorder_traversal(4);
    for (int i = 0; i < tmpinorder.size(); i++){
    	cout << i << " "<< tmpgraph.g[tmpinorder[i]].name << "\n";
    }
    cout << "\n";

    Graph::vertex_descriptor tv = tmpgraph.get_LCA(tmpgraph.root, tmpgraph.popname2tip["YRI"], tmpgraph.popname2tip["JPT"]);
    cout << tmpgraph.g[tv].index << "\n";
    cout<< tmpgraph.get_newick_format()<< "\n";
    tmpgraph.randomize_tree(r);
    cout<< tmpgraph.get_newick_format()<< "\n";
    return 0;
}
