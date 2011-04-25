/*
 * PopGraph.h
 *
 *  Created on: Apr 21, 2011
 *      Author: pickrell
 */

#ifndef POPGRAPH_H_
#define POPGRAPH_H_
#include "Settings.hpp"


struct Node
{
	int index;
	string name;
	double height;
	bool is_tip;
	bool is_root;
	bool rev;
};

struct Dist
{
	float weight;
	float len;
};

typedef adjacency_list<vecS, listS, boost::bidirectionalS, Node, Dist> Graph;
typedef pair<int, int> Edge;
typedef boost::property_map<Graph, int Node::*>::type IndexMap;
typedef boost::graph_traits<Graph>::vertex_iterator vertex_iter;

class PopGraph{
public:

	// initialize the graph
	PopGraph();

	// initialize from a Newick string
	PopGraph(string);
	Graph g;
	bool istree; // is this a tree? if so, allow shortcuts

	vector<Graph::vertex_descriptor> index2father;
	map<string, Graph::vertex_descriptor> popname2tip;
	Graph::vertex_descriptor root;

	inline double get_dist_to_root(Graph::vertex_descriptor);
	inline void flip_sons(Graph::vertex_descriptor, gsl_rng* );

	//get the inorder traversal
	inline vector< Graph::vertex_descriptor > get_inorder_traversal( int);

	//helper functions for inorder traversal
	inline void inorder_traverse(Graph::vertex_descriptor, int*, vector<Graph::vertex_descriptor >*);
};
#endif /* POPGRAPH_H_ */
