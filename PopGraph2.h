/*
 * PopGraph2.h
 *
 *  Created on: Jun 17, 2011
 *      Author: pickrell
 */

#ifndef POPGRAPH2_H_
#define POPGRAPH2_H_

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

class PopGraph2{
	public:
	// initialize the graph
	PopGraph2(string, vector<string>);
	Graph::vertex_descriptor add_tip(Graph::vertex_descriptor, string);
	void remove_tip(Graph::vertex_descriptor);


	//copy from another location
	void copy(PopGraph2*);

	//print
	void print();

	double get_dist_to_root(Graph::vertex_descriptor); //get the distance to the root for any node in the tree

	//inorder traversal functions
	vector< Graph::vertex_descriptor > get_inorder_traversal( int);
	void inorder_traverse(Graph::vertex_descriptor, int*, vector<Graph::vertex_descriptor >*);
	vector< Graph::vertex_descriptor > get_inorder_traversal_noroot( int);

	//get the LCA of two tips in the tree
	Graph::vertex_descriptor get_LCA(Graph::vertex_descriptor, Graph::vertex_descriptor, Graph::vertex_descriptor);

	//get the path to the root for a vertex
	set<Graph::vertex_descriptor> get_path_to_root(Graph::vertex_descriptor);

	Graph g;
	bool istree; // is this a tree? if so, allow shortcuts
	vector<string> popnames;
	Graph::vertex_descriptor root;
	int indexcounter;

};
#endif /* POPGRAPH2_H_ */
