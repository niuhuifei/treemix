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
	float mig_frac;
	bool is_tip;
	bool is_root;
	bool is_mig;
	bool rev;
};

struct Dist
{
	float weight;
	float len;
	bool is_mig;
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
	void set_graph(string); //or set it to a Newick string

	// initialize with three populations
	PopGraph(vector<string>);

	//copy from another location
	void copy(PopGraph*);
	//void copy_helper(Graph::vertex_descriptor, Graph::vertex_descriptor, PopGraph *);

	Graph g;
	bool istree; // is this a tree? if so, allow shortcuts
	int indexcounter;

	vector<Graph::vertex_descriptor> index2father;
	vector<string> popnames;
	Graph::vertex_descriptor root;


	Graph::vertex_descriptor add_tip(Graph::vertex_descriptor, string);
	void remove_tip(Graph::vertex_descriptor);


	void set_root(Graph::vertex_descriptor);
	set<Graph::vertex_descriptor> get_root_adj();
	void print();
	map<string, Graph::vertex_descriptor> get_tips( Graph::vertex_descriptor);

	double get_height(Graph::vertex_descriptor);
	double get_dist_to_root(Graph::vertex_descriptor); //get the distance to the root for any node in the tree
	void flip_sons(Graph::vertex_descriptor, gsl_rng* ); //flip the inorder traversal from a node

	//get the inorder traversal
	vector< Graph::vertex_descriptor > get_inorder_traversal( int);
	void inorder_traverse(Graph::vertex_descriptor, int*, vector<Graph::vertex_descriptor >*);
	vector< Graph::vertex_descriptor > get_inorder_traversal_noroot( int);

	//get the LCA of two tips in the tree
	Graph::vertex_descriptor get_LCA(Graph::vertex_descriptor, Graph::vertex_descriptor, Graph::vertex_descriptor);

	//get the path to the root for a vertex
	set<Graph::vertex_descriptor> get_path_to_root(Graph::vertex_descriptor);

	//local rearrangement
	void local_rearrange(Graph::vertex_descriptor, int);
	void move_root(gsl_rng*);

	//global rearrangement
	void global_rearrange(Graph::vertex_descriptor, Graph::vertex_descriptor);

	//Newick format
	string get_newick_format();
	void newick_helper( Graph::vertex_descriptor, string*);
	string get_newick_format(map<string, double>*);
	void newick_helper( Graph::vertex_descriptor, string*, map<string, double>*);

	//MIGRATION GRAPH
	Graph::vertex_descriptor add_mig_edge(Graph::vertex_descriptor, Graph::vertex_descriptor);
	set<pair<double, set<Graph::vertex_descriptor> > > get_paths_to_root(Graph::vertex_descriptor);
	set<pair<double, set<Graph::edge_descriptor> > > get_paths_to_root_edge(Graph::vertex_descriptor);

	pair<Graph::vertex_descriptor, Graph::vertex_descriptor> get_child_nodes(Graph::vertex_descriptor); // skips over migration nodes
	Graph::vertex_descriptor get_child_node_mig(Graph::vertex_descriptor); //input a migration node, get the next non-migration node (along non-migration edges)
	pair<Graph::vertex_descriptor, double> get_parent_node(Graph::vertex_descriptor); //skips over migration nodes
	pair<Graph::vertex_descriptor, double> get_parent_node_wmig(Graph::vertex_descriptor); // skips migration edges, but will return a migration node. returns the node and distance
	bool does_mig_exist(Graph::vertex_descriptor, Graph::vertex_descriptor); //does a migration edge already exist between these two nodes?
	bool is_legal_migration(Graph::vertex_descriptor, Graph::vertex_descriptor); //is migration between these two nodes legal? (no migration already, not the same parent, not creating loop)
};

#endif /* POPGRAPH_H_ */