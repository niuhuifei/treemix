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
	PopGraph();
	PopGraph(string);
	Graph g;
	vector<Graph::vertex_descriptor> index2father;
};
#endif /* POPGRAPH_H_ */
