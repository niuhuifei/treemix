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
	string testnewick = "((YRI:0.8,LWD:0.4):0.1,(JPT:0.2,CHB:0.3):0.5);";
	PopGraph tmpgraph(testnewick);

	IndexMap index = get(&Node::index, tmpgraph.g);


	cout << "vertices(g) = ";
	pair<vertex_iter, vertex_iter> vp;
	for (vp = vertices(tmpgraph.g); vp.first != vp.second; ++vp.first)
		std::cout << index[*vp.first] <<  ":"<< tmpgraph.g[*vp.first].name << " ";
	cout << "\n";

    std::cout << "edges(g) = ";
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(tmpgraph.g); ei != ei_end; ++ei)
        std::cout << "(" << index[source(*ei, tmpgraph.g)]
                  << "," << index[target(*ei, tmpgraph.g)] << ","<< tmpgraph.g[*ei].len <<  ") ";
    std::cout << std::endl;
 // std::cout << "edges(g) = ";
 // boost::graph_traits<Graph>::edge_iterator ei, ei_end;
 // for (boost::tie(ei, ei_end) = edges(tmpgraph.g); ei != ei_end; ++ei)
//	  std::cout << "(" << index[source(*ei, tmpgraph.g)]
//	                            << "," << index[target(*ei, tmpgraph.g)] << ") ";
 // std::cout << std::endl;
  return 0;
}
