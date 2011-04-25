/*
 * PopGraph.cpp
 *
 *  Created on: Apr 21, 2011
 *      Author: pickrell
 */
#include "PopGraph.h"

PopGraph::PopGraph(){
	enum {A, B, C, D, E, N};
	const int num_vertices = N;
	const char * name = "ABCDE";
	Edge edge_array[] =
	{ Edge(A,B), Edge(A,D), Edge(C,A), Edge(D,C),
	      Edge(C,E), Edge(B,D), Edge(D,E) };
	const int num_edges = sizeof(edge_array)/sizeof(edge_array[0]);
	g = Graph(num_vertices);
	//for (int i = 0; i < num_edges; ++i) add_edge(edge_array[i].first, edge_array[i].second, g);
}

PopGraph::PopGraph(string p_newickString){
	g = Graph(1);
	index2father.clear();
	int i = 0;
	istree = true;
	Graph::vertex_descriptor v = *vertices(g).first;
	Graph::vertex_descriptor v2 = *vertices(g).first;
	Graph::edge_descriptor e;
	g[v].index = i;
	g[v].name = "NA";
	g[v].height = 0;
	index2father.push_back(NULL);
	g[v].is_tip = false;
	g[v].is_root = true;
	g[v].rev = false;
	root = v;
	for(string::const_iterator I = p_newickString.begin();
		I != p_newickString.end(); ++I)
      {
		//cout << *I << "\n";
		if ( *I == '(' )//Begin new subtree
		{
			v2 = add_vertex(g);
			i++;
			g[v2].index = i;
			g[v2].name = "NA";
			g[v2].height = 0;
			index2father.push_back(v);
			g[v2].is_tip = false;
			g[v2].is_root = false;
			g[v2].rev = false;
			e = add_edge( v, v2, g).first;
			v = v2;
		}
		else if( *I == ')' )// Subtree finished, get back
		{

			v = index2father[g[v].index];
			//e = add_edge( index2father[g[v].index], v, g).first;

		}
		else if( *I == ',' )// insert brother
		{
			v = index2father[g[v].index];
			v2 = add_vertex(g);
			i++;
			g[v2].index = i;
			g[v2].name = "NA";
			g[v2].height = 0;
			index2father.push_back(v);
			g[v2].is_tip = false;
			g[v2].is_root = false;
			g[v2].rev = false;
			e = add_edge( v, v2, g).first;
			v = v2;
      }
      else if( *I == ':' )// treelength
      {
              std::string length = "";
              ++I;
              while( *I!=',' && *I!=')' && *I!=':' && *I!='(' && *I!=';')
              {
                      length += *I;
                      ++I;
              }
              --I;
              graph_traits<Graph>::in_edge_iterator in_i;
              graph_traits<Graph>::in_edge_iterator in_end;
              for ( tie(in_i, in_end)= in_edges(v, g); in_i != in_end; ++in_i ){
            	  g[*in_i].len = atof(length.c_str());
            	  g[*in_i].weight = 1;
              }
      }
      else if( *I == ';' )
      {
              break;
      }
      else// name
      {
              std::string name = "";
              do
              {
                      name += *I;
                      ++I;
              }
              while( *I!=',' && *I!=')' && *I!=':' && *I!='(' && *I!=';');
              --I;
              g[v2].name = name;
              g[v2].is_tip = true;
              popname2tip.insert(make_pair(name, v2));
      }
     }
}

double PopGraph::get_dist_to_root(Graph::vertex_descriptor v){
	double toreturn = 0;
	while (g[v].is_root == false){
		graph_traits<Graph>::in_edge_iterator in_i = in_edges(v, g).first;
		toreturn += g[*in_i].len;
		v = source(*in_i, g);
	}
	return toreturn;
}

void PopGraph::flip_sons(Graph::vertex_descriptor v, gsl_rng* r){
	if (g[v].is_tip == false){
		//cout << "here\n"; cout.flush();
		double ran = gsl_rng_uniform(r);
		//cout << ran <<"\n"; cout.flush();
		if (ran <0.5){
			if (g[v].rev == true) g[v].rev = false;
			else g[v].rev = true;
		}
		graph_traits<Graph>::out_edge_iterator out_i = out_edges(v, g).first;
		//flip_sons(dest(g[*out_i]), r);
		++out_i;
		//flip_sons(destination(g[*out_i]), r);
	}
}

vector<Graph::vertex_descriptor > PopGraph::get_inorder_traversal(int nodes){
        	vector<Graph::vertex_descriptor> toreturn(2*nodes-1);
        	int count = 0;
        	inorder_traverse(root, &count, &toreturn);
        	return toreturn;
}

void PopGraph::inorder_traverse(Graph::vertex_descriptor rootIterator, int* i, vector<Graph::vertex_descriptor >* v){
	graph_traits<Graph>::out_edge_iterator out_i = out_edges(v, g).first;
	if (g[rootIterator].is_tip == false){
		if (g[v].rev == true) out_i++;
 		//inorder_traverse(*out_i,i , v);
 	}
 	v->at(*i) = rootIterator;
 	*i = *i+1;

 	if (g[rootIterator].is_tip == false){
 		if (g[v].rev ==true) out_i--;
 		else ++out_i;
 		//inorder_traverse(*out_i, i, v);
 	}
}
