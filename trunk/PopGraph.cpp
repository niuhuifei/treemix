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
              popnames.push_back(name);
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
		double ran = gsl_rng_uniform(r);
		if (ran <0.5){
			if (g[v].rev == true) g[v].rev = false;
			else if (g[v].rev == false) g[v].rev = true;
		}
		graph_traits<Graph>::out_edge_iterator out_i = out_edges(v, g).first;
		flip_sons(target(*out_i, g), r);
		++out_i;
		flip_sons(target(*out_i, g), r);
	}
}

vector<Graph::vertex_descriptor > PopGraph::get_inorder_traversal(int nodes){
        	vector<Graph::vertex_descriptor> toreturn(2*nodes-1);
        	int count = 0;
        	inorder_traverse(root, &count, &toreturn);
        	return toreturn;
}

void PopGraph::inorder_traverse(Graph::vertex_descriptor rootIterator, int* i, vector<Graph::vertex_descriptor >* v){
	graph_traits<Graph>::out_edge_iterator out_i = out_edges(rootIterator, g).first;
	if (g[rootIterator].is_tip == false){
		if (g[rootIterator].rev == true) out_i++;
 		inorder_traverse(target(*out_i, g), i , v);
 	}
 	v->at(*i) = rootIterator;
 	*i = *i+1;

 	if (g[rootIterator].is_tip == false){
 		if (g[rootIterator].rev ==true) out_i--;
 		else ++out_i;
 		inorder_traverse(target(*out_i, g), i, v);
 	}
}

Graph::vertex_descriptor PopGraph::get_LCA(Graph::vertex_descriptor root_v,
		Graph::vertex_descriptor tip1_v, Graph::vertex_descriptor tip2_v){
		if (g[root_v].is_tip == true) return NULL;
		graph_traits<Graph>::out_edge_iterator out_i = out_edges(root_v, g).first;
		bool found = false;
       	if (g[target(*out_i, g)].index == g[tip1_v].index || g[target(*out_i, g)].index == g[tip2_v].index) found = true;
       	++out_i;
       	if (g[target(*out_i, g)].index == g[tip1_v].index || g[target(*out_i, g)].index == g[tip2_v].index) found = true;
       	if (found)  return root_v;
       	else{
       		Graph::vertex_descriptor firstit = get_LCA(target(*out_i, g), tip1_v, tip2_v);
       		--out_i;
       		Graph::vertex_descriptor lastit = get_LCA(target(*out_i, g), tip1_v, tip2_v);
       		if (firstit && lastit) return root_v;
       		else if (firstit) return firstit;
       		else return lastit;
       	}
}

void PopGraph::set_node_heights(vector<Graph::vertex_descriptor> trav){
        	for (int i = 0; i < trav.size(); i++){
        		g[trav[i]].height = get_dist_to_root(trav[i]);
        	}
}

map<string, Graph::vertex_descriptor> PopGraph::get_tips(Graph::vertex_descriptor p_rootIterator){
	map<string, Graph::vertex_descriptor> toreturn;
  	if (out_degree(p_rootIterator, g) ==0){
 		toreturn.insert(make_pair(g[p_rootIterator].name, p_rootIterator));
 	}
 	else{
 		graph_traits<Graph>::out_edge_iterator out_i = out_edges(p_rootIterator, g).first;
 		map<string, Graph::vertex_descriptor> t1 = get_tips(target(*out_i, g));


 		++out_i;
 		map<string, Graph::vertex_descriptor> t2 = get_tips(target(*out_i, g));
 		for (map<string, Graph::vertex_descriptor >::iterator it1 = t1.begin(); it1 != t1.end(); it1++){
 			//pair<int, iterator<NODEDATA> > = std::make_pair(it1->first, it1->second);
 			toreturn.insert(make_pair(it1->first, it1->second));
 		}
 		for (map<string, Graph::vertex_descriptor>::iterator it2 = t2.begin(); it2 != t2.end(); it2++){
 			toreturn.insert(make_pair(it2->first, it2->second));
 		}
 	}
 	return toreturn;
}


void PopGraph::perturb_node_heights(vector< Graph::vertex_descriptor > trav, double epsilon, gsl_rng *r){

  	// perturb the tip heights, they're in even positions
  	for(int i = 0 ; i < trav.size(); i+=2){
  		double d1 = 0;
  		double d2 = 0;
  		double max = 0;
  		if (i-1 >=0)	d1 = g[trav[i-1]].height;
  		if (i+1 < trav.size())	d2 = g[trav[i+1]].height;

  		max = d1;
  		if (d2 > max ) max = d2;
  		double toadd = (2* gsl_rng_uniform(r) - 1)*epsilon;
  		double newheight = g[trav[i]].height + toadd;
  		if (newheight < max) newheight = max+ (max-newheight);
  		g[trav[i]].height = newheight;
  	}

  	//perturb the heights of the interior nodes
  	for (int i = 1 ; i < trav.size(); i+=2){
     		double d1 = 10000;
     		double d2 = 10000;
     		double min = 0;
     		if (i-1 >=0)	d1 = g[trav[i-1]].height;
      		if (i+1 < trav.size())	d2 = g[trav[i+1]].height;

    		min = d1;
    		if (d2 < min ) min = d2;
    		double toadd = (2* gsl_rng_uniform(r) - 1)*epsilon;
    		double newheight = fabs(g[trav[i]].height + toadd);
    		if (newheight > min) newheight = min - (newheight-min);
    		if(newheight > 0) g[trav[i]].height = newheight;
  	}
  }



void PopGraph::build_tree(vector< Graph::vertex_descriptor > trav){
 	//
 	//first find the new root (the minimum height)
 	//

 	Graph::vertex_descriptor newroot = trav[0];
 	double minheight = g[trav[0]].height;
 	int minpos = 0;
 	for(int i = 0 ; i < trav.size(); i ++){
 		if (g[trav[i]].height < minheight) {
 			newroot = trav[i];
 			minheight = g[trav[i]].height;
 			minpos =i;
 		}
 	}
 	set_root(newroot);

 	build_tree_helper(&trav, minpos);


 }

void PopGraph::set_root(Graph::vertex_descriptor v){
	g[root].is_root = false;
	g[v].is_root = true;
	root = v;
}
void PopGraph::build_tree_helper(vector< Graph::vertex_descriptor >* trav, int index){
 	//
 	// look left
 	//
 	bool leftborder = false;
 	bool rightborder = false;
 	double leftmin = 10000;
 	double rightmin = 10000;
 	int minindex = 0;
 	int i = index-1;
 	bool foundleft = false;
 	while (i >= 0 && leftborder ==false){
 		if (g[trav->at(i)].height < g[trav->at(index)].height){
 			leftborder = true;
 			continue;
 		}
 		if (g[trav->at(i)].height < leftmin){
 			minindex = i;
 			leftmin = g[trav->at(i)].height;
 			foundleft = true;
 		}
 		i--;
 	}
  	//
 	// fiddle with the adjacency list if necessary
 	//
 	if (foundleft){
 		clear_out_edges(trav->at(index), g);
 		add_edge(trav->at(index), trav->at(minindex), g);
 		build_tree_helper(trav, minindex);
 	}

 	//
 	// look right
 	//
 	minindex = 0;
 	i = index+1;
 	bool foundright = 0;
   	while (i < trav->size() && rightborder ==false){
     		if (g[trav->at(i)].height < g[trav->at(index)].height){
     			rightborder = true;
     			continue;
     		}
     		if (g[trav->at(i)].height < rightmin){
     			minindex = i;
     			rightmin = g[trav->at(i)].height;
     			foundright = true;
     		}
     		i++;
   	}
   	//
   	// fiddle with the adjacency list if necessary
   	//
   	if (foundright){
		add_edge(trav->at(index), trav->at(minindex), g);
 		build_tree_helper(trav, minindex);
   	}
 }


void PopGraph::update_branch_lengths(Graph::vertex_descriptor root_v){
		graph_traits<Graph>::out_edge_iterator out_i = out_edges(root_v, g).first;
		graph_traits<Graph>::out_edge_iterator out_end = out_edges(root_v, g).second;
		graph_traits<Graph>::in_edge_iterator in_i = in_edges(root_v, g).first;
		if (out_i != out_end) update_branch_lengths(target(*out_i, g));
		if (g[root_v].is_root == false) {
			g[*in_i].len = g[root_v].height - g[source(*in_i, g)].height;
		}
 		if (out_i != out_end) {
 			++out_i;
 			update_branch_lengths(target(*out_i, g));
 		}
  }


void PopGraph::randomize_tree(gsl_rng* r){
	map<string, Graph::vertex_descriptor> popname2tip = get_tips(root);
  	vector< Graph::vertex_descriptor > trav = get_inorder_traversal( popname2tip.size());
 	double tipheight = gsl_rng_uniform(r)*0.5+0.5;

 	// set the tip heights
 	for(int i = 0 ; i < trav.size(); i+=2){
 		g[trav[i]].height = tipheight;
 	}
 	for(int i = 1; i < trav.size(); i+=2){
 		g[trav[i]].height = gsl_rng_uniform(r)*tipheight;
 	}

 	build_tree(trav);
 	update_branch_lengths(root);
 	trav = get_inorder_traversal(popname2tip.size());
 	set_node_heights(trav);
 }

string PopGraph::get_newick_format(){
 	string toreturn = "";
 	newick_helper(root, &toreturn);
 	return toreturn;
 }



void PopGraph::print(){
	IndexMap index = get(&Node::index, g);

	cout << "vertices(g) = ";
	pair<vertex_iter, vertex_iter> vp;
	for (vp = vertices(g); vp.first != vp.second; ++vp.first)
		std::cout << index[*vp.first] <<  ":"<< g[*vp.first].name << ":"<< get_dist_to_root(*vp.first)<< " ";
	cout << "\n";

    std::cout << "edges(g) = ";
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei)
        std::cout << "(" << index[source(*ei, g)]
                  << "," << index[target(*ei, g)] << ","<< g[*ei].len <<  ") ";
    std::cout << std::endl;

}


void PopGraph::newick_helper(Graph::vertex_descriptor node, string* s){
	if (node != NULL){
		graph_traits<Graph>::out_edge_iterator out_i = out_edges(node, g).first;
		graph_traits<Graph>::out_edge_iterator out_end = out_edges(node, g).second;
		if (out_i != out_end){
			s->append("(");
			newick_helper(target(*out_i, g), s);
			s->append(",");
			++out_i;
			newick_helper(target(*out_i, g), s);
			s->append(")");
			if (g[node].is_root == false){
				graph_traits<Graph>::in_edge_iterator in_i = in_edges(node, g).first;
				s->append(":");
				stringstream ss;
				ss<< g[*in_i].len;
				s->append(ss.str());
			}
			else s->append(";");
		}
		else{
			graph_traits<Graph>::in_edge_iterator in_i = in_edges(node, g).first;
			stringstream ss;
			ss << g[node].name;
			ss << ":";
			ss<< g[*in_i].len;
			s->append(ss.str());
		}
	}
}

double PopGraph::get_height(Graph::vertex_descriptor v){
	return g[v].height;
}
