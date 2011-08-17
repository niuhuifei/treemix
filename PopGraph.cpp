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


PopGraph::PopGraph(vector<string> first3pops){

	g = Graph(1);
	istree = true;
	Graph::vertex_descriptor v = *vertices(g).first;
	Graph::vertex_descriptor v2, v3, v4, v5, v6;
	Graph::edge_descriptor e;
	g[v].index = 0;
	g[v].name = "NA";
	g[v].height = 0;
	g[v].is_tip = false;
	g[v].is_root = true;
	g[v].rev = false;
	g[v].is_mig = false;
	g[v].mig_frac = 0;
	root = v;

	// add the first three populations to the tree
	// tree will look like this:
	//
	//           root
	//            /\
	//           /  \
	//          /   /\
	//         1    2 3
	if (first3pops.size() != 3) {
		std::cout << "Error initializing PopGraph: need 3 pops, not "<< first3pops.size() << "\n";
		exit(1);
	}
	v3 = add_vertex(g);
	g[v3].index = 1;
	g[v3].name = first3pops[0];
	g[v3].height = 0;
	g[v3].is_tip = true;
	g[v3].is_root = false;
	g[v3].rev = false;
	g[v3].is_mig = false;
	g[v].mig_frac = 0;
	e = add_edge( v, v3, g).first;
	g[e].weight = 1;
	g[e].len = 1;
	g[e].is_mig = false;


	v4 = add_vertex(g);
	g[v4].index = 2;
	g[v4].name = "NA";
	g[v4].height = 0;
	g[v4].is_tip = false;
	g[v4].is_root = false;
	g[v4].rev = false;
	g[v4].is_mig = false;
	g[v4].mig_frac = 0;
	e = add_edge( v, v4, g).first;
	g[e].weight = 1;
	g[e].len = 1;
	g[e].is_mig = false;


	v5 = add_vertex(g);
	g[v5].index = 3;
	g[v5].name = first3pops[1];
	g[v5].height = 0;
	g[v5].is_tip = true;
	g[v5].is_root = false;
	g[v5].rev = false;
	g[v5].is_mig = false;
	g[v5].mig_frac = 0;
	e = add_edge( v4, v5, g).first;
	g[e].weight = 1;
	g[e].len = 1;
	g[e].is_mig = false;

	v6 = add_vertex(g);
	g[v6].index = 4;
	g[v6].name = first3pops[2];
	g[v6].height = 0;
	g[v6].is_tip = true;
	g[v6].is_root = false;
	g[v6].rev = false;
	g[v6].is_mig = false;
	g[v6].mig_frac = 0;
	e = add_edge( v4, v6, g).first;
	g[e].weight = 1;
	g[e].len = 1;
	g[e].is_mig = false;
	indexcounter = 5;
}



Graph::vertex_descriptor PopGraph::add_tip(Graph::vertex_descriptor v, string name){
	// add a tip attached above v, halfway to the parent node
	Graph::vertex_descriptor vtipnew, vintnew, vparent;
	double len, newlen;
	Graph::edge_descriptor e;

	if (g[v].is_root == false){
		// remove edge from parent to target
		graph_traits<Graph>::in_edge_iterator init = in_edges(v, g).first;
		len = g[*init].len;
		newlen = len/2;
		vparent = source(*init, g);
		remove_edge(*init, g);

		//add new vertices
		vtipnew = add_vertex(g);
		vintnew = add_vertex(g);

		g[vtipnew].index = indexcounter;
		g[vtipnew].name = name;
		g[vtipnew].height = 0;
		g[vtipnew].is_tip = true;
		g[vtipnew].is_root = false;
		g[vtipnew].rev = false;
		g[vtipnew].is_mig = false;
		g[vtipnew].mig_frac = 0;
		indexcounter++;

		g[vintnew].index = indexcounter;
		g[vintnew].name = "NA";
		g[vintnew].height = 0;
		g[vintnew].is_tip = false;
		g[vintnew].is_root = false;
		g[vintnew].rev = false;
		g[vintnew].is_mig = false;
		g[vintnew].mig_frac = 0;
		indexcounter++;

		// add new edges
		e = add_edge( vparent, vintnew, g).first;
		g[e].weight = 1;
		g[e].len = newlen;
		g[e].is_mig = false;
		e = add_edge( vintnew, v, g).first;
		g[e].weight = 1;
		g[e].len = newlen;
		g[e].is_mig = false;
		e = add_edge( vintnew, vtipnew, g).first;
		g[e].weight = 1;
		g[e].len = 1;
		g[e].is_mig = false;
		return vtipnew;
	}
	else{

		//add a new root
		g[v].is_root = false;
		vintnew = add_vertex(g);
		vtipnew = add_vertex(g);
		root = vintnew;
		g[vtipnew].index = indexcounter;
		g[vtipnew].name = name;
		g[vtipnew].height = 0;
		g[vtipnew].is_tip = true;
		g[vtipnew].is_root = false;
		g[vtipnew].rev = false;
		g[vtipnew].is_mig = false;
		g[vtipnew].mig_frac = 0;
		indexcounter++;

		g[vintnew].index = indexcounter;
		g[vintnew].name = "NA";
		g[vintnew].height = 0;
		g[vintnew].is_tip = false;
		g[vintnew].is_root = true;
		g[vintnew].rev = false;
		g[vintnew].is_mig = false;
		g[vintnew].mig_frac = 0;
		indexcounter++;

		// add new edges
		e = add_edge( vintnew, v, g).first;
		g[e].weight = 1;
		g[e].len = 1;
		g[e].is_mig = false;
		e = add_edge( vintnew, vtipnew, g).first;
		g[e].weight = 1;
		g[e].len = 1;
		g[e].is_mig = false;
		return vtipnew;
	}
}

Graph::edge_descriptor PopGraph::add_mig_edge(Graph::vertex_descriptor st, Graph::vertex_descriptor sp){
	//
	//
	//
	//
	//
	// from      p1
	//          /         \
	//         /           \
	//        st            sp
	// to
	//
	//          p1
    //         /
	//        /
	//       p2      /
	//      /  \    /
	//     /    \  /
	//    st     sp
 	//
	//
/*
	Graph::edge_descriptor e;
	istree = false;
	e = add_edge(st, sp, g).first;
	g[e].weight = 0.25;
	g[e].len = 1;
	g[e].is_mig = true;
	return e;
*/

	Graph::vertex_descriptor p1, p2;
	Graph::edge_descriptor e;
	Graph::in_edge_iterator in_it = in_edges(st, g).first;
	p1 = source(*in_it, g);
	e = *in_it;
	double oldlen = g[e].len;
	//remove edge, insert node
	remove_edge(e, g);
	p2 = add_vertex(g);
	g[p2].index = indexcounter;
	g[p2].name = "NA";
	g[p2].height = 0;
	g[p2].is_tip = false;
	g[p2].is_root = false;
	g[p2].rev = false;
	g[p2].is_mig = true;
	g[p2].mig_frac = 0.5; //default to migration halfway between nodes
	indexcounter++;

	//add new edges
	e = add_edge(p1, p2, g).first;
	g[e].weight = 1;
	g[e].is_mig = false;
	if (g[p1].is_mig) g[e].len = 0;
	else g[e].len = oldlen/2;

	e = add_edge(p2, st, g).first;
	g[e].weight = 1;
	g[e].is_mig = false;
	if (g[p1].is_mig) g[e].len = oldlen;
	else g[e].len = oldlen/2;

	e = add_edge(p2, sp, g).first;
	g[e].weight = 0.25;
	g[e].len = 0;
	g[e].is_mig = true;

	return e;

}

set<pair<double, set<Graph::vertex_descriptor> > > PopGraph::get_paths_to_root(Graph::vertex_descriptor v){

	set<pair<double, set<Graph::vertex_descriptor> > > toreturn;


	double wsum = 0;

	Graph::edge_descriptor nonmig;
	for (Graph::in_edge_iterator it = in_edges(v, g).first; it != in_edges(v, g).second; it++){
		if (g[*it].is_mig == false) nonmig = *it;
		else wsum+= g[*it].weight;
	}
	g[nonmig].weight = 1-wsum;
	Graph::in_edge_iterator init = in_edges(v, g).first;
	while (init != in_edges(v, g).second){

		double w = g[*init].weight;
		//cout << "here2 " << w << " "<< g[v].index << "\n"; cout.flush();
		cout.flush();
		set<Graph::vertex_descriptor> tmpset;
		tmpset.insert(v);
		Graph::vertex_descriptor parent = source(*init, g);
		tmpset.insert(parent);
		if (g[parent].is_root){
			toreturn.insert(make_pair(w, tmpset));
		}
		else{
			set<pair<double, set<Graph::vertex_descriptor> > > p_path = get_paths_to_root(parent);
			for (set<pair<double, set<Graph::vertex_descriptor> > >::iterator it = p_path.begin(); it != p_path.end(); it++){
				double w2 = w*it->first;
				set<Graph::vertex_descriptor> tmpset2 = tmpset;
				for (set<Graph::vertex_descriptor>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) tmpset2.insert(*it2);
				toreturn.insert(make_pair(w2, tmpset2));
			}
		}
		init++;
	}
	return toreturn;
}


set<pair<double, set<Graph::edge_descriptor> > > PopGraph::get_paths_to_root_edge(Graph::vertex_descriptor v){

	set<pair<double, set<Graph::edge_descriptor> > > toreturn;


	double wsum = 0;
	Graph::edge_descriptor nonmig;
	for (Graph::in_edge_iterator it = in_edges(v, g).first; it != in_edges(v, g).second; it++){
		if (g[*it].is_mig == false) nonmig = *it;
		else wsum+= g[*it].weight;
	}
	g[nonmig].weight = 1-wsum;

	Graph::in_edge_iterator init = in_edges(v, g).first;
	while (init != in_edges(v, g).second){

		double w = g[*init].weight;
		//scout << "here2 " << w << " "<< g[v].index << "\n"; cout.flush();
		cout.flush();
		set<Graph::edge_descriptor> tmpset;
		tmpset.insert(*init);
		Graph::vertex_descriptor parent = source(*init, g);
		if (g[parent].is_root){
			toreturn.insert(make_pair(w, tmpset));
		}
		else{
			set<pair<double, set<Graph::edge_descriptor> > > p_path = get_paths_to_root_edge(parent);
			for (set<pair<double, set<Graph::edge_descriptor> > >::iterator it = p_path.begin(); it != p_path.end(); it++){
				double w2 = w*it->first;
				set<Graph::edge_descriptor> tmpset2 = tmpset;
				for (set<Graph::edge_descriptor>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++) tmpset2.insert(*it2);
				toreturn.insert(make_pair(w2, tmpset2));
			}
		}
		init++;
	}
	return toreturn;
}

//pair<double, set<Graph::vertex_descriptor> > PopGraph::get_paths_helper(Graph::vertex_descriptor v){
//	pair<double, set<Graph::vertex_descriptor> > toreturn;
//	return toreturn;
//}

void PopGraph::remove_tip(Graph::vertex_descriptor v){
	Graph::vertex_descriptor vparent, vparent2, vnewdest;
	double newlen;
	Graph::edge_descriptor e;


	graph_traits<Graph>::in_edge_iterator init = in_edges(v, g).first;
	vparent = source(*init, g);

	//if the node isn't coming from the root
	if (g[vparent].is_root == false){

		/*
		 *     vp2
		 *     / \
		 *    /   \
		 *         vp
		 *         / \
		 *        /   \
		 *       v     vnewdest
		 *
		 *
		 *  delete v, vp, draw edge from vp2 to vnewdest
		 */
		graph_traits<Graph>::in_edge_iterator init2 = in_edges(vparent, g).first;
		vparent2 = source(*init2, g);

		remove_edge(*init, g);

		graph_traits<Graph>::out_edge_iterator outit = out_edges(vparent, g).first;
		vnewdest = target(*outit, g);
		newlen = g[*outit].len+ g[*init2].len;

		e = add_edge(vparent2, vnewdest, g).first;
		g[e].len = newlen;
		g[e].weight = 1;

		clear_vertex(vparent, g);
		remove_vertex(vparent, g);
		remove_vertex(v, g);
	}

	//otw

	else{

		/*
		 *
		 *         vp
		 *        /  \
		 *       /    \
		 *      v      vnewdest
		 *
		 * delete v, vp, make vnewdest the new root
		 *
		 */
		remove_edge(*init, g);
		graph_traits<Graph>::out_edge_iterator outit = out_edges(vparent, g).first;
		vnewdest = target(*outit, g);
		g[vnewdest].is_root = true; //this is the new root
		root = vnewdest;
		clear_vertex(vparent, g);
		remove_vertex(vparent, g);
		remove_vertex(v, g);
 	}
}


void PopGraph::copy(PopGraph * s){
	istree = s->istree;
	popnames = s->popnames;


	//cout << "copying graph\n"; cout.flush();
	g.clear();
	map<int, Graph::vertex_descriptor> tmpmap;
	Graph::vertex_descriptor vd;
	pair<Graph::vertex_iterator, Graph::vertex_iterator> v = vertices(s->g);
	pair<Graph::edge_iterator, Graph::edge_iterator> e = edges(s->g);
	for (Graph::vertex_iterator it = v.first; it != v.second; ++it){
		vd = add_vertex(g);
		g[vd].index = s->g[*it].index;
		g[vd].name = s->g[*it].name;
		g[vd].height = s->g[*it].height;
		g[vd].rev = s->g[*it].rev;
		g[vd].is_root = s->g[*it].is_root;
		g[vd].is_tip = s->g[*it].is_tip;
		g[vd].is_mig = s->g[*it].is_mig;
		tmpmap.insert(make_pair( g[vd].index, vd));
		if (s->g[*it].is_root == true) root = vd;
	}

	for (Graph::edge_iterator it = e.first; it != e.second; ++it){
		Graph::vertex_descriptor tmpsource = source(*it, s->g);
		Graph::vertex_descriptor tmptarget = target(*it, s->g);
		Graph::edge_descriptor ed = add_edge( tmpmap[ s->g[tmpsource].index ], tmpmap[ s->g[tmptarget].index ], g ).first;
		g[ed].weight = s->g[*it].weight;
		g[ed].len = s->g[*it].len;
		g[ed].is_mig = s->g[*it].is_mig;
	}
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
	g[v].is_mig = false;
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
			g[v2].is_mig = false;
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
			g[v2].is_mig = false;
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
            	  g[*in_i].is_mig = false;
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



void PopGraph::set_graph(string p_newickString){
	g.clear();
	add_vertex(g);
	index2father.clear();
	popnames.clear();
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
	g[v].is_mig = false;
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
			g[v2].is_mig = false;
			e = add_edge( v, v2, g).first;
			v = v2;
		}
		else if( *I == ')' )// Subtree finished, get back
		{
			v = index2father[g[v].index];
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
			g[v2].is_mig = false;
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
            	  g[*in_i].is_mig = false;
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
	indexcounter = i+1;
}


set<Graph::vertex_descriptor> PopGraph::get_root_adj(){
	set<Graph::vertex_descriptor> toreturn;
	pair<Graph::vertex_descriptor, Graph::vertex_descriptor> adj = get_child_nodes(root);
	toreturn.insert( adj.first);
	toreturn.insert( adj.second);
	return toreturn;
}

double PopGraph::get_dist_to_root(Graph::vertex_descriptor v){
	double toreturn = 0;
	if (g[v].is_root) return 0;
	set<pair<double, set<Graph::edge_descriptor> > > paths = get_paths_to_root_edge(v);
	for (set<pair<double, set<Graph::edge_descriptor> > >::iterator it = paths.begin(); it!= paths.end(); it++){
		double w = it->first;
		double l = 0;
		for (set<Graph::edge_descriptor>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++){
			l += g[*it2].len;
		}
		toreturn += w*l;
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

vector<Graph::vertex_descriptor > PopGraph::get_inorder_traversal_noroot(int nodes){
	vector<Graph::vertex_descriptor> toreturn;
	vector<Graph::vertex_descriptor> tmp = get_inorder_traversal(nodes);
	for (int i = 0; i < tmp.size(); i++){
		if (g[tmp[i]].is_root == false) toreturn.push_back(tmp[i]);
	}
	return toreturn;
}
void PopGraph::inorder_traverse(Graph::vertex_descriptor rootIterator, int* i, vector<Graph::vertex_descriptor >* v){
	pair<Graph::vertex_descriptor, Graph::vertex_descriptor> child_nodes = get_child_nodes(rootIterator);
	//graph_traits<Graph>::out_edge_iterator out_i = out_edges(rootIterator, g).first;
	//cout << g[rootIterator].index << " "<< g[rootIterator].name << " "<< g[rootIterator].is_tip << " "<< g[rootIterator].rev << "\n";
	if (g[rootIterator].is_tip == false){
		//inorder_traverse(target(*out_i, g), i , v);
		inorder_traverse(child_nodes.first, i , v);
 	}

 	v->at(*i) = rootIterator;
 	*i = *i+1;

 	if (g[rootIterator].is_tip == false){
 		inorder_traverse(child_nodes.second, i, v);
 	}
}

vector<Graph::edge_descriptor> PopGraph::get_mig_edges(){
	vector<Graph::edge_descriptor> toreturn;
	for (Graph::edge_iterator it = edges(g).first; it != edges(g).second; it++){
		if (g[*it].is_mig == true) toreturn.push_back(*it);
	}
	return toreturn;

}
Graph::vertex_descriptor PopGraph::get_LCA(Graph::vertex_descriptor root_v,
		Graph::vertex_descriptor tip1_v, Graph::vertex_descriptor tip2_v){
		//cout << "in get_LCA "<< g[root_v].index << " "<< g[tip1_v].index << " "<< g[tip2_v].index << "\n"; cout.flush();
		if (g[root_v].is_tip == true) return NULL;
		pair<Graph::vertex_descriptor, Graph::vertex_descriptor> children = get_child_nodes(root_v);
		Graph::vertex_descriptor c1 = children.first;
		Graph::vertex_descriptor c2 = children.second;
		graph_traits<Graph>::out_edge_iterator out_i = out_edges(root_v, g).first;
		bool found = false;
       	if (g[c1].index == g[tip1_v].index || g[c1].index == g[tip2_v].index) found = true;
       	if (g[c2].index == g[tip1_v].index || g[c2].index == g[tip2_v].index) found = true;
       	if (found)  return root_v;
       	else{
       		Graph::vertex_descriptor firstit = get_LCA(c2, tip1_v, tip2_v);
       		Graph::vertex_descriptor lastit = get_LCA(c1, tip1_v, tip2_v);
       		if (firstit && lastit) return root_v;
       		else if (firstit) return firstit;
       		else return lastit;
       	}
}


map<string, Graph::vertex_descriptor> PopGraph::get_tips(Graph::vertex_descriptor p_rootIterator){
	map<string, Graph::vertex_descriptor> toreturn;
  	if (out_degree(p_rootIterator, g) ==0){
 		toreturn.insert(make_pair(g[p_rootIterator].name, p_rootIterator));
 	}
 	else{
 		pair<Graph::vertex_descriptor, Graph::vertex_descriptor> children =  get_child_nodes(p_rootIterator);
 		map<string, Graph::vertex_descriptor> t1 = get_tips(children.first);


 		map<string, Graph::vertex_descriptor> t2 = get_tips(children.second);
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

void PopGraph::move_root(gsl_rng* r){
	double ran = gsl_rng_uniform(r);
	Graph::out_edge_iterator it = out_edges(root, g).first;
	Graph::vertex_descriptor target;
	if (ran < 0.5) it++;

	//target = target(*it, g);


}


void PopGraph::local_rearrange(Graph::vertex_descriptor v3, int i){
	/*
	 *
	 * go from    1
	 *          /  \*
	 *         /    3
	 *        /    / \
	 *       2    4   5
	 *
	 *  if i == 1, go to
	 *
	 *           1
	 *          / \
	 *         3   \
	 *        / \   \
	 *       2   4   5
	 *
	 *  if i == 2, go to
	 *
	 *           1
	 *          / \
	 *         3   \
	 *        / \   \
	 *       2   5   4
	 *
	 *
	 */
	if (!g[v3].is_root){

		// set up descriptors and edge lengths
		Graph::in_edge_iterator tmp = in_edges(v3, g).first;
		Graph::edge_descriptor ed = *tmp;
		Graph::vertex_descriptor v1, v2, v4, v5;
		double d12, d13, d34, d35;
		double l2, l4, l5;
		d13 = g[ed].len;
		v1 = source(ed, g);
		Graph::out_edge_iterator out3 = out_edges(v3, g).first;
		v4 = target(*out3, g);
		d34 = g[*out3].len;

		++out3;
		v5 = target(*out3, g);
		d35 = g[*out3].len;

		Graph::out_edge_iterator out1 = out_edges(v1, g).first;
		if (target(*out1, g) == v3) {
			out1++;
			v2 = target(*out1, g);
			d12 = g[*out1].len;
		}
		else {
			v2 = target(*out1, g);
			d12 = g[*out1].len;
		}
		l2 = d12;
		l4 = d13+d34;
		l5 = d13+d35;
		double min = l2;
		if (l4 < min) min = l4;
		if (l5 < min) min = l5;

		if (i == 1){
		// go to ((2,4),5)

			remove_edge(v1, v2, g);
			remove_edge(v3, v5, g);
			Graph::edge_descriptor e = add_edge(v3, v2, g).first;
			g[e].weight = 1;
			g[e].len = d12;
			g[e].is_mig = false;
			e = add_edge(v1, v5, g).first;
			g[e].weight = 1;
			g[e].len = d35;
			g[e].is_mig = false;

		}
		else{
			//go to (4,(2,5))

			remove_edge(v1, v2, g);
			remove_edge(v3, v4, g);
			Graph::edge_descriptor e = add_edge(v3, v2, g).first;
			g[e].weight = 1;
			g[e].len = d12;
			g[e].is_mig = false;

			e = add_edge(v1, v4, g).first;
			g[e].weight = 1;
			g[e].len = d34;
			g[e].is_mig = false;
		}
	}

}

void PopGraph::move_root(int i){
	/*
	 *  go from
	 *        root
	 *         / \
	 *        /   \
	 *       v1    v2
	 *      / \    / \
	 *     /   \  v5  v6
	 *    v3   v4
	 *

	 *
	 */

	Graph::vertex_descriptor v1, v2, v3, v4, v5, v6;
	pair<Graph::vertex_descriptor, Graph::vertex_descriptor> children = get_child_nodes(root);
	Graph::edge_descriptor er1, er2;
	double lr1, lr2;
	v1 = children.first;
	v2 = children.second;
	if ( !g[v1].is_tip ){
		children = get_child_nodes(v1);
		v3 = children.first;
		v4 = children.second;
	}
	if ( !g[v2].is_tip ){
		children = get_child_nodes(v2);
		v5 = children.first;
		v6 = children.second;
	}
	er1 = edge(root, v1, g).first;
	er2 = edge(root, v2, g).first;
	lr1 = g[er1].len;
	lr2 = g[er2].len;

	 /*      to
		 *
		 *   i == 1
		 *
		 *    root
		 *     / \
		 *    /   \
		 *   v3    v1
		 *        / \
		 *       /   \
		 *      v4    v2
*/

	if ( i == 1 && !g[v1].is_tip){
		double l13;
		Graph::edge_descriptor e13 = edge(v1, v3, g).first;
		l13 = g[e13].len;
		remove_edge( er1, g);
		remove_edge( er2, g);
		remove_edge( e13, g);
		Graph::edge_descriptor e = add_edge(v1, v2, g).first;
		g[e].len = lr1+lr2;
		g[e].is_mig = false;
		g[e].weight = 1;

		e = add_edge(root, v3, g).first;
		g[e].len = l13/2;
		g[e].is_mig = false;
		g[e].weight = 1;

		e = add_edge(root, v1, g).first;
		g[e].len = l13/2;
		g[e].is_mig = false;
		g[e].weight = 1;
	}

	 /*
	 *   i == 2
	 *
	 *    root
	 *     / \
	 *    /   \
	 *   v4    v1
	 *        / \
	 *       /   \
	 *      v3    v2
	 */
	else if ( i == 2 && !g[v1].is_tip){
		double l14;
		Graph::edge_descriptor e14 = edge(v1, v4, g).first;
		l14 = g[e14].len;
		remove_edge( er1, g);
		remove_edge( er2, g);
		remove_edge( e14, g);
		Graph::edge_descriptor e = add_edge(v1, v2, g).first;
		g[e].len = lr1+lr2;
		g[e].is_mig = false;
		g[e].weight = 1;

		e = add_edge(root, v4, g).first;
		g[e].len = l14/2;
		g[e].is_mig = false;
		g[e].weight = 1;

		e = add_edge(root, v1, g).first;
		g[e].len = l14/2;
		g[e].is_mig = false;
		g[e].weight = 1;

	}

	/*   i == 3
	 *
	 *    root
	 *     / \
	 *    /   \
	 *   v5    v2
	 *        / \
	 *       /   \
	 *      v1    v6
	 */
	else if ( i == 3 && !g[v2].is_tip){
		double l25;
		Graph::edge_descriptor e25 = edge(v2, v5, g).first;
		l25 = g[e25].len;
		remove_edge( er1, g);
		remove_edge( er2, g);
		remove_edge( e25, g);
		Graph::edge_descriptor e = add_edge(v2, v1, g).first;
		g[e].len = lr1+lr2;
		g[e].is_mig = false;
		g[e].weight = 1;

		e = add_edge(root, v5, g).first;
		g[e].len = l25/2;
		g[e].is_mig = false;
		g[e].weight = 1;

		e = add_edge(root, v2, g).first;
		g[e].len = l25/2;
		g[e].is_mig = false;
		g[e].weight = 1;

	}


	 /*   i == 4
	 *
	 *    root
	 *     / \
	 *    /   \
	 *   v6    v2
	 *        / \
	 *       /   \
	 *      v1    v5
	 */

	else if ( i == 4 && !g[v2].is_tip){
		double l26;
		Graph::edge_descriptor e26 = edge(v2, v6, g).first;
		l26 = g[e26].len;
		remove_edge( er1, g);
		remove_edge( er2, g);
		remove_edge( e26, g);
		Graph::edge_descriptor e = add_edge(v2, v1, g).first;
		g[e].len = lr1+lr2;
		g[e].is_mig = false;
		g[e].weight = 1;

		e = add_edge(root, v6, g).first;
		g[e].len = l26/2;
		g[e].is_mig = false;
		g[e].weight = 1;

		e = add_edge(root, v2, g).first;
		g[e].len = l26/2;
		g[e].is_mig = false;
		g[e].weight = 1;

	}
}

void PopGraph::local_rearrange_wmig(Graph::vertex_descriptor v){
	/*
	 *             vppp
	 *            /  \
	 *           /
	 *          vpp
	 *    \    /  \
	 *     \  /    \
	 *       vp    vu
	 *  \   / \
	 *   \ /   \
	 *    v     vs
	 *
	 *
	 *    to
	 *
	 *           vppp
	 *            /
	 *           /
	 *         vnew
	 *    \    / \
	 *     \  /   \
	 *      v     vpp
	 *         \   /
	 *          \ /
	 *           vs
	 */

	Graph::vertex_descriptor vp, vpp, vppp, vs;
}


void PopGraph::global_rearrange(Graph::vertex_descriptor v1, Graph::vertex_descriptor v2){
	/*
	 *   take the tree below v1p, attach it above v2
	 *   if v2 is the root, there's a special case
	 *   do not do anything is v1 or v1p is the root
	 *
	 *   otherwise:
	 *
	 *             v1pp
	 *             /  \
	 *            /    \
	 *          v1p                  v2p
	 *          / \                  / \
	 *         /   \                /   \
	 *        v1   v1s            v2
	 *       / \
	 *      /   \
	 *
	 *      goes to
	 *
	 *
	 *      v1pp                     v2p
	 *      / \                      /  \
	 *     /   \                    /    \
	 *    v1s                      v1p
	 *                            / \
	 *                           /   \
	 *                          v2   v1
	 */
	 Graph::vertex_descriptor v1p, v1pp, v1s, v2p;
	 Graph::edge_descriptor e;
	 if (!g[v1].is_root && !g[v2].is_root && !g[get_parent_node(v1).first].is_root){
		 v1p = get_parent_node(v1).first;
		 v1pp = get_parent_node(v1p).first;
		 v2p = get_parent_node(v2).first;
		 pair<Graph::vertex_descriptor, Graph::vertex_descriptor> d_v1p = get_child_nodes(v1p);
		 if (d_v1p.first == v1) v1s = d_v1p.second;
		 else v1s = d_v1p.first;

		 // take care of the tree with v1
		 remove_edge(v1pp, v1p, g);
		 remove_edge(v1p, v1s, g);
		 e = add_edge(v1pp, v1s, g).first;
		 g[e].weight = 1;
		 g[e].len = 1;
		 g[e].is_mig = false;

		 //and the tree with v2
		 remove_edge(v2p, v2, g);
		 e = add_edge(v2p, v1p, g).first;
		 g[e].weight = 1;
		 g[e].len = 1;
		 g[e].is_mig = false;

		 e = add_edge(v1p, v2, g).first;
		 g[e].weight = 1;
		 g[e].len = 1;
		 g[e].is_mig = false;
	 }
	 else if (g[v2].is_root){
		 /*
		  *      go to
		  *
		  *
		  *      v1pp
		  *      / \
		  *     /   \
		  *    v1s                      v1p
		  *                            / \
		  *                           /   \
		  *                          v2   v1
		  *
		  *
		  *
		  */
		 v1p = get_parent_node(v1).first;
		 v1pp = get_parent_node(v1p).first;
		 pair<Graph::vertex_descriptor, Graph::vertex_descriptor> d_v1p = get_child_nodes(v1p);
		 if (d_v1p.first == v1) v1s = d_v1p.second;
		 else v1s = d_v1p.first;

		 // take care of the tree with v1
		 remove_edge(v1pp, v1p, g);
		 remove_edge(v1p, v1s, g);
		 e = add_edge(v1pp, v1s, g).first;
		 g[e].weight = 1;
		 g[e].len = 1;
		 g[e].is_mig = false;

		 //and the tree with v2
		 e = add_edge(v1p, v2, g).first;
		 g[e].weight = 1;
		 g[e].len = 1;
		 g[e].is_mig = false;
		 set_root(v1p);
	 }



}

void PopGraph::set_root(Graph::vertex_descriptor v){
	g[root].is_root = false;
	g[v].is_root = true;
	root = v;
}

string PopGraph::get_newick_format(){
 	string toreturn = "";
 	newick_helper(root, &toreturn);
 	return toreturn;
 }


string PopGraph::get_newick_format(Graph::vertex_descriptor v){
 	string toreturn = "";
 	newick_helper(v, &toreturn);
 	return toreturn;
 }


string PopGraph::get_newick_subtrees(Graph::vertex_descriptor v1, Graph::vertex_descriptor v2){
 	string toreturn1 = "";
 	string toreturn2 = "";
 	stringstream ss;
 	set<Graph::vertex_descriptor> p1 = get_path_to_root(v1);
 	set<Graph::vertex_descriptor> p2 = get_path_to_root(v2);

 	if ( p1.find(v2) != p1.end()){
 	 	newick_helper(v2, &toreturn2, g[v2].index);
 	 	ss<< toreturn2 << ";";
 	 	return ss.str();
 	}
 	else if ( p2.find(v1) != p2.end()){
 	 	newick_helper(v1, &toreturn1, g[v1].index);
 	 	ss<< toreturn1 << ";";
 	 	return ss.str();
 	}
 	else{
 		double d1 = get_dist_to_root(v1);
 		double d2 = get_dist_to_root(v2);
 		double d3 = get_dist_to_root(get_LCA(root, v1, v2));
 		newick_helper(v1, &toreturn1, g[v1].index);
 		newick_helper(v2, &toreturn2, g[v2].index);
 		ss << "(" << toreturn1 << ":"<< d1-d3 << ","<< toreturn2<< ":"<< d2-d3 << ");";
 		//cout << ss.str() <<"\n";
 		return ss.str();
 	}
 }


void PopGraph::print(){
	IndexMap index = get(&Node::index, g);

	cout << "vertices(g) = ";
	pair<vertex_iter, vertex_iter> vp;
	for (vp = vertices(g); vp.first != vp.second; ++vp.first){
		std::cout << index[*vp.first] <<  ":"<< g[*vp.first].name << ":"<< get_dist_to_root(*vp.first);
		if (g[*vp.first].is_root) std::cout << ":ROOT";
		if (g[*vp.first].is_tip) std::cout << ":TIP";
		cout << " ";

	}
	cout << "\n";

    std::cout << "edges(g) = ";
    graph_traits<Graph>::edge_iterator ei, ei_end;
    for (tie(ei, ei_end) = edges(g); ei != ei_end; ++ei){
        std::cout << "(" << index[source(*ei, g)]
                  << "," << index[target(*ei, g)] << ","<< g[*ei].len << ","<< g[*ei].weight;
		if (g[*ei].is_mig) cout << ",MIG";
		cout << ") ";
    }
    std::cout << std::endl;

}


string PopGraph::get_newick_format(map<string, double>* trim){
 	string toreturn = "";
 	newick_helper(root, &toreturn, trim);
 	return toreturn;
 }



void PopGraph::newick_helper(Graph::vertex_descriptor node, string* s){
	if (out_degree(node, g) !=  0){
		pair<Graph::vertex_descriptor, Graph::vertex_descriptor> children = get_child_nodes(node);
		s->append("(");
		newick_helper(children.first, s);
		s->append(",");
		newick_helper(children.second, s);
		s->append(")");
		if (g[node].is_root == false){

			s->append(":");

			stringstream ss;
			ss<< get_parent_node(node).second;
			s->append(ss.str());
		}
		else s->append(";");
	}
	else{
		stringstream ss;
		ss << g[node].name;
		ss << ":";
		ss<< get_parent_node(node).second;
		s->append(ss.str());
	}

}


void PopGraph::newick_helper(Graph::vertex_descriptor node, string* s, int rootindex){
	if (out_degree(node, g) !=  0){
		pair<Graph::vertex_descriptor, Graph::vertex_descriptor> children = get_child_nodes(node);
		s->append("(");
		newick_helper(children.first, s, rootindex);
		s->append(",");
		newick_helper(children.second, s, rootindex);
		s->append(")");
		if (g[node].index != rootindex){

			s->append(":");

			stringstream ss;
			ss<< get_parent_node(node).second;
			s->append(ss.str());
		}
		//else s->append(";");
	}
	else{
		stringstream ss;
		ss << g[node].name;
		if (g[node].index != rootindex) {
			ss << ":";
			ss<< get_parent_node(node).second;
		}
		s->append(ss.str());
	}

}


void PopGraph::newick_helper(Graph::vertex_descriptor node, string* s, map<string, double>* trim){
	if (out_degree(node, g) !=  0){
		pair<Graph::vertex_descriptor, Graph::vertex_descriptor> children = get_child_nodes(node);
		s->append("(");
		newick_helper(children.first, s);
		s->append(",");
		newick_helper(children.second, s);
		s->append(")");
		if (g[node].is_root == false){

			s->append(":");

			stringstream ss;
			ss<< get_parent_node(node).second;
			s->append(ss.str());
		}
		else s->append(";");
	}
	else{
		stringstream ss;
		ss << g[node].name;
		ss << ":";
		double trimmedlen = get_parent_node(node).second - trim->find( g[node].name )->second;
		if (trimmedlen < 0) trimmedlen = 0;
		ss<< trimmedlen;
		s->append(ss.str());
	}
}

double PopGraph::get_height(Graph::vertex_descriptor v){
	return g[v].height;
}


vector<Graph::vertex_descriptor> PopGraph::get_path_to_root_vec(Graph::vertex_descriptor v){
	vector<Graph::vertex_descriptor> toreturn;
	while (g[v].is_root == false){
		toreturn.push_back(v);
		v = get_parent_node(v).first;
	}
	return toreturn;
}


set<Graph::vertex_descriptor> PopGraph::get_path_to_root(Graph::vertex_descriptor v){
	set<Graph::vertex_descriptor> toreturn;
	while (g[v].is_root == false){
		toreturn.insert(v);
		v = get_parent_node(v).first;
	}
	return toreturn;
}

pair<Graph::vertex_descriptor, Graph::vertex_descriptor> PopGraph::get_child_nodes(Graph::vertex_descriptor v){
	Graph::vertex_descriptor c1 = v;
	Graph::vertex_descriptor c2 = v;
	if (g[v].is_mig == true){
		cerr << "Calling get_child_nodes on a migration node\n";
		exit(1);
	}
	if (out_degree(v, g) != 2){
		return make_pair(c1, c2);
	}
	graph_traits<Graph>::out_edge_iterator out_i = out_edges(v, g).first;
	c1 = target(*out_i, g);
	while (g[c1].is_mig){
		graph_traits<Graph>::out_edge_iterator out_i2 = out_edges(c1, g).first;
		bool found = false;
		while (!found){
			if (g[*out_i2].is_mig == false) {
				c1 = target(*out_i2, g);
				found = true;
			}
			else out_i2++;
		}
	}
	out_i++;
	c2 = target(*out_i, g);
	while (g[c2].is_mig){
		graph_traits<Graph>::out_edge_iterator out_i2 = out_edges(c2, g).first;
		bool found = false;
		while (!found){
			if (g[*out_i2].is_mig == false) {
				c2 = target(*out_i2, g);
				found = true;
			}
			else out_i2++;
		}
	}
	return make_pair(c1, c2);
}

Graph::vertex_descriptor PopGraph::get_child_node_mig(Graph::vertex_descriptor v){
	Graph::vertex_descriptor toreturn = v;
	if (g[v].is_mig == false){
		cerr << "Calling get_child_node_mig on a non-migration node\n";
		exit(1);
	}
	if (out_degree(v, g) != 2){
		return toreturn;
	}
	while( g[toreturn].is_mig){
		graph_traits<Graph>::out_edge_iterator out_i = out_edges(v, g).first;
		if ( g[*out_i].is_mig == false) toreturn = target(*out_i, g);
		else{
			out_i++;
			toreturn = target(*out_i, g);
		}
	}
	return toreturn;
}

pair<Graph::vertex_descriptor, double> PopGraph::get_parent_node_wmig(Graph::vertex_descriptor v){
	pair<Graph::vertex_descriptor, double> toreturn;
	toreturn.first = v;
	toreturn.second = 0;
	graph_traits<Graph>::in_edge_iterator in_i = in_edges(v, g).first;
	while (in_i != in_edges(v, g).second){
		if (g[*in_i].is_mig == false) {
			toreturn.first = source(*in_i, g);
			toreturn.second = g[*in_i].len;
			return toreturn;
		}
		in_i++;
	}
	return toreturn;
}

pair<Graph::vertex_descriptor, double> PopGraph::get_parent_node(Graph::vertex_descriptor v){
	pair<Graph::vertex_descriptor, double> toreturn;
	if (g[v].is_root) {
		toreturn.first = v;
		toreturn.second = 0;
		return toreturn;
	}
	toreturn = get_parent_node_wmig(v);
	while (g[toreturn.first].is_mig == true) {
		pair<Graph::vertex_descriptor, double> tmp = get_parent_node_wmig(toreturn.first);
		toreturn.first = tmp.first;
		toreturn.second += tmp.second;
	}
	return toreturn;
}

bool PopGraph::does_mig_exist(Graph::vertex_descriptor st, Graph::vertex_descriptor sp){
	if (g[st].is_root) return false;
	Graph::vertex_descriptor v2 = get_parent_node_wmig(st).first;
	while (g[v2].is_mig == true){
		if ( edge(v2, sp, g).second == true) return true;
		v2 = get_parent_node_wmig(v2).first;
	}
	return false;
}

bool PopGraph::is_legal_migration(Graph::vertex_descriptor st, Graph::vertex_descriptor sp){
	if (g[st].is_mig || g[sp].is_mig) return false;
	if (st == sp) return false;
	if (g[st].is_root || g[sp].is_root) return false;
	if ( get_parent_node(st) == get_parent_node(sp)) return false;
	if ( does_mig_exist( st, sp)) return false;
	set<pair<double, set<Graph::vertex_descriptor> > > paths = get_paths_to_root(st);
	for (set<pair<double, set<Graph::vertex_descriptor> > >::iterator it = paths.begin(); it!= paths.end(); it++){
		if (it->second.find(sp) != it->second.end()) return false;
	}
	return true;
}

void PopGraph::remove_mig_edge(Graph::edge_descriptor e){
	if ( !g[e].is_mig){
		cerr << "ERROR: Calling remove_mig_edge on a non-migrations edge\n";
		exit(1);
	}

	//remove the edge
	Graph::vertex_descriptor v = source(e, g);
	remove_edge(e, g);

	//remove the node (ensure that everything remains linked up)

	Graph::vertex_descriptor parent, child;
	double len;
	Graph::in_edge_iterator init = in_edges(v, g).first;
	Graph::out_edge_iterator outit = out_edges(v, g).first;
	len = g[*init].len +g[*outit].len;

	parent = source(*init, g);
	child = target(*outit, g);

	Graph::edge_descriptor newe = add_edge(parent, child, g).first;

	g[newe].weight = 1;
	g[newe].len = len;
	g[newe].is_mig = false;
	clear_vertex(v, g);

	remove_vertex(v, g);
}
