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
		indexcounter++;

		g[vintnew].index = indexcounter;
		g[vintnew].name = "NA";
		g[vintnew].height = 0;
		g[vintnew].is_tip = false;
		g[vintnew].is_root = false;
		g[vintnew].rev = false;
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
		indexcounter++;

		g[vintnew].index = indexcounter;
		g[vintnew].name = "NA";
		g[vintnew].height = 0;
		g[vintnew].is_tip = false;
		g[vintnew].is_root = true;
		g[vintnew].rev = false;
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
	// from     p1   \
	//          /     \
	//         /       sp
	//        st
	// to
	//
	//
    //
	//
	//       p1     /
	//      /      /
	//     /      /
	//    st-----sp
 	//
	//

	Graph::edge_descriptor e;
	istree = false;
	e = add_edge(st, sp, g).first;
	g[e].weight = 0.25;
	g[e].len = 1;
	g[e].is_mig = true;
	return e;

	/*
	Graph::vertex_descriptor p1, p2;
	Graph::edge_descriptor e;
	Graph::in_edge_iterator in_it = in_edges(st, g).first;

	//find p1
	bool foundparent = false;
	while (!foundparent){
		if ( g[*in_it].is_mig == false){
			foundparent = true;
			p1 = source(*in_it, g);
			e = *in_it;
		}
		in_it++;
	}
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
	indexcounter++;

	//add new edges
	e = add_edge(p1, p2, g).first;
	g[e].weight = 1;
	g[e].len = oldlen/2;
	g[e].is_mig = false;

	e = add_edge(p2, st, g).first;
	g[e].weight = 1;
	g[e].len = oldlen/2;
	g[e].is_mig = false;

	e = add_edge(p2, sp, g).first;
	g[e].weight = 0.25;
	g[e].len = 1;
	g[e].is_mig = true;

	return p2;
	*/
}

set<pair<double, set<Graph::vertex_descriptor> > > PopGraph::get_paths_to_root(Graph::vertex_descriptor v){

	set<pair<double, set<Graph::vertex_descriptor> > > toreturn;
	Graph::in_edge_iterator init = in_edges(v, g).first;
	while (init != in_edges(v, g).second){

		double w = g[*init].weight;
		//scout << "here2 " << w << " "<< g[v].index << "\n"; cout.flush();
		set<Graph::vertex_descriptor> tmpset;
		tmpset.insert(v);
		Graph::vertex_descriptor parent = source(*init, g);
		if (g[parent].is_root){
			tmpset.insert(parent);
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
	//cout << "done\n";
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


set<Graph::vertex_descriptor> PopGraph::get_root_adj(){
	set<Graph::vertex_descriptor> toreturn;
	Graph::out_edge_iterator it = out_edges(root, g).first;
	toreturn.insert( target(*it, g));
	it++;
	toreturn.insert( target(*it, g));
	return toreturn;
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

vector<Graph::vertex_descriptor > PopGraph::get_inorder_traversal_noroot(int nodes){
	vector<Graph::vertex_descriptor> toreturn;
	vector<Graph::vertex_descriptor> tmp = get_inorder_traversal(nodes);
	for (int i = 0; i < tmp.size(); i++){
		if (g[tmp[i]].is_root == false) toreturn.push_back(tmp[i]);
	}
	return toreturn;
}
void PopGraph::inorder_traverse(Graph::vertex_descriptor rootIterator, int* i, vector<Graph::vertex_descriptor >* v){
	graph_traits<Graph>::out_edge_iterator out_i = out_edges(rootIterator, g).first;
	//cout << g[rootIterator].index << " "<< g[rootIterator].name << " "<< g[rootIterator].is_tip << " "<< g[rootIterator].rev << "\n";
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
		//cout << "in get_LCA "<< g[root_v].index << " "<< g[tip1_v].index << " "<< g[tip2_v].index << "\n"; cout.flush();
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

void PopGraph::local_update(Graph::vertex_descriptor v3, gsl_rng* r){
	/*
	 *
	 * go from    1
	 *          /  \*
	 *         /    3
	 *        /    / \
	 *       2    4   5
	 *
	 * do local rearrangements.
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
	//cout << min << " "<< l2 << " "<< l4 << " "<< l5 << "\n";
	//choose which rearrangement to do
	// either 4 or 5 slides with 3
	/*double ran = gsl_rng_uniform(r);
	ran = 0.7;
	if (ran < 0.5){
		// 4 slides
		double ran2 = gsl_rng_uniform(r);
		double newpos = ran2*min;
		remove_edge(v1, v2, g);
		remove_edge(v3, v5, g);
		Graph::edge_descriptor e = add_edge(v3, v2, g).first;
		g[e].weight = 1;
		g[e].len = min - newpos;
		e = add_edge(v1, v5, g).first;
		g[e].weight = 1;
		g[e].len = l5;
		e = edge(v1, v3, g).first;
		g[e].len = newpos;
		e = edge(v3, v4, g).first;
		g[e].len = l4 - newpos;
		g[v3].height = g[v1].height+ newpos;
	}
	else{
		//5 slides
		double ran2 = gsl_rng_uniform(r);
		double newpos = ran2*min;
		remove_edge(v1, v2, g);
		remove_edge(v3, v4, g);
		Graph::edge_descriptor e = add_edge(v3, v2, g).first;
		g[e].weight = 1;
		g[e].len = min - newpos;
		e = add_edge(v1, v4, g).first;
		g[e].weight = 1;
		g[e].len = l4;
		e = edge(v1, v3, g).first;
		g[e].len = newpos;
		e = edge(v3, v5, g).first;
		g[e].len = l5 - newpos;
		g[v3].height = g[v1].height+ newpos;
	}
*/
	//choose which rearrangement to do
	double ran = gsl_rng_uniform(r);
	ran = 0.2;
	if (ran < 0.5){
		// go to ((2,4),5)

		remove_edge(v1, v2, g);
		remove_edge(v3, v5, g);
		Graph::edge_descriptor e = add_edge(v3, v2, g).first;
		g[e].weight = 1;
		g[e].len = d12;
		e = add_edge(v1, v5, g).first;
		g[e].weight = 1;
		g[e].len = d35;
		update_heights_local( v5, -d13);
		update_heights_local( v2, d13);

	}
	else{
		//go to (4,(2,5))

		remove_edge(v1, v2, g);
		remove_edge(v3, v4, g);
		Graph::edge_descriptor e = add_edge(v3, v2, g).first;
		g[e].weight = 1;
		g[e].len = d12;

		e = add_edge(v1, v4, g).first;
		g[e].weight = 1;
		g[e].len = d34;
		update_heights_local( v4, -d13);
		update_heights_local( v2, d13);
	}
	}

}


void PopGraph::update_heights_local(Graph::vertex_descriptor v, double adj){
	if (out_degree(v, g) > 0){
		Graph::out_edge_iterator it = out_edges(v, g).first;
		update_heights_local( target(*it, g), adj);
		it++;
		update_heights_local( target(*it, g), adj);
	}
	g[v].height += adj;
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
		e = add_edge(v1, v5, g).first;
		g[e].weight = 1;
		g[e].len = d35;

	}
	else{
		//go to (4,(2,5))

		remove_edge(v1, v2, g);
		remove_edge(v3, v4, g);
		Graph::edge_descriptor e = add_edge(v3, v2, g).first;
		g[e].weight = 1;
		g[e].len = d12;

		e = add_edge(v1, v4, g).first;
		g[e].weight = 1;
		g[e].len = d34;

	}
	}

}


void PopGraph::local_update_branches(Graph::vertex_descriptor v3, gsl_rng* r, double epsilon){
	/*
	 *
	 * go from    1
	 *          /  \*
	 *         /    3
	 *        /    / \
	 *       2    4   5
	 *
	 * do updates of branches in vicinity. input is v3
	 *
	 */
	if (!g[v3].is_root){

		// set up descriptors and edge lengths
		Graph::in_edge_iterator tmp = in_edges(v3, g).first;
		Graph::edge_descriptor ed = *tmp;
		Graph::edge_descriptor e12, e34, e35;
		Graph::vertex_descriptor v1, v2, v4, v5;
		double d12, d13, d34, d35;

		d13 = g[ed].len;
		v1 = source(ed, g);
		Graph::out_edge_iterator out3 = out_edges(v3, g).first;
		v4 = target(*out3, g);
		d34 = g[*out3].len;
		e34 = *out3;

		++out3;
		v5 = target(*out3, g);
		d35 = g[*out3].len;
		e35 = *out3;
		Graph::out_edge_iterator out1 = out_edges(v1, g).first;
		if (target(*out1, g) == v3) {
			out1++;
			v2 = target(*out1, g);
			d12 = g[*out1].len;
			e12 = *out1;
		}
		else {
			v2 = target(*out1, g);
			d12 = g[*out1].len;
			e12 = *out1;
		}

		single_branch_update( e12, r, epsilon);
		single_branch_update( ed, r, epsilon);
		single_branch_update( e34, r, epsilon);
		single_branch_update( e35, r, epsilon);
	}
}

void PopGraph::single_branch_update(Graph::edge_descriptor e, gsl_rng* r, double epsilon){
	double toadd = (2* gsl_rng_uniform(r) - 1)*epsilon;
	double oldlen = g[e].len;
	double newlen = fabs(oldlen+toadd);
	g[e].len = newlen;
	update_heights_local( target(e, g), newlen-oldlen);
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
    		g[trav[i]].height = newheight; //removed if (newheight > 0), this now allows re-rooting
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
 	double tipheight = gsl_rng_uniform(r)*0.01+0.01;

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
	for (vp = vertices(g); vp.first != vp.second; ++vp.first){
		std::cout << index[*vp.first] <<  ":"<< g[*vp.first].name << ":"<< get_dist_to_root(*vp.first);
		if (g[*vp.first].is_root) std::cout << ":ROOT";
		if (g[*vp.first].is_tip) std::cout << ":TIP";
		cout << " ";

	}
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


set<Graph::vertex_descriptor> PopGraph::get_path_to_root(Graph::vertex_descriptor v){
	set<Graph::vertex_descriptor> toreturn;
	while (g[v].is_root == false){
		toreturn.insert(v);
		graph_traits<Graph>::in_edge_iterator in_i = in_edges(v, g).first;
		v = source(*in_i, g);
	}
	return toreturn;
}
