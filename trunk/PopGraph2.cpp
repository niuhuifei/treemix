/*
 * PopGraph2.cpp
 *
 *  Created on: Jun 17, 2011
 *      Author: pickrell
 */

#include "PopGraph2.h"

PopGraph2::PopGraph2(string meanpop, vector<string> first3pops){

	g = Graph(1);
	istree = true;
	Graph::vertex_descriptor v = *vertices(g).first;
	Graph::vertex_descriptor v2, v3, v4, v5, v6;
	Graph::edge_descriptor e;
	g[v].index = 0;
	g[v].name = meanpop;
	g[v].height = 0;
	g[v].is_tip = false;
	g[v].is_root = true;
	g[v].rev = false;
	root = v;


	v2 = add_vertex(g);
	g[v2].index = 1;
	g[v2].name = "NA";
	g[v2].height = 0;
	g[v2].is_tip = false;
	g[v2].is_root = false;
	g[v2].rev = false;
	e = add_edge( v, v2, g).first;
	g[e].weight = 1;
	g[e].len = 1;

	// add the first three populations to the psuedo-tree
	// tree will look like this:
	//         meanpop
	//            |
	//        pseudoroot
	//            /\
	//           /  \
	//          /   /\
	//         1    2 3
	if (first3pops.size() != 3) {
		std::cout << "Error initializing PopGraph: need 3 pops, not "<< first3pops.size() << "\n";
		exit(1);
	}
	v3 = add_vertex(g);
	g[v3].index = 2;
	g[v3].name = first3pops[0];
	g[v3].height = 0;
	g[v3].is_tip = true;
	g[v3].is_root = false;
	g[v3].rev = false;
	e = add_edge( v2, v3, g).first;
	g[e].weight = 1;
	g[e].len = 1;


	v4 = add_vertex(g);
	g[v4].index = 3;
	g[v4].name = "NA";
	g[v4].height = 0;
	g[v4].is_tip = true;
	g[v4].is_root = false;
	g[v4].rev = false;
	e = add_edge( v2, v4, g).first;
	g[e].weight = 1;
	g[e].len = 1;


	v5 = add_vertex(g);
	g[v5].index = 4;
	g[v5].name = first3pops[1];
	g[v5].height = 0;
	g[v5].is_tip = true;
	g[v5].is_root = false;
	g[v5].rev = false;
	e = add_edge( v4, v5, g).first;
	g[e].weight = 1;
	g[e].len = 1;

	v6 = add_vertex(g);
	g[v6].index = 5;
	g[v6].name = first3pops[2];
	g[v6].height = 0;
	g[v6].is_tip = true;
	g[v6].is_root = false;
	g[v6].rev = false;
	e = add_edge( v4, v6, g).first;
	g[e].weight = 1;
	g[e].len = 1;
	indexcounter = 6;
}


void PopGraph2::copy(PopGraph2 * s){
	istree = s->istree;
	popnames = s->popnames;

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
	}
}



void PopGraph2::print(){
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





double PopGraph2::get_dist_to_root(Graph::vertex_descriptor v){
	double toreturn = 0;
	while (g[v].is_root == false){
		graph_traits<Graph>::in_edge_iterator in_i = in_edges(v, g).first;
		toreturn += g[*in_i].len;
		v = source(*in_i, g);
	}
	return toreturn;
}

Graph::vertex_descriptor PopGraph2::add_tip(Graph::vertex_descriptor v, string name){
	// add a tip attached above v, halfway to the parent node
	Graph::vertex_descriptor vtipnew, vintnew, vparent;
	double len, newlen;
	Graph::edge_descriptor e;

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
	e = add_edge( vintnew, v, g).first;
	g[e].weight = 1;
	g[e].len = newlen;
	e = add_edge( vintnew, vtipnew, g).first;
	g[e].weight = 1;
	g[e].len = 1;
	return vtipnew;
}


void PopGraph2::remove_tip(Graph::vertex_descriptor v){
	Graph::vertex_descriptor vparent, vparent2, vnewdest;
	double newlen;
	Graph::edge_descriptor e;

	graph_traits<Graph>::in_edge_iterator init = in_edges(v, g).first;
	vparent = source(*init, g);

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



vector<Graph::vertex_descriptor > PopGraph2::get_inorder_traversal(int npop){
        	vector<Graph::vertex_descriptor> toreturn(2*npop);
        	graph_traits<Graph>::out_edge_iterator outit = out_edges(root, g).first;
        	Graph::vertex_descriptor pseudoroot = target(*outit, g);
        	toreturn[0] = root;

        	int count = 1;
        	//cout << g[pseudoroot].index << "\n";
        	inorder_traverse(pseudoroot, &count, &toreturn);
        	return toreturn;
}

void PopGraph2::inorder_traverse(Graph::vertex_descriptor rootIterator, int* i, vector<Graph::vertex_descriptor >* v){
	graph_traits<Graph>::out_edge_iterator out_i = out_edges(rootIterator, g).first;

	if (out_degree(rootIterator, g) == 2) inorder_traverse(target(*out_i, g), i , v);

 	v->at(*i) = rootIterator;
 	*i = *i+1;

 	if (out_degree(rootIterator, g) == 2){
 		++out_i;
 		inorder_traverse(target(*out_i, g), i, v);
 	}

}


vector<Graph::vertex_descriptor > PopGraph2::get_inorder_traversal_noroot(int nodes){
	vector<Graph::vertex_descriptor> toreturn;
	vector<Graph::vertex_descriptor> tmp = get_inorder_traversal(nodes);
	for (int i = 0; i < tmp.size(); i++){
		if (g[tmp[i]].is_root == false) toreturn.push_back(tmp[i]);
	}
	return toreturn;
}



Graph::vertex_descriptor PopGraph2::get_LCA(Graph::vertex_descriptor root_v,
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


set<Graph::vertex_descriptor> PopGraph2::get_path_to_root(Graph::vertex_descriptor v){
	set<Graph::vertex_descriptor> toreturn;
	while (g[v].is_root == false){
		toreturn.insert(v);
		graph_traits<Graph>::in_edge_iterator in_i = in_edges(v, g).first;
		v = source(*in_i, g);
	}
	return toreturn;
}
