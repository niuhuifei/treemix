/*
 * GraphState3.cpp
 *
 *  Created on: Aug 12, 2011
 *      Author: pickrell
 */

#include "GraphState3.h"

GraphState3::GraphState3( GraphState2 * g, PhyloPop_params * p){
	state = g;
	params = p;
	scratch = new GraphState2( g->countdata, p);
}


void GraphState3::add_mig_targeted(){
	pair<string, string> target = state->get_max_resid();
	cout << "Targeting migration to "<< target.first << " "<< target.second << "\n"; cout.flush();

}
