/*
 * GraphState3.h
 *
 *  Created on: Aug 12, 2011
 *      Author: pickrell
 */

#ifndef GRAPHSTATE3_H_
#define GRAPHSTATE3_H_

#include "GraphState2.h"

class GraphState3{
public:
	GraphState3(GraphState2*, PhyloPop_params *);

	GraphState2 * state;
	GraphState2 * scratch;
	PhyloPop_params * params;

	bool add_mig_targeted();
};

#endif /* GRAPHSTATE3_H_ */
