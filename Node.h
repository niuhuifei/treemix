/*
 * Node.h
 *
 *  Created on: Mar 28, 2011
 *      Author: pickrell
 */

#ifndef NODE_H_
#define NODE_H_

namespace PhyloPop_Tree
{
template<class NODEDATA>

	class Node{
	public:
		// contructor
		Node();

		// copy constructor
		// does not copy pointers
		Node (const Node&);

		// initializes pointers to 0, copies node data
		Node(const NODEDATA&);

		// assignment operator
		// does not copy pointers
		Node& operator=(const Node&);

		~Node();
		inline Node* getFather() const;
		inline Node* getRightSon() const;
		inline Node* getLeftSon() const;

		inline void setFather(Node*   p_pNewFather);
		inline void setRightSon(Node* p_pNewRightSon);
		inline void setLeftSon(Node*  p_pNewLeftSon);
		inline NODEDATA&        getNodeData();

	private:
		Node*     m_pFather;      //pointer to father
		Node*     m_pRightSon;    // pointer to right son
		Node*     m_pLeftSon;     //Pointer to the left son
		NODEDATA                m_NodeData;             /**< NodeData */
	}; //Node

} //PhyloPop_Tree

#include "Node.cpp"

#endif /* NODE_H_ */
