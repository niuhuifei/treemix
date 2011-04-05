/*
 * Node.cpp
 *
 *  Created on: Mar 28, 2011
 *      Author: pickrell
 */
#pragma once
#include "Node.h"

namespace PhyloPop_Tree{

	template<class NODEDATA>
	Node<NODEDATA>::Node()
	: m_pFather(0)
	, m_pRightSon(0)
	, m_pLeftSon(0)
	, m_NodeData()
	{
		m_NodeData.m_id = -1;
		m_NodeData.m_len = -1;
		m_NodeData.m_time = -1;
	}

	template<class NODEDATA>
	Node<NODEDATA>::Node(const Node<NODEDATA >& p_Node)
    : m_pFather(0)
    , m_pRightSon(0)
    , m_pLeftSon(0)
    , m_NodeData(p_Node.m_NodeData){
	}

    template<class NODEDATA>
       Node<NODEDATA>::Node(const NODEDATA& p_NodeData)
       : m_pFather  (0)
       , m_pRightSon(0)
       , m_pLeftSon (0)
       , m_NodeData (p_NodeData){
       }

	template<class NODEDATA>
	Node<NODEDATA>& Node<NODEDATA>::operator=(const Node<NODEDATA>& p_Node){
		m_pFather   = 0;
		m_pRightSon = 0;
		m_pLeftSon  = 0;
		m_NodeData = p_Node.m_NodeData;
		return *this;
	}

	template<class NODEDATA>
	Node<NODEDATA>::~Node(){}

	template<class NODEDATA>
	inline Node<NODEDATA>* Node<NODEDATA>::getFather() const{
         return m_pFather;
	}

	template<class NODEDATA>
	inline Node<NODEDATA>* Node<NODEDATA>::getRightSon() const{
         return m_pRightSon;
	}

	template<class NODEDATA>
	inline Node<NODEDATA>* Node<NODEDATA>::getLeftSon() const{
         return m_pLeftSon;
}
	template<class NODEDATA>
	inline void Node<NODEDATA>::setFather(Node<NODEDATA>* p_pNewFather){
         m_pFather = p_pNewFather;
}
	template<class NODEDATA>
	inline void Node<NODEDATA>::setRightSon(Node<NODEDATA>* p_pNewRightSon){
         m_pRightSon = p_pNewRightSon;
}
	template<class NODEDATA>
	inline void Node<NODEDATA>::setLeftSon(Node<NODEDATA>* p_pNewLeftSon){
         m_pLeftSon = p_pNewLeftSon;
}
    template<class NODEDATA>
     inline NODEDATA& Node<NODEDATA>::getNodeData()
     {
             return m_NodeData;
     }



}//PhyloPop_Tree
