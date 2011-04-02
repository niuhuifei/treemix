/*
Copyright (C) 2006-2008 Lin Himmelmann

This file is part of TreeTime.
See the NOTICE file distributed with this work for additional
information regarding copyright ownership and licensing.

TreeTime is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

TreeTime is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with TreeTime.  If not, see <http://www.gnu.org/licenses/>.
*/

#pragma once
#include "iterator.h"
namespace PhyloPop_Tree
{

	template<class NODEDATA>
	iterator<NODEDATA>::iterator(Node<NODEDATA>* p_pNode)
	: m_pNode(p_pNode)
	{
	}

	template<class NODEDATA>
	iterator<NODEDATA>::iterator(const iterator<NODEDATA>& p_iterator)
	: m_pNode(p_iterator.m_pNode)
	{
	}

	template<class NODEDATA>
	iterator<NODEDATA>& iterator<NODEDATA>::operator=(const iterator<NODEDATA>& p_iterator)
	{
		m_pNode = p_iterator.m_pNode;
		return *this;
	}

	template<class NODEDATA>
	iterator<NODEDATA>::~iterator()
	{
	}

	template<class NODEDATA>
	inline NODEDATA& iterator<NODEDATA>::operator*()
	{
		return m_pNode->getNodeData();
	}

	template<class NODEDATA>
	inline NODEDATA* iterator<NODEDATA>::operator->()
	{
		return &(m_pNode->getNodeData());
	}

	template<class NODEDATA>
	inline iterator<NODEDATA>::operator bool()
	{
		if(m_pNode)
			return true;
		return false;
	}

	template<class NODEDATA>
	inline bool iterator<NODEDATA>::operator==(const iterator<NODEDATA>& p_iterator) const
	{
		return m_pNode == p_iterator.m_pNode;
	}

	template<class NODEDATA>
	inline bool iterator<NODEDATA>::operator<(const iterator<NODEDATA>& p_iterator) const
	{
		return (m_pNode < p_iterator.m_pNode);
	}

	template<class NODEDATA>
	inline iterator<NODEDATA> iterator<NODEDATA>::getFather()
	{
		return iterator<NODEDATA>(m_pNode->getFather());
	}

	template<class NODEDATA>
	inline iterator<NODEDATA> iterator<NODEDATA>::getFirstChild()
	{
		return iterator<NODEDATA>(m_pNode->getLeftSon());
	}

	template<class NODEDATA>
	inline iterator<NODEDATA> iterator<NODEDATA>::getNextChild()
	{
		if(m_pNode->getFather()->getLeftSon() == m_pNode)
			return iterator<NODEDATA>(m_pNode->getRightSon());
		return iterator<NODEDATA>();
	}

	template<class NODEDATA>
	inline iterator<NODEDATA> iterator<NODEDATA>::getLastChild()
	{
		return iterator<NODEDATA>(m_pNode->getRightSon());
	}

	template<class NODEDATA>
	inline void iterator<NODEDATA>::setFather(const iterator<NODEDATA>& p_iterator)
	{
		m_pNode->setFather(p_iterator.m_pNode);
	}

	template<class NODEDATA>
	inline void iterator<NODEDATA>::setFirstChild(const iterator<NODEDATA>& p_iterator)
	{
		m_pNode->setLeftSon(p_iterator.m_pNode);
	}

	template<class NODEDATA>
	inline void iterator<NODEDATA>::setLastChild(const iterator<NODEDATA>& p_iterator)
	{
		m_pNode->setRightSon(p_iterator.m_pNode);
	}

	template<class NODEDATA>
	inline void iterator<NODEDATA>::appendChild(const iterator<NODEDATA>& p_iterator)
	{
		if(m_pNode->getLeftSon() == 0)
			m_pNode->setLeftSon(p_iterator.m_pNode);
		else if(m_pNode->getRightSon()  == 0)
			m_pNode->setRightSon(p_iterator.m_pNode);
	}

	template<class NODEDATA>
	inline void iterator<NODEDATA>::removeChild(const iterator<NODEDATA>& p_iterator)
	{
		if(m_pNode->getLeftSon() == p_iterator.m_pNode)
			m_pNode->setLeftSon(0);
		if(m_pNode->getRightSon() == p_iterator.m_pNode)
			m_pNode->setRightSon(0);
	}

	template<class NODEDATA>
	inline bool iterator<NODEDATA>::hasFather()
	{
		if(m_pNode->getFather() != 0)
			return true;
		return false;
	}

	template<class NODEDATA>
	inline bool iterator<NODEDATA>::hasChildNodes()
	{
		if(m_pNode->getLeftSon() == 0 && m_pNode->getRightSon() == 0)
			return false;
		return true;
	}

	template<class NODEDATA>
	 inline iterator<NODEDATA> iterator<NODEDATA>::createNewNode()
	 {
	 	Node<NODEDATA>* p = new Node<NODEDATA>;
		return iterator<NODEDATA>(p);
	 }

	template<class NODEDATA>
	 inline iterator<NODEDATA> iterator<NODEDATA>::createNewNode(iterator<NODEDATA> p_iterator)
	 {
		Node<NODEDATA>* p = new Node<NODEDATA>(*p_iterator);
		return iterator<NODEDATA>(p);
	 }

	template<class NODEDATA>
	 inline void iterator<NODEDATA>::deleteNode()
	 {
	 	delete m_pNode;
	 }

}//PhyloPop_Tree
