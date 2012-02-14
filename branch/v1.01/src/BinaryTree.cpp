/*
 * BinaryTree.cpp
 *
 *  Created on: Mar 29, 2011
 *      Author: pickrell
 */

#pragma once
#include "BinaryTree.h"

namespace PhyloPop_Tree
{
        template<class NODEDATA>
        BinaryTree<NODEDATA>::BinaryTree(const string& p_newickString,
                        map<string,int>& p_mPopToId)
        {
                Node<NODEDATA>* pRootNode  = new Node<NODEDATA>;
                m_iRoot = iterator<NODEDATA>(pRootNode);
                Node<NODEDATA>* aZ = pRootNode;//Initialise Pointer to the Root
                for(string::const_iterator I = p_newickString.begin();
                                I != p_newickString.end(); ++I)
                {
                        if ( *I == '(' )//Begin new subtree
                        {
                                Node<NODEDATA>* p = new Node<NODEDATA>;
                                p->setFather(aZ);
                                if(aZ->getLeftSon() == 0)
                                        aZ->setLeftSon(p);
                                else
                                        aZ->setRightSon(p);
                                //aZ->setLeftSon(p);
                                aZ = p;
                        }
                else if( *I == ')' )// Subtree finished, get back
                {
                        aZ = aZ->getFather();
                }
                else if( *I == ',' )// insert brother
                {
                        Node<NODEDATA>* p = new Node<NODEDATA>;
                        aZ = aZ->getFather();
                        p->setFather(aZ);
                        if(aZ->getLeftSon() == 0)
                                aZ->setLeftSon(p);
                        else
                                aZ->setRightSon(p);
                        //aZ->setRightSon(p);
                        aZ = p;
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
                        aZ->getNodeData().m_len = atof(length.c_str());
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
                        aZ->getNodeData().m_id = p_mPopToId[name];
                }
               }
        }


        template<class NODEDATA>
        BinaryTree<NODEDATA>::BinaryTree(const BinaryTree<NODEDATA>& p_BinaryTree)
        {
                iterator<NODEDATA> iMaster = p_BinaryTree.getRoot();
                if(iMaster)
                {
                        Node<NODEDATA>* pRootNode  = new Node<NODEDATA>(*iMaster);
                        m_iRoot = iterator<NODEDATA>(pRootNode);
                        copy(iMaster.m_pNode,pRootNode);
                }
        }

        template<class NODEDATA>
        void BinaryTree<NODEDATA>::copy(Node<NODEDATA>* p_pMaster, Node<NODEDATA>* p_pReplica)
        {
                Node<NODEDATA>* pMasterChild = p_pMaster->getLeftSon();
                if(pMasterChild != 0)
                {
                        Node<NODEDATA>* pReplicaChild = new Node<NODEDATA>(*pMasterChild);
                        pReplicaChild->setFather(p_pReplica);
                        p_pReplica->setLeftSon(pReplicaChild);
                        copy(pMasterChild,pReplicaChild);
                }
                pMasterChild = p_pMaster->getRightSon();
                if(pMasterChild != 0)
                {
                        Node<NODEDATA>* pReplicaChild = new Node<NODEDATA>(*pMasterChild);
                        pReplicaChild->setFather(p_pReplica);
                        p_pReplica->setRightSon(pReplicaChild);
                        copy(pMasterChild,pReplicaChild);
                }
        }

        template<class NODEDATA>
        void BinaryTree<NODEDATA>::remove(Node<NODEDATA>* p_pNode)
        {
                Node<NODEDATA>* pChildNode = p_pNode->getLeftSon();
                if(pChildNode != 0)
                {
                        remove(pChildNode);
                }
                pChildNode = p_pNode->getRightSon();
                if(pChildNode != 0)
                {
                        remove(pChildNode);
                }
                delete p_pNode;
        }

        template<class NODEDATA>
        BinaryTree<NODEDATA>::~BinaryTree()
        {
                if(m_iRoot)
                        remove(m_iRoot.m_pNode);
        }

        template<class NODEDATA>
        iterator<NODEDATA> BinaryTree<NODEDATA>::BinaryTree::getRoot() const
        {
                return m_iRoot;
        }

        template<class NODEDATA>
        void BinaryTree<NODEDATA>::setRoot(iterator<NODEDATA> p_rootIterator)
        {
                m_iRoot = p_rootIterator;
        }
        template<class NODEDATA>
         Tree<NODEDATA>* BinaryTree<NODEDATA>::copy()
         {
                 return new BinaryTree<NODEDATA>(*this);
         }
}
/*        template<class NODEDATA>
        void BinaryTree<NODEDATA>::print_inorder(iterator<NODEDATA> p_rootIterator){
        	if (p_rootIterator != NULL){
        		print_inorder(p_rootIterator.getFirstChild());
        		cout << p_rootIterator->m_id <<" "<< p_rootIterator->m_len<< "\n";
        		print_inorder(p_rootIterator.getLastChild());
        	}

        }


        template<class NODEDATA>
        map<int, iterator<NODEDATA> > BinaryTree<NODEDATA>::get_tips(iterator<NODEDATA> p_rootIterator){
        	map<int, iterator<NODEDATA> > toreturn;
        	if (p_rootIterator.hasChildNodes() == false){
        		toreturn.insert(make_pair(p_rootIterator->m_id, p_rootIterator));
        	}
        	else{
         		map<int, iterator<NODEDATA> > t1 = get_tips(p_rootIterator.getFirstChild());
          		map<int, iterator<NODEDATA> > t2 = get_tips(p_rootIterator.getLastChild());
          		for (typename map<int, iterator<NODEDATA> >::const_iterator it1 = t1.begin(); it1 != t1.end(); it1++){
          			//pair<int, iterator<NODEDATA> > = std::make_pair(it1->first, it1->second);
          			toreturn.insert(make_pair(it1->first, it1->second));
          		}
        		for (typename map<int, iterator<NODEDATA> >::iterator it2 = t2.begin(); it2 != t2.end(); it2++){
          			toreturn.insert(std::make_pair(it2->first, it2->second));
          		}
         	}
        	return toreturn;
         }
/*
        template<class NODEDATA>
        iterator<NODEDATA> BinaryTree<NODEDATA>::get_LCA(iterator<NODEDATA> p_rootIterator,
        		iterator<NODEDATA> p_tip1Iterator, iterator<NODEDATA> p_tip2Iterator){
			if (!p_rootIterator) return NULL;
			if (p_rootIterator.getFirstChild() == p_tip1Iterator || p_rootIterator.getFirstChild() == p_tip2Iterator
					|| p_rootIterator.getLastChild() ==p_tip2Iterator  || p_rootIterator.getLastChild() ==p_tip1Iterator){
				return p_rootIterator;
			}
			else{
				iterator<NODEDATA> firstit = get_LCA(p_rootIterator.getFirstChild(), p_tip1Iterator, p_tip2Iterator);
				iterator<NODEDATA> lastit = get_LCA(p_rootIterator.getLastChild(), p_tip1Iterator, p_tip2Iterator);
				if (firstit && lastit) return p_rootIterator;
				else if (firstit) return firstit;
				else return lastit;
			}

        }

        template<class NODEDATA>
        double BinaryTree<NODEDATA>::get_dist_to_root(iterator<NODEDATA> input_it){
        	double toreturn = 0;
        	while (input_it.hasFather()){
        		toreturn += input_it->m_len;
        		input_it = input_it.getFather();
        	}
        	return toreturn;
        }
}
*/
