/*
 * BinaryTree.h
 *
 *  Created on: Mar 29, 2011
 *      Author: pickrell
 */

#ifndef BINARYTREE_H_
#define BINARYTREE_H_

#include "Settings.hpp"
#include "iterator.h"
#include "Node.h"
#include "Tree.h"

namespace PhyloPop_Tree
{
	template<class NODEDATA>
	class BinaryTree: public Tree<NODEDATA>{
	public:
            /**
             * The Constructor of BinaryTree is build from the Newick String
             * given in the the first Argument,
             * which contains topological and branchlength Indformation.
             * The Names in the Newick String corresponding to a Nodes,
             * were mapped to a unique id, given by the the Map in the
             * second Argument, and assigned to the Nodes attribute.
             * @param Newick String
             * @param Map LeafIdentifier String to unique Id
             * @param Reference to MemoryAllocator
             */
            BinaryTree(const string& p_newickString,
                            map<string,int>& p_mPopToId);
            /**
             * The Copy Constructor of BinaryTree clones the whole
             * TreeStructure. The Data at the Nodes were copied
             * by the Standard CopyConstructor of NodeData.
             */
            BinaryTree(const BinaryTree& p_Tree);
            /**
             * The Destructor for BinaryTree recursively destructs all
             * NodeElements.
             */
            virtual ~BinaryTree();
            /**
             * Returns an iterator to the first Node (Root) in the BinaryTree.
             * @return Iterator to the Root Node
             */
            virtual iterator<NODEDATA> getRoot() const;
            /**
             * Sets the Root Node of the Tree to the Node in the Argument
             * Iterator.
             * @param Iterator to the new Root Node
             */
            virtual void setRoot(iterator<NODEDATA> p_rootIterator);
            /**
             * Clones the whole TreeStructure,
             * The Data at the Nodes were copyed by the Standard
             * Copy Constructor of NodeData.
             */
            virtual Tree<NODEDATA>* copy();
            /**
             *  print the tree structure in NewickFormat
             */
            //virtual void print_inorder(iterator<NODEDATA> p_rootIterator);
            /**
             * map of the population id to pointer to that node

           // inline map<int, iterator<NODEDATA> > get_tips(iterator<NODEDATA>);

           /*
            *  get lowest common ancestor of two nodes (this definitely works if the two are tips; must be modified otherwise)


            inline iterator<NODEDATA> get_LCA(iterator <NODEDATA>, iterator<NODEDATA>, iterator<NODEDATA>);

            /*
             *  takes an iterator, calculates distance to root


            inline double get_dist_to_root(iterator<NODEDATA>);
             */
    private:
            /**
             * The Assignment Operator is ambivalent,
             * so declared private and not implemented.
             */
            BinaryTree& operator=(const BinaryTree& p_Tree);
            /**
             * Helper Function for recursive copying and cloning the Tree
             * Structure
             */
            void copy(Node<NODEDATA>* p_pMaster, Node<NODEDATA>* p_pReplica);
            /**
             * Helper Function for recursive deleting the Tree Structure
             */
            void remove(Node<NODEDATA>* p_pNode);


    private:
            iterator<NODEDATA>      m_iRoot;/**<Iterator to the Root*/
	}; //BinaryTree


} //PhyloPop_Tree

#include "BinaryTree.cpp"
#endif /* BINARYTREE_H_ */
