/*
 * iterator.hpp
 *
 *  Created on: Mar 29, 2011
 *      Author: pickrell
 */

#ifndef ITERATOR_HPP_
#define ITERATOR_HPP_

#include "Node.h"
#include "NodeData.hpp"

namespace PhyloPop_Tree
{
        // Foreward declaration
        template<class NODEDATA>
        class BinaryTree;

//      template<class NODEDATA>
//      class BinaryNode;

        /// Iterator
        /**
         * Iterator which walks through the Nodes of the Tree.
         * The dereferenced Iterator points to the NodeData.
         */
        template<class NODEDATA>
        class iterator
        {

                friend class BinaryTree<NODEDATA>;

        public:
                /**
                 * The Copy Constructs copys the Pointer to a Node
                 * @param Pointer to Node
                 */
                iterator(const iterator& p_iterator);
                /**
                 * The Assignment Operator copys the Pointer to Node from the
                 * Argument
                 * @param iterator to copy
                 */
                iterator& operator=(const iterator& p_iterator);
                /**
                 * The Destructor does nothing else.
                 */
                ~iterator();
                /**
                 * The operator returns the dereferenced Pointer to Node.
                 * @return Reference to NodeData
                 */
                inline NODEDATA& operator*();
                /**
                 * The operator returns the Pointer to Node.
                 * @return Pointer to NodeData
                 */
                inline NODEDATA* operator->();
                /**
                 * False, if the Pointer to Node points to 0, else True
                 * @return bool
                 */
                inline operator bool();
                /**
                 * True, if the argument iterator points to the same location
                 * @return bool
                 */
                inline bool operator ==(const iterator& p_iterator) const;
                /**
                 * Compares the Pointer Locations
                 * @return bool
                 */
                inline bool operator <(const iterator& p_iterator) const;
                /**
                 * Returns an Iterator to the Father Node
                 * @return iterator
                 */
                inline iterator<NODEDATA> getFather();
                /**
                 * Returns an Iterator to the First Child Node
                 * @return iterator
                 */
                inline iterator<NODEDATA> getFirstChild();
                /**
                 * Returns an Iterator to the Next Child Node
                 * @return iterator
                 */
                inline iterator<NODEDATA> getNextChild();
                /**
                 * Returns an Iterator to the Last Child Node
                 * @return iterator
                 */
                inline iterator<NODEDATA> getLastChild();
                /**
                 * Sets the Father of the Iterator, to the Iterator given in the
                 * Argument
                 * @param iterator
                 */
                inline void setFather(const iterator<NODEDATA>& p_iterator);
                /**
                 * Sets the First Child of the Iterator, to the Iterator given in
                 * the Argument
                 * @param iterator&
                 */
                inline void setFirstChild(const iterator<NODEDATA>& p_iterator);
                /**
                 * Sets the Last Child of the Iterator, to the Iterator given in
                 * the Argument
                 * @param iterator&
                 */
                inline void setLastChild(const iterator<NODEDATA>& p_iterator);
                /**
                 * Appends the Iterator given in the Argument, as a Child Node
                 * if the Pointer to the Left Child is 0, the Iterator is inserte
                 *                         * else nothing happens.
                 * @param iterator&
                 */
                inline void appendChild(const iterator<NODEDATA>& p_iterator);
                /**
                 * Removes the Child Node given in the Argument from the Child Nodes
                 * @param iterator&
                 */
                inline void removeChild(const iterator<NODEDATA>& p_iterator);
                /**
                 * False, if the Pointer to Father points to 0, else True
                 * @return bool
                 */
                inline bool hasFather();
                /**
                 * False, if the Pointer to Right and Left Child points to 0,
                 * else True
                 * @return bool
                 */
                inline bool hasChildNodes();
                /**
                 * Creates a new allocated Node, with the Standard Constructor of
                 * Node
                 * @return iterator which points to the new allocated Node
                 */
                inline iterator<NODEDATA> createNewNode();
                /**
                 * Creates a new allocated Node, with the Copy Constructor of Node,
                 * which takes the Node where p_iterator points to as the prototype
                 * @return iterator which points to the new allocated Node
                 */
                inline iterator<NODEDATA> createNewNode(iterator<NODEDATA> p_iterator);
                /**
                 * Deletes a new allocated Node and frees the Memory via the Node
                 * Destructor
                 */
                inline void  deleteNode();

                /*
                 *  jiggles the height of the node from a U(currentheight-epsilon, currentheight+epsilon), reflecting excess
                 */
                //inline void jiggleheight(double, gsl_rng*);
        //private:
        	/**
        	 * The Constructor expects a Pointer to a Node
        	 * @param Pointer to Node
        	 */
        	iterator(Node<NODEDATA>* p_pNode = 0);

        private:
        	Node<NODEDATA>* m_pNode;

        };//iterator

}//TreeTime_MolecularTree

#include "iterator.cpp"
#endif /* ITERATOR_HPP_ */
