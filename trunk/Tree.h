/*
 * Tree.hpp
 *
 *  Created on: Mar 29, 2011
 *      Author: pickrell
 */

#ifndef TREE_HPP_
#define TREE_HPP_


#include "Settings.hpp"
#include "iterator.h"
#include "NodeData.hpp"

namespace PhyloPop_Tree
{

//      class iterator;// Forward declaration

        /**
         * @brief Tree Interface.
         *
         * Tree Interface.
         * Access to the NodeData via the Iterator iterator.
         */
        template<class NODEDATA>
        class Tree
        {
                public:
                        /**
                         * Copys a Tree
                         * @return Pointer to the new allocated Tree
                         */
                        virtual Tree* copy() = 0;
                        /**
                         * Returns an Iterator to the firts Element (Root) of the Tree
                         * @return Iterator to the first Element
                         */
                        virtual iterator<NODEDATA> getRoot() const = 0;
                        /**
                         * Sets the Root Node of the Tree, given in the Argument.
                         * @param Iterator to the new Root Node
                         */
                        virtual void setRoot(iterator<NODEDATA> p_rootIterator)=0;

                        virtual ~Tree();
                        /**
                          *  print the tree structure in NewickFormat
                          */
                         virtual void print_inorder(iterator<NODEDATA> p_rootIterator);
                         /**
                          * map of the population id to pointer to that node
                          */
                        virtual inline map<int, iterator<NODEDATA> > get_tips(iterator<NODEDATA>);

                        /*
                         *  get lowest common ancestor of two nodes (this definitely works if the two are tips; must be modified otherwise)
                         */

                         virtual inline iterator<NODEDATA> get_LCA(iterator <NODEDATA>, iterator<NODEDATA>, iterator<NODEDATA>);

                         /*
                          *  takes an iterator, calculates distance to root
                          */

                         virtual inline double get_dist_to_root(iterator<NODEDATA>);


        };//Tree

}//PhyloPop_Tree

#include "Tree.cpp"
#endif /* TREE_HPP_ */
