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

                         /*
                          * flip order of child nodes with probability 0.5
                          */
                         virtual inline void flip_sons(iterator<NODEDATA>, gsl_rng* );

                         /*
                          * make a vector of pointers to the nodes in order of traversal (input the number of populations)
                          */
                         virtual inline vector<iterator<NODEDATA> > get_inorder_traversal( int);

                         /*
                          * helper function does the inorder traversal
                          */

                         virtual inline void inorder_traverse(iterator<NODEDATA>, int*, vector<iterator<NODEDATA> >*);
                         /*
                          *  set heights of nodes
                          */
                         virtual inline void set_node_heights(vector<iterator<NODEDATA> >);
                         /*
                          * perturb node heights by adding a U(-epsilon, epsilon) to each (reflecting excess if necessary)
                          */
                         virtual inline void perturb_node_heights(vector<iterator<NODEDATA> >, double, gsl_rng*);
                         /*
                          *  reconstruct a tree from a traversal order and heights
                          */
                         virtual inline void build_tree(vector<iterator<NODEDATA> >);
                         /*
                          * recursive function to build trees
                          */
                         virtual inline void build_tree_helper(vector<iterator<NODEDATA> >*, int);
                         /*
                          * reset branch length from new tree
                          */
                         virtual inline void update_branch_lengths(iterator<NODEDATA>);

                         /*
                          * randomize a tree
                          */
                         virtual inline void randomize_tree(gsl_rng*);
                         /*
                          * get newick format
                          */
                         virtual inline string get_newick_format();
                         virtual inline void newick_helper(iterator<NODEDATA>, string*);
        };//Tree

}//PhyloPop_Tree

#include "Tree.cpp"
#endif /* TREE_HPP_ */
