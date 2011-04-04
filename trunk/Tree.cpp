/*
 * Tree.cpp
 *
 *  Created on: Mar 29, 2011
 *      Author: pickrell
 */
#pragma once
#include "Tree.h"

namespace PhyloPop_Tree
{
        template<class NODEDATA>
        Tree<NODEDATA>::~Tree()
        {
        }

        template<class NODEDATA>
        void Tree<NODEDATA>::print_inorder(iterator<NODEDATA> p_rootIterator){
        	if (p_rootIterator != NULL){
        		print_inorder(p_rootIterator.getFirstChild());
        		cout << p_rootIterator->m_id <<" "<< p_rootIterator->m_len<< "\n";
        		print_inorder(p_rootIterator.getLastChild());
        	}

        }


        template<class NODEDATA>
        map<int, iterator<NODEDATA> > Tree<NODEDATA>::get_tips(iterator<NODEDATA> p_rootIterator){
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
        			toreturn.insert(make_pair(it2->first, it2->second));
        		}
        	}
        	return toreturn;
        }

        template<class NODEDATA>
        iterator<NODEDATA> Tree<NODEDATA>::get_LCA(iterator<NODEDATA> p_rootIterator,
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
        double Tree<NODEDATA>::get_dist_to_root(iterator<NODEDATA> input_it){
        	double toreturn = 0;
        	while (input_it.hasFather()){
        		toreturn += input_it->m_len;
        		input_it = input_it.getFather();
        	}
        	return toreturn;
        }
        template<class NODEDATA>
         void Tree<NODEDATA>::flip_sons(iterator<NODEDATA> root_it, gsl_rng* r){
         	if (root_it.hasChildNodes()== true){
         		//cout << "here\n"; cout.flush();
         		double ran = gsl_rng_uniform(r);
         		//cout << ran <<"\n"; cout.flush();
         		if (ran <0.5){
         			iterator<NODEDATA> tmp = root_it.getFirstChild();
         			root_it.setFirstChild(root_it.getLastChild());
         			root_it.setLastChild(tmp);
         		}
         		flip_sons(root_it.getFirstChild(), r);
         		flip_sons(root_it.getLastChild(), r);
         	}
         }
}

