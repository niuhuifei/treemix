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
        		//cout << input_it->m_id << " "<< input_it->m_len << "\n";
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
        template< class NODEDATA>
        vector<iterator<NODEDATA> > Tree<NODEDATA>::get_inorder_traversal(int nodes){
        	vector<iterator<NODEDATA> > toreturn(2*nodes-1);
        	int count = 0;
        	inorder_traverse(getRoot(), &count, &toreturn);
        	return toreturn;
        }

        template<class NODEDATA>
        void Tree<NODEDATA>::inorder_traverse(iterator<NODEDATA> rootIterator, int* i, vector<iterator<NODEDATA> >* v){
        	if (rootIterator != NULL){
        		//cout << rootIterator->m_id <<" "<< rootIterator->m_len<<  " "<< *i << "\n"; cout.flush();
        		inorder_traverse(rootIterator.getFirstChild(),i , v);
        		v->at(*i) = rootIterator;

        		//cout << "before " << *i << "\n";
        		*i = *i+1;
        		//cout << *i << "\n"; cout.flush();
        		inorder_traverse(rootIterator.getLastChild(), i, v);
        	}
        }
        template<class NODEDATA>
        void Tree<NODEDATA>::set_node_heights(vector<iterator<NODEDATA> > trav){
        	for (int i = 0; i < trav.size(); i++){
        		trav[i]->m_time = get_dist_to_root(trav[i]);
        	}
        }
        template<class NODEDATA>
        void Tree<NODEDATA>::perturb_node_heights(vector<iterator<NODEDATA> > trav, double epsilon, gsl_rng *r){

        	// perturb the node heights, they're in even positions
        	for(int i = 0 ; i < trav.size(); i+=2){
        		double d1 = 0;
        		double d2 = 0;
        		double max = 0;
        		if (i-1 >=0){
        			d1 = trav[i-1]->m_time;
        		}
        		if (i+1 < trav.size()){
        			d2 = trav[i+1]->m_time;
        		}
        		max = d1;
        		if (d2 > max ) max = d2;
        		double toadd = (2* gsl_rng_uniform(r) - 1)*epsilon;
        		double newheight = trav[i]->m_time + toadd;
        		if (newheight < max) newheight = max+ (max-newheight);
        		trav[i]->m_time = newheight;
        	}

        	//perturb the heights of the interior nodes
        	for (int i = 1 ; i < trav.size(); i+=2){
           		double d1 = 10000;
           		double d2 = 10000;
           		double min = 0;
          		if (i-1 >=0){
          			d1 = trav[i-1]->m_time;
          		}
          		if (i+1 < trav.size()){
          			d2 = trav[i+1]->m_time;
          		}
          		min = d1;
          		if (d2 < min ) min = d2;
          		double toadd = (2* gsl_rng_uniform(r) - 1)*epsilon;
          		double newheight = trav[i]->m_time + toadd;
          		if (newheight > min) newheight = min - (newheight-min);
          		trav[i]->m_time = newheight;
        	}
        }
        template<class NODEDATA>
        void Tree<NODEDATA>::build_tree(vector<iterator<NODEDATA> > trav){
        	//
        	//first find the new root (the minimum height)
        	//
        	iterator<NODEDATA> newroot = trav[0];
        	double minheight = trav[0]->m_time;
        	int minpos = 0;
        	for(int i = 0 ; i < trav.size(); i ++){
        		if (trav[i]->m_time < minheight) {
        			newroot = trav[i];
        			minheight = trav[i]->m_time;
        			minpos =i;
        		}
        	}
        	newroot.setFather(0);
        	setRoot(newroot);

        	build_tree_helper(&trav, minpos);
        	// now
        }

        template<class NODEDATA>
        void Tree<NODEDATA>::build_tree_helper(vector<iterator<NODEDATA> >* trav, int index){
        	//
        	// look left
        	//
        	//cout << index << "\n"; cout.flush();
        	bool leftborder = false;
        	bool rightborder = false;
        	double leftmin = 10000;
        	double rightmin = 10000;
        	double minindex = 0;
        	int i = index-1;
        	bool foundleft = false;
        	while (i >= 0 && leftborder ==false){
        		if (trav->at(i)->m_time < trav->at(index)->m_time){
        			leftborder = true;
        			continue;
        		}
        		if (trav->at(i)->m_time < leftmin){
        			minindex = i;
        			leftmin = trav->at(i)->m_time;
        			foundleft = true;
        		}
        		i--;
        	}
        	//cout <<"here "<< minindex << " "<< foundleft << "\n"; cout.flush();
        	//
        	// set the pointers
        	//
        	if (foundleft){
        		//cout <<"here2 \n"; cout.flush();
        		trav->at(index).setFirstChild(trav->at(minindex));

        		trav->at(minindex).setFather(trav->at(index));
        		build_tree_helper(trav, minindex);
        	}
        	else{
        		trav->at(index).setFirstChild(0);
        	}

        	//
        	// look right
        	//
        	minindex = 0;
        	i = index+1;
        	bool foundright = 0;
          	while (i < trav->size() && rightborder ==false){
            		if (trav->at(i)->m_time < trav->at(index)->m_time){
            			rightborder = true;
            			continue;
            		}
            		if (trav->at(i)->m_time < rightmin){
            			minindex = i;
            			rightmin = trav->at(i)->m_time;
            			foundright = true;
            		}
            		i++;
          	}
           	//
          	// set the pointers
          	//
          	if (foundright){
          		trav->at(index).setLastChild(trav->at(minindex));
          		trav->at(minindex).setFather(trav->at(index));
          		build_tree_helper(trav, minindex);
			}
          	else{
          		trav->at(index).setLastChild(0);
          	}

        }
        template<class NODEDATA>
        void Tree<NODEDATA>::update_branch_lengths(iterator<NODEDATA> rootIterator){
        	if (rootIterator != NULL){
        		update_branch_lengths(rootIterator.getFirstChild());
        		if (rootIterator.hasFather() == false) rootIterator->m_len = -1;
        		else	rootIterator->m_len = rootIterator->m_time - rootIterator.getFather()->m_time;
        		update_branch_lengths(rootIterator.getLastChild());
        	}
        }
}

