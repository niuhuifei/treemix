/*
 * test.cpp
 *
 *  Created on: Mar 29, 2011
 *      Author: pickrell
 */

//#include "Tree.h"
//#include "BinaryTree.h"
#include "State.h"
#include "CountData.h"

int main(){
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_ranlxs2;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, 100);
	string testnewick = "(YRI:0.8,(CEU:0.4,(JPT:0.2,CHB:0.3):0.2):0.5);";

	CountData counts("testin.gz");
	ofstream tmp("tmp");
	MCMC_params p;
	State teststate(testnewick, &counts, &p);
	vector<PhyloPop_Tree::iterator< PhyloPop_Tree::NodeData> > vec = teststate.tree->get_inorder_traversal(counts.npop);
	for (int i = 0; i < vec.size(); i++){
		cout << vec[i]->m_id << " "<< vec[i]->m_len << " "<< vec[i]->m_time<< "\n";
	}
	teststate.tree->set_node_heights(vec);
	cout << "\n";
	for (int i = 0; i < vec.size(); i++){
			cout << vec[i]->m_id << " "<< vec[i]->m_len << " "<< vec[i]->m_time<< "\n";
		}
	cout <<"\n";

	//for(int i = 0; i < 10; i++){
	//	teststate.propose_tree(r);
	//}
	//vec = teststate.tree->get_inorder_traversal(counts.npop);
	//for (int i = 0; i < vec.size(); i++){
	//	cout << vec[i]->m_id << " "<< vec[i]->m_len << " "<< vec[i]->m_time<< "\n";
	//}
	cout << "\n";
	teststate.compute_sigma();
	for (int i = 0 ; i < 10000; i++){
	teststate.update(r);
	tmp << gsl_vector_get(teststate.means, 0) << "\n";
	}
	cout << "\n";
	vec = teststate.tree->get_inorder_traversal(counts.npop);
	for (int i = 0; i < vec.size(); i++){
			cout << vec[i]->m_id << " "<< vec[i]->m_len << " "<< vec[i]->m_time<< "\n";
		}

	cout << "llik: "<< teststate.llik() << "\n";
	/*
	map<int, PhyloPop_Tree::iterator<PhyloPop_Tree::NodeData> > tmp = teststate.tree->get_tips(teststate.tree->getRoot());
	for (map<int, PhyloPop_Tree::iterator<PhyloPop_Tree::NodeData> >::iterator it = tmp.begin(); it != tmp.end(); it++){
		cout << it->first << " "<< it->second->m_len<< "\n";
	}
	cout << "\n";
	State teststate2(teststate);
	teststate2.tree->print_inorder(teststate2.tree->getRoot());
	cout <<"\n";
	teststate.compute_sigma();
	cout << "here\n"; cout.flush();
	teststate.print_sigma();
	double tmplik = teststate.llik();
	cout << "llik: "<< tmplik << "\n";
	*/
};
