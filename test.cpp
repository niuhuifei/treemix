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
	string testnewick = "(YRI:0.8,(CEU:0.4,(JPT:0.2,CHB:0.3):0.2):0.5);";
	/*map<string, int> testmap;
	testmap["YRI"] = 0;
	testmap["CEU"] = 1;
	testmap["JPT"] = 2;
	testmap["CHB"] = 3;
	PhyloPop_Tree::Node<PhyloPop_Tree::NodeData> test();
	PhyloPop_Tree::BinaryTree<PhyloPop_Tree::NodeData> testtree(testnewick, testmap);
	//testtree.print_inorder(testtree.getRoot());
	map<int, PhyloPop_Tree::iterator<PhyloPop_Tree::NodeData> > tmp = testtree.get_tips(testtree.getRoot());
	for (map<int, PhyloPop_Tree::iterator<PhyloPop_Tree::NodeData> >::iterator it = tmp.begin(); it != tmp.end(); it++){
		cout << it->first << " "<< it->second->m_len<< "\n";
	}

	double tmpdouble = testtree.get_dist_to_root(tmp[3]);

	//PhyloPop_Tree::iterator<PhyloPop_Tree::NodeData> it = tmp[1];


	PhyloPop_Tree::iterator<PhyloPop_Tree::NodeData> tmpnode = testtree.get_LCA(testtree.getRoot(), tmp[3], tmp[2]);

	cout << tmpdouble << " "<< testtree.get_dist_to_root(tmpnode) <<  "\n";
	//cout << tmpnode->m_id << " "<< tmpnode->m_len<< "\n";

	//PhyloPop_Tree::iterator<PhyloPop_Tree::NodeData> testit = testtree.getRoot();
	//while (testit.hasChildNodes()){ testit = testit.getFirstChild();}
	//cout << testit->m_id << "\n" << testit->m_len << "\n";
	return 0;

	*/
	CountData counts("testin.gz");
	State teststate(testnewick, &counts);
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
};
