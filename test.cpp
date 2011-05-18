/*
 * test.cpp
 *
 *  Created on: Mar 29, 2011
 *      Author: pickrell
 */

//#include "Tree.h"
//#include "BinaryTree.h"
#include "MCMC.h"
#include "CountData.h"
#include "WishartState.h"
#include "GraphState.h"

int main(){
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();
	T = gsl_rng_ranlxs2;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, 100);


	//string testnewick = "((YRI:0.8,LWD:0.4):0.1,(CEU:0.4,(JPT:0.2,CHB:0.3):0.2):0.5);";
	//string testnewick = "((((pop12:0.00906342,pop11:0.00939995):0.000409853,(pop9:0.0106216,pop10:0.010027):0.00037873):0.000125612,((pop7:0.0121508,pop8:0.0111988):0.000449598,((pop5:0.0126136,(((pop1:0.0156482,pop2:0.0141466):0.000919802,pop3:0.0139073):0.000560655,pop4:0.0133937):0.000572252):0.000359304,pop6:0.012505):0.000776528):0.000163241):4.31837e-07,((pop14:0.00787589,((pop16:0.00661613,(((pop19:0.00494392,pop20:0.00515132):0.000711535,pop18:0.00547178):0.000608346,pop17:0.00587652):0.000647406):0.000533695,pop15:0.0073544):0.000446108):0.000422077,pop13:0.0088081):0.000607757)";
	//string testnewick = "(pop2:0.01405,(pop3:0.01374,(pop4:0.01314,(pop5:0.01239,(pop6:0.01213,(pop7:0.01165,(pop8:0.01089,(pop9:0.01020,(pop10:0.00970,(pop11:0.00916,(pop12:0.00868,(pop13:0.00830,(pop14:0.00761,(pop15:0.00706,(pop16:0.00644,(pop17:0.00587,(pop18:0.00545,(pop19:0.00494,pop20:0.00514):0.00074):0.00065):0.00084):0.00069):0.00069):0.00077):0.00075):0.00073):0.00074):0.00081):0.00076):0.00079):0.00086):0.00077):0.00080):0.00079):0.00109,pop1:0.01576)";
	string testnewick = "(pop1:0.01576,(pop2:0.01405,(pop3:0.01374,(pop4:0.01314,(pop5:0.01239,(pop6:0.01213,(pop7:0.01165,(pop8:0.01089,(pop9:0.01020,(pop10:0.00970,(pop11:0.00916,(pop12:0.00868,(pop13:0.00830,(pop14:0.00761,(pop15:0.00706,(pop16:0.00644,(pop17:0.00587,(pop18:0.00545,(pop19:0.00494,pop20:0.00514):0.00074):0.00065):0.00084):0.00069):0.00069):0.00077):0.00075):0.00073):0.00074):0.00081):0.00076):0.00079):0.00086):0.00077):0.00080):0.00079):0.00109):0.000001)";
	ogzstream treeout("treeout.gz");
	//ogzstream mout("meanout.gz");
	CountData counts("/Users/pickrell/Desktop/rosenburg_model_20pop2_phylopop_in.gz", 1);
	//CountData counts("testin_counts_scale.gz", 1);
	MCMC_params p;

	GraphState teststate(testnewick, &counts, &p);
	cout << teststate.llik() <<"\n";
	teststate.set_branches_ls();
	teststate.compute_sigma();
	cout << teststate.llik() << "\n";
	cout << teststate.tree->get_newick_format() << "\n";
	//teststate.tree->print();
	//vector< Graph::vertex_descriptor > vec = teststate.tree->get_inorder_traversal(5);
	//teststate.tree->local_update(vec[1], r);
	//teststate.tree->print();
	//teststate.compute_sigma();
	vector< Graph::vertex_descriptor > vec = teststate.tree->get_inorder_traversal(20);
	cout << teststate.llik() << "\n";
	//cout << teststate.tree->g[teststate.tree->root].index;
	//teststate.tree->set_root(vec[1]);
	///teststate.compute_sigma();
	//cout << teststate.llik() << "\n";
	//teststate.tree->set_root(vec[5]);
	//teststate.compute_sigma();
	testnewick = "((((pop12:0.00906342,pop11:0.00939995):0.000409853,(pop9:0.0106216,pop10:0.010027):0.00037873):0.000125612,((pop7:0.0121508,pop8:0.0111988):0.000449598,((pop5:0.0126136,(((pop1:0.0156482,pop2:0.0141466):0.000919802,pop3:0.0139073):0.000560655,pop4:0.0133937):0.000572252):0.000359304,pop6:0.012505):0.000776528):0.000163241):4.31837e-07,((pop14:0.00787589,((pop16:0.00661613,(((pop19:0.00494392,pop20:0.00515132):0.000711535,pop18:0.00547178):0.000608346,pop17:0.00587652):0.000647406):0.000533695,pop15:0.0073544):0.000446108):0.000422077,pop13:0.0088081):0.000607757)";
	GraphState teststate2(testnewick, &counts, &p);
	cout << teststate2.llik() << " 2\n";
	teststate2.set_branches_ls();
	teststate2.compute_sigma();
	cout << teststate2.llik() << " 2\n";

	//cout << teststate.tree->get_newick_format() << "\n";
//
	//for (int i = 0; i < vec.size(); i++){
	//	cout << i << " "<< teststate.tree->g[vec[i]].name << " " << teststate.tree->g[vec[i]].index<< "\n";
	//}
	//teststate.tree->local_update(vec[5], r);
	//teststate.tree->print();
	//teststate.compute_sigma();
	//cout << teststate.llik() <<"\n";
	//cout << teststate.tree->get_newick_format() << "\n";
/*
 * teststate.tree->local_update(vec[1], r);
	teststate.compute_sigma();
	cout << teststate.llik() <<"\n";
	cout << teststate.tree->get_newick_format() << "\n";
	vec = teststate.tree->get_inorder_traversal(20);
	for (int i = 0; i < vec.size(); i++){
		cout << i << " "<< teststate.tree->g[vec[i]].name << " "<< teststate.tree->g[vec[i]].index << "\n";
	}
	cout << teststate.tree->get_newick_format() << "\n";
	teststate.tree->print();
	Graph::edge_descriptor ed = edge( vec[9], vec[19], teststate.tree->g).first;
	teststate.tree->g[ed].len = teststate.tree->g[ed].len - 0.00043;
	ed = edge( vec[21], vec[22], teststate.tree->g).first;
	teststate.tree->g[ed].len = teststate.tree->g[ed].len + 0.00043;
	teststate.compute_sigma();
	cout << teststate.llik() <<"\n";
	teststate.tree->print();
*/
};
