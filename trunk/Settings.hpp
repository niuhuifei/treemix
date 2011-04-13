/*
 * Settings.hpp
 *
 *  Created on: Mar 29, 2011
 *      Author: pickrell
 */

#ifndef SETTINGS_HPP_
#define SETTINGS_HPP_

#include <string>
#include <list>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <sstream>
#include <fstream>
#include <ostream>
#include <cstdlib>
#include <cmath>
#include <stack>
#include <valarray>
#include <utility>
#include "gzstream.h"
#include <sys/stat.h>
#include <boost/tokenizer.hpp>
#include "mvnorm.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#include <math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_randist.h>
#include <omp.h>
#include "CmdLine.h"


using std::string;
using std::vector;
using std::list;
using std::stack;
using std::map;
using std::set;
using std::multiset;
using std::valarray;
using std::cout;
using std::cin;
using std::endl;
using std::ostream;
using std::ofstream;
using std::stringstream;
using std::pair;
using std::iterator;
using std::pair;
using std::make_pair;
using std::fstream;
using std::ifstream;

typedef vector<vector<double> >                         vector2d;
typedef vector<vector<vector<double> > >        vector3d;
#endif /* SETTINGS_HPP_ */
