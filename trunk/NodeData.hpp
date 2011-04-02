/*
 * NodeData.hpp
 *
 *  Created on: Mar 29, 2011
 *      Author: pickrell
 */

#ifndef NODEDATA_HPP_
#define NODEDATA_HPP_


#include "Settings.hpp"

namespace PhyloPop_Tree
{
        ///NodeData is the Storage for the Tree Container.
        /**
         * NodeData is the Storage for the Tree Container.
         * All Attributes are public.
         */
        struct NodeData
        {
                virtual ~NodeData(){};

                int m_id; //id
                double m_len; //length to father node
                double m_time; //node timepoint
                //vector<float> m_sigma; // sigmas at node (if tip)
                //vector<pair<int, int> > m_counts; // allele counts at node (if tip)

        };//NodeData

}//PhyloPop_Tree

#endif /* NODEDATA_HPP_ */
