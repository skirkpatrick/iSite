/*
    iSite data structure. Models the network.
*/

#ifndef ISITE_GRAPH_HPP
#define ISITE_GRAPH_HPP

#include <boost/graph/adjacency_list.hpp>
#include <vector>
#include <iostream>

#include "graph_properties.hpp"
#include "Random.hpp"


class isite_graph
    : public boost::adjacency_list<boost::vecS,
                                   boost::vecS,
                                   boost::undirectedS,
                                   vertexsites>
{
    private:
        double prob_loss;
        double prob_asym;
        double prob_self;
        double prob_fusion;
        double prob_fission;

        void _progenation();
        void _fusion();
        void _fission();
        void _print_edge();
        void _add_age();

    public:
        isite_graph();
        ~isite_graph();
        void generation();
        void domain_distribution(std::vector<int>&);
        isite_graph get_large_component();
        void simplify();
        int triangles();
        int triples();
        int components();
        int order();
        void print(std::ostream&);
        
};


#endif
