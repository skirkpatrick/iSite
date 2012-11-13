#ifndef GRAPH_PROPERTIES_HPP
#define GRAPH_PROPERTIES_HPP

#include <boost/graph/adjacency_list.hpp>
#include <boost/property_map/property_map.hpp>
#include <vector>
#include <map>
#include <utility>
#include <string>

typedef boost::adjacency_list_traits<boost::listS,boost::listS,boost::undirectedS>::edge_descriptor edge_descriptor;
typedef boost::adjacency_list_traits<boost::listS,boost::listS,boost::undirectedS>::vertex_descriptor vertex_descriptor;


struct isite
{
    std::vector<edge_descriptor> edges;
    unsigned int age;
    //string site_name;

    isite() : /*site_name(""),*/ edges(), age(0)
    {
    }
    isite(const isite& is) : /*site_name(is.site_name),*/ edges(is.edges), age(is.age)
    {
    }
    isite(edge_descriptor e, unsigned int a, std::string name) : /*site_name(name),*/ age(a)
    {
        edges.push_back(e);
    }
};

struct vertexsites
{
    std::map<edge_descriptor, std::pair<int, int> > selfLoops;
    std::vector<isite> sites;
    std::map<edge_descriptor,int> edgeToSite;
    std::map<std::string, int> site_name_to_index;
};


typedef boost::adjacency_list<boost::listS,boost::listS,boost::undirectedS,boost::property<boost::vertex_index_t,int,vertexsites> > Graph;

typedef boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;

typedef boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
typedef boost::graph_traits<Graph>::edge_iterator edge_iterator;
typedef boost::graph_traits<Graph>::adjacency_iterator AdjIter;
typedef boost::property_map<Graph, boost::vertex_index_t>::type vimap;

#endif
