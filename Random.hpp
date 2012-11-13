#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/graph/random.hpp>

#include "graph_properties.hpp"


class Random
{
    private:
        boost::mt19937 rng;
        boost::uniform_real<> real_dist;
        boost::variate_generator<boost::mt19937&, boost::uniform_real<> > rand_real;
        unsigned int counter;

    public:
        Random(unsigned int sd)
            : rng(sd), real_dist(0.0, 1.0), rand_real(rng, real_dist), counter(0) {}
            
        Graph::vertex_descriptor random_vertex(Graph& graph)
        {
            ++counter;
            return boost::random_vertex(graph, rng);
        }
        
        double rand()
        {
            ++counter;
            return rand_real();
        }
};

#endif
