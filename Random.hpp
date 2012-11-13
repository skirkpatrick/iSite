#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <boost/random.hpp>
#include <boost/generator_iterator.hpp>
#include <boost/graph/random.hpp>
#include <fstream>

#include "graph_properties.hpp"


class Random
{
    private:
        boost::mt19937 rng;
        boost::uniform_real<> real_dist;
        boost::variate_generator<boost::mt19937&, boost::uniform_real<> > rand_real;
        unsigned int counter;
        std::ifstream inf;

    public:
        Random(unsigned int sd, const char* file = "")
            : rng(sd), real_dist(0.0, 1.0), rand_real(rng, real_dist), counter(0), inf(file) {}
            
        Graph::vertex_descriptor random_vertex(Graph& graph)
        {
            ++counter;
            return boost::random_vertex(graph, rng);
        }
        
        double rand()
        {
            ++counter;
            if (inf)
            {
                double res;
                inf >> res;
                return res;
            }
            return rand_real();
        }

        unsigned int get_count() const
        {
            return counter;
        }
};

#endif
