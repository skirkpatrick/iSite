//iaddedaline!
#include<cstdio>    //fscanf can be useful if sites are read as strings
#include<cstdlib>   //for atoi
#include<iostream>  //for input and output
#include<fstream>   //for file input and output
#include<utility>   //for std::pair
#include<vector>    //for std::vector
#include<map>       //for std::map
#include<cstring>   //for std::string

#include<boost/graph/adjacency_list.hpp>
#include<boost/graph/connected_components.hpp>

using namespace std;
using namespace boost;


typedef adjacency_list_traits<listS,vecS,undirectedS>::edge_descriptor edge_descriptor;
typedef adjacency_list_traits<listS,vecS,undirectedS>::vertex_descriptor vertex_descriptor;

struct isite
{
    vector<edge_descriptor> edges;
    unsigned int age;
};

struct vertexsites
{
    vector<isite> sites;
};



typedef adjacency_list<listS, vecS, undirectedS, vertexsites> Graph;

/*
adjacency_list<OutEdgeList,         <-listS=time, vecS=space
                VertexList,         <-listS=time, vecS=space
                Directed,           <-iSite graphs are undirected
                VertexProperties,   <-Possibly not needed?
                EdgeProperties,     <-Needed for iSites
                GraphProperties,    <-Possibly not needed?
                EdgeList>           <-Not needed
                */

typedef graph_traits<Graph>::vertex_iterator vertex_iterator;
typedef graph_traits<Graph>::edge_iterator edge_iterator;
typedef graph_traits<Graph>::adjacency_iterator AdjIter;


/*
    Counts the number of triangles in a graph
*/
int triangles(Graph& graph)
{
    vertex_iterator vi, viend;
    AdjIter adji, adjiend, adji2;
    int numTriangles=0;

    //Cycle through all vertices
    for (tie(vi, viend)=vertices(graph); vi!=viend; vi++)
    {
        //Cycle through adjacent vertices
        for (tie(adji, adjiend)=adjacent_vertices(*vi, graph); adji!=adjiend; adji++)
        {
            //Check for edges between vertices adjacent to vi aka triangle
            for (adji2=adji/*works with self-loops?*/; adji2!=adjiend; ++adji2)
            {
                if (edge(*adji, *adji2, graph).second)
                numTriangles++;
            }
        }
    }
    return numTriangles/3;
}


/*
    Counts the number of triples in a graph
*/
int triples(Graph& graph)
{
    vertex_iterator vi, viend;
    int numTriples=0;

    for (tie(vi, viend)=vertices(graph); vi!=viend; vi++)
        numTriples+=((degree(*vi, graph)*(degree(*vi,graph)-1))/2);
    return numTriples;
}


/*
    Counts the number of components in a graph
*/
int components(Graph& graph)
{
    vector<int> c(num_vertices(graph));
    return connected_components(graph, &c[0]);//this is an easier way I found
}



//Stores input parameters
struct parameters
{
    string infile;
    double prob_loss;
    double prob_asym;
    int end_order;
    int iterations;
} param;

pair<int, int> duplication(Graph& graph)
{
    srand(time(NULL));
    int parent = (int)((float)rand()/(RAND_MAX)*(boost::num_vertices(graph)-1) );
//the graph -1 and node +1 make it so the range is between 1 to last_node 
    parent+=1;
    cout << "parent: " << parent << endl;

    //get a new vertex
    Graph::vertex_descriptor v_description = *child;
    child = boost::add_vertex(graph);
    cout << "added vertex: " << *child << endl;

    //give all of the selectedNode edges to the cloned node
    out_edge_iterator oei, oeiend;
    for( tie(oei, oeiend) = out_edges(parent, graph); oei != oeiend; ++oei )
    {
        boost::add_edge(target(*oei,graph), *child, graph);
        cout << "adding: " << target(*oei,graph) << " to: " << *child << endl;
    }
//return parent, child
return pair<parent, *child>
}

//Perform iSite algorithm
void isite(Graph& graph)
{
    pair<int,int> vertices;
    vertices=duplication(graph);

    
}


int main(int argc, char* argv[])
{
    if (argc!=6)
    {
        cerr<<"Usage: ./iSite <seed-graph> <probability of loss> "
                "<probability of assymetry> <end order> <iterations>"<<endl;
        exit(1);
    }

    Graph graph;                        //Protein network
    param.infile=argv[1];               //Seed graph
    param.prob_loss=atof(argv[2]);      //Probability of loss of redundancy
    param.prob_asym=atof(argv[3]);      //Probability of assymetry
    param.end_order=atoi(argv[4]);      //Order at which to stop
    param.iterations=atoi(argv[5]);     //Number of times to run algorithm

    ifstream infile(param.infile.c_str());
    if (!infile)
    {
        cerr<<"Error opening specified seed graph!"<<endl;
        exit(1);
    }

    map<string,int> nodes;
    int counter=0;
    string v1,v2;

    //Building seed graph
    while(!(infile>>v1>>v2).eof())
	{
		if(nodes[v1]==0)
		{
			nodes[v1]=++count;
			add_vertex(graph);
		}
		if(nodes[v2]==0)
		{
			nodes[v2]=++count;
			add_vertex(graph);
		}

		add_edge(nodes[v1],nodes[v2],graph);
		//loop to add 1 to all member vectices' age.
		
	}

    while (num_vertices(graph)!=param.end_order)
    {
        
    }



    //Output
    ofstream outfile("result");
    if (!outfile)
    {
        cerr<<"Error opening output file: result"<<endl;
        exit(1);
    }
	
    outfile<<param.prob_loss<<" ";
    outfile<<param.prob_asym<<" ";
    outfile<<num_vertices(graph)<<" ";
    outfile<<num_edges(graph)<<" ";
    int triangles=triangles(graph);
    int triples=triples(graph);
    outfile<<triangles<<" ";
    outfile<<triples<<" ";
    outfile<<(3*triangles)/triples<<" ";
    outifle<<components(graph);



    outfile.close();
    infile.close();
}
