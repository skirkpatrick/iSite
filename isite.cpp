#include<cstdio>    //fscanf can be useful if sites are read as strings
#include<cstdlib>   //for atoi
#include<iostream>  //for input and output
#include<fstream>   //for file input and output
#include<utility>   //for std::pair
#include<vector>    //for std::vector
#include<map>       //for std::map
#include<cstring>   //for std::string
#include<algorithm> //for copying vectors


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
    isite() : age(0), edges()
    {
    }
    isite(edge_descriptor e, unsigned int a) : age(a)
    {
        edges.push_back(e);
    }
};

struct vertexsites
{
    vector<isite> sites;
};

typedef adjacency_list<listS, vecS, undirectedS, vertexsites> Graph;

typedef boost::graph_traits<Graph>::out_edge_iterator out_edge_iterator;

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
int countTriples(Graph& graph)
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

pair<Graph::vertex_descriptor, Graph::vertex_descriptor> duplicate(Graph& graph)
{
    srand(time(NULL));
    int parent = (int)((float)rand()/(RAND_MAX)*(boost::num_vertices(graph)-1) );
//the (graph -1) and (node +1) make it so the range is between 1 to last_node 
    parent+=1;
    cout << "parent: " << parent << endl;

//get a new vertex
    vertex_iterator child;
    child = boost::add_vertex(graph);
    Graph::vertex_descriptor child_description = *child;
    cout << "added vertex: " << *child << endl;

//give all of the parent edges to the child
    out_edge_iterator oei, oeiend;
    for( tie(oei, oeiend) = out_edges(parent, graph); oei != oeiend; ++oei )
    {
        boost::add_edge(target(*oei,graph), *child, graph);
        cout << "adding: " << target(*oei,graph) << " to: " << *child << endl;
    }

//get vertex_descriptor to parent node
    tie(oei, oeiend) = out_edges(parent, graph);
    Graph::vertex_descriptor parent_description = source(*oei, graph);

//clone each of the parent's isites and give it to the child
    int numSites = graph[parent_description].sites.size();
    for( int i=0; i < numSites; i++)
    {
    /*
        isite* newSite = new isite(); 
        newSite->age = 0;
        copy( graph[parent_description].sites[i].edges.begin(), graph[parent_description].sites[i].edges.end(), newSite->edges.begin() ); 
        graph[child_description].sites.push_back(*newSite);
        */
        isite newSite;
        newSite.age=0;
        newSite.edges.insert(newSite.edges.begin(),graph[parent_description].sites[i].edges.begin(),graph[parent_description].sites[i].edges.end());
        //If we'd rather use copy() instead of insert, we need to first do a
            //resize() or we get a segfault

        graph[child_description].sites.push_back(newSite);
//SHOULD DO MEMORY CLEAN UP HERE???
    //Got rid of the new operator. No worries.
    }

//return parent, child
    pair<Graph::vertex_descriptor, Graph::vertex_descriptor> tmpPair;
    tmpPair = std::make_pair(parent_description, child_description);
return tmpPair; 
}

//Perform iSite algorithm
void duplication(Graph& graph)
{
    pair<Graph::vertex_descriptor,Graph::vertex_descriptor> vertices;
    vertices=duplicate(graph);

    
}

void addAge(Graph& graph)
{
    vertex_iterator vi, viend;
    for(tie(vi, viend) = vertices(graph); vi != viend; ++vi) 
    {
//loop through each site on the vertex and add 1 to it
        Graph::vertex_descriptor v_description = *vi;
        int numSites = graph[v_description].sites.size();
        for( int i=0; i < numSites; i++)
        {
            graph[v_description].sites[i].age += 1;
        }
    }

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
        //Create nodes if don't exist
		if(nodes[v1]==0)
		{
			nodes[v1]=++counter;
            add_vertex(graph);
		}
		if(nodes[v2]==0)
		{
			nodes[v2]=++counter;
			add_vertex(graph);
		}

        //Add edge
        Graph::vertex_descriptor vd1 = vertex(nodes[v1],graph);
        Graph::vertex_descriptor vd2 = vertex(nodes[v2],graph);
        Graph::edge_descriptor ed;
        bool temp_bool;

        tie(ed, temp_bool) = add_edge(vd1,vd2,graph);

        //Add iSite
        graph[vd1].sites.push_back(isite(ed,0));
        graph[vd2].sites.push_back(isite(ed,0));
		
	}

#ifdef DEBUG
    vertex_iterator vi, viend;
    for (tie(vi,viend) = vertices(graph); vi!=viend; ++vi)
    {
        cout<<*vi;
        out_edge_iterator oei, oeiend;
        for (tie(oei,oeiend)=out_edges(*vi,graph); oei!=oeiend; ++oei)
            cout<<"->"<<target(*oei,graph);
        cout<<endl;
    }
#endif

    while (num_vertices(graph)!=param.end_order+1)
    {
        duplication(graph);
    }

#ifdef DEBUG
    for (tie(vi,viend) = vertices(graph); vi!=viend; ++vi)
    {
        cout<<*vi;
        out_edge_iterator oei, oeiend;
        for (tie(oei,oeiend)=out_edges(*vi,graph); oei!=oeiend; ++oei)
            cout<<"->"<<target(*oei,graph);
        cout<<endl;
    }
#endif



    //Output
    ofstream outfile("result");
    if (!outfile)
    {
        cerr<<"Error opening output file: result"<<endl;
        exit(1);
    }
	
    outfile<<param.prob_loss<<" ";
    outfile<<param.prob_asym<<" ";
    outfile<<num_vertices(graph)-1<<" ";
    outfile<<num_edges(graph)<<" ";
    int numTriangles=triangles(graph);
    int numTriples=countTriples(graph);
    outfile<<numTriangles<<" ";
    outfile<<numTriples<<" ";
    outfile<<(double)(3*numTriangles)/numTriples<<" ";
    outfile<<components(graph);



    outfile.close();
    infile.close();
}
