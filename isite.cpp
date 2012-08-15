#include<cstdio>    //fscanf can be useful if sites are read as strings
#include<cstdlib>   //for atoi
#include <time.h>   //for random number generation
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
    map<edge_descriptor,int> edgeToSite;
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
        for(tie(adji,adjiend)=adjacent_vertices(*vi,graph);adji!=adjiend;adji++)
        {
            //Check for edges between vertices adjacent to vi aka triangle
            for (adji2=adji; adji2!=adjiend; ++adji2)
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
    return connected_components(graph, &c[0]);
}



//Stores input parameters
struct parameters
{
    string infile;
    double prob_loss;
    double prob_asym;
    unsigned int end_order;
    unsigned int iterations;
} param;

//Given an edge and a vertex, returns the vertex on the opposite end
Graph::vertex_descriptor edgeDest(const Graph::vertex_descriptor vd,
                                  const Graph::edge_descriptor ed,
                                  const Graph& graph)
{
    return (vd==source(ed, graph)) ? target(ed, graph) : source(ed, graph);
}

pair<Graph::vertex_descriptor,Graph::vertex_descriptor> duplicate(Graph& graph)
{
    int parent = rand() % (num_vertices(graph)-1);
//the (graph -1) and (node +1) make it so the range is between 1 to last_node 
    parent+=1;
    cout << "node chosen for duplication (parent): " << parent << endl;
    Graph::vertex_descriptor parent_description = vertex(parent, graph);

    //get a new vertex
    vertex_iterator child;
    child = boost::add_vertex(graph);
    Graph::vertex_descriptor child_description = *child;
    cout << "duplicated node (child): " << *child << endl;

    //copy edges and isites
    int numSites = graph[parent_description].sites.size();
    for (int i=0; i < numSites; i++)
    {
        isite newSite;
        newSite.age=0;

        int siteEdges = graph[parent_description].sites[i].edges.size();
        for (int j=0; j < siteEdges; j++)
        {
            
            Graph::edge_descriptor ed, ed2;
            bool temp_bool;
            Graph::vertex_descriptor vd = edgeDest(parent_description,
                graph[parent_description].sites[i].edges[j], graph);

            tie(ed, temp_bool) = add_edge(child_description,vd,graph);

            //update child iSite
            newSite.edges.push_back(ed);
            graph[child_description].edgeToSite[ed]=i;

            //update other node's iSite
            tie(ed2, temp_bool) = edge(parent_description, vd, graph);
            graph[vd].sites[graph[vd].edgeToSite[ed2]].edges.push_back(ed);
            graph[vd].edgeToSite[ed]=graph[vd].edgeToSite[ed2];

            //self-loop (not yet implemented)
                //may need to be at beginning of loop
            if (target(ed, graph)==source(ed, graph));
        }

        graph[child_description].sites.push_back(newSite);
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

    //prob_asym
    int numSites = graph[vertices.first].sites.size();

#ifdef DEBUG
    cout<<"numSites: "<<numSites<<endl;
#endif

    for (int i=0; i<numSites; i++)
    {
        Graph::vertex_descriptor vertexLoss;

        double rand_res = ((double) rand()) / (double)RAND_MAX;
        if (rand_res <= param.prob_asym) //Parent loss
        {
#ifdef DEBUG
            cout<<"prob_asym result: Parent"<<endl;
#endif
            vertexLoss = vertices.first;
        }
        else //Child loss
        {
#ifdef DEBUG
            cout<<"prob_asym result: Child"<<endl;
#endif
            vertexLoss = vertices.second;
        }

        //prob_loss
        int numEdges = graph[vertexLoss].sites[i].edges.size();
        for (int j=0; j<numEdges; j++)
        {
#ifdef DEBUG
            cout<<"numEdges: "<<numEdges<<endl;
#endif
            rand_res = ((double) rand()) / (double)RAND_MAX;
            if (rand_res <= param.prob_loss) //Edge is lost
            {
#ifdef DEBUG
                cout<<"prob_loss: Yes"<<endl;
#endif
                Graph::edge_descriptor edgeLoss;
                edgeLoss = graph[vertexLoss].sites[i].edges[j];

                //Remove from vertexLoss
                graph[vertexLoss].edgeToSite.erase(edgeLoss);
                graph[vertexLoss].sites[i].edges.erase(
                    graph[vertexLoss].sites[i].edges.begin()+j);

                //Update connected vertex
                Graph::vertex_descriptor connectedVertex;
                connectedVertex = edgeDest(vertexLoss, edgeLoss, graph);
                int connectedSite = graph[connectedVertex].edgeToSite[edgeLoss];
                int connectedSiteSize = graph[connectedVertex].sites[connectedSite].edges.size();
                int k;
                for (k=0; k<connectedSiteSize; k++)
                {
                    if (graph[connectedVertex].sites[connectedSite].edges[k]==edgeLoss)
                        break;
                }
                graph[connectedVertex].edgeToSite.erase(edgeLoss);
                graph[connectedVertex].sites[connectedSite].edges.erase(
                    graph[connectedVertex].sites[connectedSite].edges.begin()+k);

                //If connected iSite is empty, remove it
                if (graph[connectedVertex].sites[connectedSite].edges.size()==0
                    && connectedVertex != vertices.first
                    && connectedVertex != vertices.second)
                {
                    graph[connectedVertex].sites.erase(
                        graph[connectedVertex].sites.begin() + connectedSite);
                }

                remove_edge(edgeLoss, graph);
                j--;
                numEdges--;
            }
#ifdef DEBUG
            else cout<<"prob_loss: No"<<endl;
#endif
        }//cycle through edges

    }//cycle through iSites

    //Check parent and child for empty iSites
    int j=numSites;
    for (int i=0; i<j; i++)
        if (graph[vertices.first].sites[i].edges.size()==0)
        {
            graph[vertices.first].sites.erase(
                graph[vertices.first].sites.begin() + i);
            i--;
            j--;
        }
    for (int i=0; i<numSites; i++)
        if (graph[vertices.second].sites[i].edges.size()==0)
        {
            graph[vertices.second].sites.erase(
                graph[vertices.second].sites.begin() + i);
            i--;
            numSites--;
        }

    //If node has no iSites (and therefore no edges), delete it
    /******************implement here***********************/
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

    //Input validation
    if (param.prob_loss<0.0 || param.prob_loss>1.0)
    {
        cerr<<"Error: Probability of loss must be between 0 and 1"<<endl;
        exit(1);
    }
    if (param.prob_asym<0.0 || param.prob_loss>1.0)
    {
        cerr<<"Error: Probability of assymetry must be between 0 and 1"<<endl;
        exit(1);
    }
    if (param.iterations<0)
    {
        cerr<<"Error: Number of iterations must be a positive integer"<<endl;
        exit(1);
    }

    ifstream infile(param.infile.c_str());
    if (!infile)
    {
        cerr<<"Error opening specified seed graph!"<<endl;
        exit(1);
    }

    map<string,int> nodes;
    int counter=0;
    string v1,v2;
    srand(time(NULL));

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
        graph[vd1].edgeToSite[ed]=graph[vd1].sites.size()-1;
        graph[vd2].sites.push_back(isite(ed,0));
        graph[vd2].edgeToSite[ed]=graph[vd2].sites.size()-1;
		
	}

#ifdef DEBUG
    cout << "original graph" << endl;
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

    //Opening output file
    ofstream outfile("result");
    if (!outfile)
    {
        cerr<<"Error opening output file: result"<<endl;
        exit(1);
    }

    while (param.iterations--)
    {
        while (num_vertices(graph)<param.end_order+1)
        {
            addAge(graph);
            duplication(graph);
#ifdef DEBUG
            cout << "graph during algorithm" << endl;
            for (tie(vi,viend) = vertices(graph); vi!=viend; ++vi)
            {
                cout<<*vi;
                out_edge_iterator oei, oeiend;
                for (tie(oei,oeiend)=out_edges(*vi,graph); oei!=oeiend; ++oei)
                    cout<<"->"<<target(*oei,graph);
                cout<<endl;
            }
#endif

        }

#ifdef DEBUG
        cout << "graph after duplication of nodes to end_order" << endl;
        for (tie(vi,viend) = vertices(graph); vi!=viend; ++vi)
        {
            cout<<*vi;
            out_edge_iterator oei, oeiend;
            for (tie(oei,oeiend)=out_edges(*vi,graph); oei!=oeiend; ++oei)
                cout<<"->"<<target(*oei,graph);
            cout<<endl;
        }
#endif

#ifdef DEBUG
        cout << "listing of all nodes. "
            "For each node - the isite and the age of that site" << endl;
        for (tie(vi,viend) = vertices(graph); vi!=viend; ++vi)
        {
            Graph::vertex_descriptor vd1 = *vi; 
            cout<<"node: " << *vi << endl;
            for (int i=0; i != graph[vd1].sites.size(); i++)
                cout<<"site: "<< i
                    <<" num_edges to site: "<<graph[vd1].sites[i].edges.size()
                    << " age of site: " << graph[vd1].sites[i].age << endl;
            cout<<endl;
        }
#endif
	
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
        outfile<<endl;

    }

    outfile.close();
    infile.close();
}
