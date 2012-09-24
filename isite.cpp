#include<cstdio>    //fscanf can be useful if sites are read as strings
#include<cstdlib>   //for atoi
#include<time.h>    //for random number generation
#include<iostream>  //for input and output
#include<fstream>   //for file input and output
#include<utility>   //for std::pair
#include<vector>    //for std::vector
#include<map>       //for std::map
#include<cstring>   //for std::string
#include<algorithm> //for copying vectors
#include<unistd.h>  //for better random seed
#include<cassert>   //for debugging
#include<stack>     //for printing predecessors


#include<boost/graph/adjacency_list.hpp>
#include<boost/graph/connected_components.hpp>
#include<boost/property_map/property_map.hpp>
#include<boost/random.hpp>
#include<boost/generator_iterator.hpp>
#include<boost/graph/random.hpp>

using namespace std;
using namespace boost;

typedef adjacency_list_traits<listS,listS,undirectedS>::edge_descriptor edge_descriptor;
typedef adjacency_list_traits<listS,listS,undirectedS>::vertex_descriptor vertex_descriptor;

vector<int> pred;                   //shows predecessor chain of nodes
                                    //index - node's id : value - parent's id it was copied from
                                    //seed graph nodes - value = -1

struct isite
{
    string site_name;
    vector<edge_descriptor> edges;
    unsigned int age;

    isite() : age(0), edges()
    {
    }
    isite(edge_descriptor e, unsigned int a, string name) : age(a), site_name( name )
    {
        edges.push_back(e);
    }
};

struct vertexsites
{
    int vertex_id;
    vector<isite> sites;
    map<edge_descriptor,int> edgeToSite;
    map<string, int> site_name_to_index;
};

typedef adjacency_list<listS, listS, undirectedS,
                        property<vertex_index_t, int, vertexsites> > Graph;

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
typedef property_map<Graph, vertex_index_t>::type vimap;
typedef mt19937 RNGType;


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
int components(Graph& graph, vimap& indexmap)
{
    dummy_property_map p;
    int count=0;
    vertex_iterator vi, viend;
    for (tie(vi, viend)=vertices(graph); vi!=viend; vi++)
    {
        put(indexmap, *vi, count++);
    }
    return connected_components(graph, p, vertex_index_map(indexmap));
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

//For vertex index mapping
unsigned int counter;
static RNGType rng;

/*
    Given an edge and a vertex, returns the vertex on the opposite end
*/
Graph::vertex_descriptor edgeDest(const Graph::vertex_descriptor vd,
                                  const Graph::edge_descriptor ed,
                                  const Graph& graph)
{
    return (vd==source(ed, graph)) ? target(ed, graph) : source(ed, graph);
}

/*
    Completely duplicates a random node along with iSites and edges
*/
pair<Graph::vertex_descriptor,Graph::vertex_descriptor> duplicate(Graph& graph,
                                                               vimap& indexmap)
{
    Graph::vertex_descriptor parent_description = random_vertex(graph, rng);

#ifdef DEBUG
    cout << "**Duplication**" << endl;
    cout << "Parent: " << indexmap(parent_description)<<endl;
#endif

    //get a new vertex
    Graph::vertex_descriptor child_description = add_vertex(graph);
    put(indexmap, child_description, counter++);

#ifdef DEBUG
    cout << "Child: " << indexmap(child_description) << endl;
#endif

    //copy edges and isites
    int numSites = graph[parent_description].sites.size();
    for (int i=0; i < numSites; i++) //cycle through sites
    {
        isite newSite;
        newSite.age=0;

        int siteEdges = graph[parent_description].sites[i].edges.size();
        for (int j=0; j < siteEdges; j++) //cycle through edges
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

    //update pred[] for child
    //index - child's id. value - parent's id
    pred.push_back(graph[parent_description].vertex_id);
    graph[child_description].vertex_id = pred.size()-1; 

    //return parent, child
    pair<Graph::vertex_descriptor, Graph::vertex_descriptor> tmpPair;
    tmpPair = std::make_pair(parent_description, child_description);
    return tmpPair; 
}

/*
    Performs iSite algorithm
    1)Duplicate random node
    2)For each iSite, use prob_asym to choose parent or child
    3)For each edge on chosen site, use prob_loss to determine if lost
    4)Remove isolated vertices (vertex with degree=0)
*/
void duplication(Graph& graph, vimap& indexmap)
{
    uniform_real<> real_dist(0, 1);
    variate_generator<RNGType&, uniform_real<> > rand_real(rng, real_dist);
    pair<Graph::vertex_descriptor,Graph::vertex_descriptor> vertices;
    vertices=duplicate(graph, indexmap);

    //prob_asym
    int numSites = graph[vertices.first].sites.size();

#ifdef DEBUG
    cout<<"iSites: "<<numSites<<endl;
#endif

    for (int i=0; i<numSites; i++)
    {
        Graph::vertex_descriptor vertexLoss;

        double rand_res = rand_real();
        if (rand_res <= param.prob_asym) //Parent loss
        {
#ifdef DEBUG
            cout<<"\nAsymetry: Parent"<<endl;
#endif
            vertexLoss = vertices.first;
        }
        else //Child loss
        {
#ifdef DEBUG
            cout<<"\nAsymetry: Child"<<endl;
#endif
            vertexLoss = vertices.second;
        }

        //prob_loss
        int numEdges = graph[vertexLoss].sites[i].edges.size();
#ifdef DEBUG
        cout<< "Edges: " << numEdges << endl;
#endif
        for (int j=0; j<numEdges; j++)
        {
            rand_res = rand_real();
            if (rand_res <= param.prob_loss) //Edge is lost
            {
#ifdef DEBUG
                cout<<"\tLoss: Yes"<<endl;
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

                //*******Need to check if connected iSite is empty

                remove_edge(edgeLoss, graph);
                j--;
                numEdges--;
            }
#ifdef DEBUG
            else cout<<"\tLoss: No"<<endl;
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
            map<edge_descriptor,int>::iterator mit;
            for (mit=graph[vertices.first].edgeToSite.begin();
                 mit!=graph[vertices.first].edgeToSite.end();
                 mit++)
            {
                if (mit->second > i) mit->second--;
            }
            i--;
            j--;
        }
    for (int i=0; i<numSites; i++)
        if (graph[vertices.second].sites[i].edges.size()==0)
        {
            graph[vertices.second].sites.erase(
                graph[vertices.second].sites.begin() + i);
            map<edge_descriptor,int>::iterator mit;
            for (mit=graph[vertices.second].edgeToSite.begin();
                 mit!=graph[vertices.second].edgeToSite.end();
                 mit++)
            {
                if (mit->second > i) mit->second--;
            }
            i--;
            numSites--;
        }

    //If node has no iSites (and therefore no edges), delete it
    if (graph[vertices.first].sites.empty())
    {
        remove_vertex(vertices.first, graph);
    }
    if (graph[vertices.second].sites.empty())
    {
        remove_vertex(vertices.second, graph);
    }
}

/*
    Increases age of every iSite by 1
*/
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


/*
    Prints graph to standard out
*/
void printGraph(Graph& graph, vimap& indexmap)
{
    vertex_iterator vi, viend;
    for (tie(vi,viend) = vertices(graph); vi!=viend; ++vi)
    {
        //Current node
        cout<<indexmap(*vi)<<":";
        out_edge_iterator oei, oeiend;
        //All connected nodes
        for (tie(oei,oeiend)=out_edges(*vi,graph); oei!=oeiend; ++oei)
        {
            Graph::vertex_descriptor vd = target(*oei,graph);
            assert(indexmap(vd) == graph[vd].vertex_id);
            int site = graph[*vi].edgeToSite[*oei];
            int connectedSite = graph[vd].edgeToSite[*oei];
            cout<<" "<<graph[*vi].sites[site].site_name<<"->"
            <<graph[vd].sites[connectedSite].site_name<<":"<<indexmap(vd);
        }
        cout<<endl;
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
    vimap indexmap = get(vertex_index, graph);
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

    bool newNode1 = false;
    bool newNode2 = false;
    map<string,int> nodes;
    counter=0;
    string v1,v2, s1, s2, tmpStr;
    rng.seed(time(NULL)*getpid());
    
#ifdef DEBUG
    cout << "***Building graph from file***" << endl;
#endif

    //Building seed graph
    while( !(getline(infile, tmpStr)).eof() )
    {
//expected input "v1 s1:s2 v2"
//Parse input
        newNode1 = false;
        newNode2 = false;

        v1 = tmpStr.substr(0,tmpStr.find(" "));
        tmpStr.erase(0,v1.size()+1);
        s1 = tmpStr.substr(0,tmpStr.find(":"));
        tmpStr.erase(0,s1.size()+1);
        s2 = tmpStr.substr(0,tmpStr.find(" "));
        tmpStr.erase(0,s2.size()+1);
        v2 = tmpStr;

#ifdef DEBUG
        cout << v1 << ":" << s1 << "<->" << s2 << ":" << v2 << endl;
#endif

        //Create nodes if don't exist
        if(nodes[v1]==0)
        {
            newNode1 = true;
            nodes[v1]=counter+1;
            put(indexmap, add_vertex(graph), counter++);
        }
        if(nodes[v2]==0)
        {
            newNode2 = true;
            nodes[v2]=counter+1;
            put(indexmap, add_vertex(graph), counter++);
        }

        //Add edge
        Graph::vertex_descriptor vd1 = vertex(nodes[v1]-1,graph);
        Graph::vertex_descriptor vd2 = vertex(nodes[v2]-1,graph);
        Graph::edge_descriptor ed;

//set id's to the nodes - but only if they have not been added before
        if (newNode1)
        {
            pred.push_back(-1);
            graph[vd1].vertex_id = pred.size()-1;
        }
        if (newNode2)
        {
            pred.push_back(-1);
            graph[vd2].vertex_id = pred.size()-1;
        }

        bool temp_bool;
        tie(ed, temp_bool) = add_edge(vd1,vd2,graph);

        //Add iSite
        map<string, int>::iterator known_site_it;
        known_site_it = graph[vd1].site_name_to_index.find(s1);
        if (known_site_it == graph[vd1].site_name_to_index.end())
        {
            //if site has not been added to this node before
            graph[vd1].sites.push_back(isite(ed,0, s1));
            graph[vd1].edgeToSite.insert(make_pair(ed, graph[vd1].sites.size()-1));
            graph[vd1].site_name_to_index.insert(make_pair(s1, graph[vd1].sites.size()-1));
        }
        else {
            //if site is already on the node - just add edge to it
            graph[vd1].sites[known_site_it->second].edges.push_back(ed);
        }

        known_site_it = graph[vd2].site_name_to_index.find(s2);
        if (known_site_it == graph[vd2].site_name_to_index.end())
        {
            //if site has not been added to this node before
            graph[vd2].sites.push_back(isite(ed,0, s2));
            graph[vd2].edgeToSite.insert(make_pair(ed, graph[vd2].sites.size()-1));
            graph[vd2].site_name_to_index.insert(make_pair(s2, graph[vd2].sites.size()-1));
        }
        else {
            //if site is already on the node - just add edge to it
            graph[vd2].sites[known_site_it->second].edges.push_back(ed);
        }
		
    }

#ifdef DEBUG
    cout << endl << "***Original graph***" << endl;
    printGraph(graph, indexmap);
    cout << endl;
    int iterations=param.iterations;
#endif

    //Opening output file
    ofstream outfile("result");
    if (!outfile)
    {
        cerr<<"Error opening output file: result"<<endl;
        exit(1);
    }

    while (param.iterations--) //needs to encompass building seed graph, too
    {
#ifdef DEBUG
        cout << "Iteration: " << iterations-param.iterations << endl;
#endif
        while (num_vertices(graph)<param.end_order)
        {
            addAge(graph);
            duplication(graph, indexmap);
#ifdef DEBUG
            cout << endl << "***Graph during algorithm***" << endl;
            printGraph(graph, indexmap);
            cout << endl;
#endif

        }

#ifdef DEBUG
        cout << "***End graph***" << endl;
        printGraph(graph, indexmap);
        cout << endl;
#endif

#ifdef DEBUG
        cout << "***Node summary***" << endl;
        Graph::vertex_iterator vi, viend;
        for (tie(vi,viend) = vertices(graph); vi!=viend; ++vi)
        {
            Graph::vertex_descriptor vd1 = *vi; 
            cout<<"Node: " << graph[vd1].vertex_id << endl;
            for (int i=0; i != graph[vd1].sites.size(); i++)
                cout<<"\t"<<graph[vd1].sites[i].site_name
                    <<":: Age: "<<graph[vd1].sites[i].age
                    <<", Edges: "<<graph[vd1].sites[i].edges.size()<<endl;
            cout<<endl;
        }
#endif
	
        outfile<<param.prob_loss<<" ";
        outfile<<param.prob_asym<<" ";
        outfile<<num_vertices(graph)<<" ";
        outfile<<num_edges(graph)<<" ";
        int numTriangles=triangles(graph);
        int numTriples=countTriples(graph);
        outfile<<numTriangles<<" ";
        outfile<<numTriples<<" ";
        outfile<<(double)(3*numTriangles)/numTriples<<" ";
        outfile<<components(graph, indexmap);
        outfile<<endl;

    }

#ifdef DEBUG
        cout << "***Node evolution***" << endl;
        Graph::vertex_iterator vi, viend;
        for (tie(vi,viend) = vertices(graph); vi!=viend; ++vi)
        {
            stack<int> predStack;
            int curNum = graph[*vi].vertex_id;
            cout << curNum << ": ";
            while(pred[curNum] != -1) 
            {
                predStack.push(pred[curNum]);
                curNum = pred[curNum];
            }
            while (!predStack.empty())
            {
                cout << predStack.top() << "->";
                predStack.pop();
            }
            cout << graph[*vi].vertex_id << endl;
        }
#endif



    outfile.close();
    infile.close();
}


