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
#include<set>


#include<boost/graph/adjacency_list.hpp>
#include<boost/graph/connected_components.hpp>
#include<boost/property_map/property_map.hpp>
#include<boost/random.hpp>
#include<boost/generator_iterator.hpp>
#include<boost/graph/random.hpp>
#include<boost/lexical_cast.hpp>

using namespace std;
using namespace boost;

enum output_type {PRINT, NODE_SUMMARY, STATUS, EVOLUTION};

typedef adjacency_list_traits<listS,listS,undirectedS>::edge_descriptor edge_descriptor;
typedef adjacency_list_traits<listS,listS,undirectedS>::vertex_descriptor vertex_descriptor;

int nSelfLoops=0;
double asymmetry=0.0;
int nDuplications=0;

struct isite
{
    string site_name;
    vector<edge_descriptor> edges;
    unsigned int age;

    isite() : site_name(""), edges(), age(0)
    {
    }
    isite(const isite& is) : site_name(is.site_name), edges(is.edges), age(is.age)
    {
    }
    isite(edge_descriptor e, unsigned int a, string name) : site_name(name), age(a)
    {
        edges.push_back(e);
    }
};

struct vertexsites
{
    map<edge_descriptor, pair<int, int> > selfLoops;
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


Graph::vertex_descriptor edgeDest(const Graph::vertex_descriptor vd,
                                  const Graph::edge_descriptor ed,
                                  const Graph& graph);


/*
    Outputs formatted edge output to standard out
*/
void printEdge(Graph& graph,
               Graph::vertex_descriptor& vd,
               Graph::edge_descriptor& ed,
               vimap& indexmap)
{
    Graph::vertex_descriptor connectedNode = edgeDest(vd, ed, graph);
    int site = graph[vd].edgeToSite[ed];
    if (site == -1)
    {
        site = graph[vd].selfLoops[ed].first;
        int site2 = graph[vd].selfLoops[ed].second;
        cout<<" "<<graph[vd].sites[site].site_name<<"->"
            <<graph[vd].sites[site2].site_name<<":"
            <<indexmap(vd);
    }
    else
    {
        int connectedSite = graph[connectedNode].edgeToSite[ed];
        cout<<" "<<graph[vd].sites[site].site_name<<"->"
            <<graph[connectedNode].sites[connectedSite].site_name<<":"
            <<indexmap(connectedNode);
    }
}


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
    double prob_self;
    double prob_fusion;
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
    Fuses duplications of two nodes
*/
Graph::vertex_descriptor duplicate(Graph& graph, vimap& indexmap, vector<Graph::vertex_descriptor>& progenitors)
{
    progenitors.push_back(random_vertex(graph, rng));
    uniform_real<> real_dist(0.0, 1.0);
    variate_generator<RNGType&, uniform_real<> > rand_real(rng, real_dist);
    //Make 'while' to do arbitrary number of nodes
    if (rand_real() <= param.prob_fusion)
    {
        Graph::vertex_descriptor progenitor;
        while ((progenitor=random_vertex(graph, rng)) == progenitors[0]);
        progenitors.push_back(progenitor);
    }

#ifdef DEBUG
    cout << "**Duplication**" << endl;
    for (unsigned int i=0; i<progenitors.size(); ++i)
        cout << "Parent: " << indexmap(progenitors[i]) << endl;
#endif

    //get a new vertex
    Graph::vertex_descriptor progeny = add_vertex(graph);
    put(indexmap, progeny, counter++);

#ifdef DEBUG
    cout << "Child: " << indexmap(progeny) << endl;
#endif

    //copy edges and isites
    int totalSites=0;
    for (vector<Graph::vertex_descriptor>::iterator it = progenitors.begin();
            it != progenitors.end(); ++it)
    {
        int numSites = graph[*it].sites.size();
        for (int i=0; i < numSites; ++totalSites, ++i)
        {
            isite newSite;
            newSite.age = 0;
            newSite.site_name = lexical_cast<string>(indexmap(*it)) + "." + graph[*it].sites[i].site_name;

            int siteEdges = graph[*it].sites[i].edges.size();
            for (int j=0; j < siteEdges; j++) //cycle through edges
            {
                Graph::edge_descriptor original_edge;
                original_edge = graph[*it].sites[i].edges[j];
                if (edgeDest(*it, original_edge, graph) == progeny) continue;

                //self loop
                if (graph[*it].edgeToSite[original_edge] == -1)
                {
                    ++nSelfLoops;
                    int site1 = graph[*it].selfLoops[original_edge].first;
                    int site2 = graph[*it].selfLoops[original_edge].second;
                    //self loop on same node
                    if (site1 == site2)
                    {
                        Graph::edge_descriptor ed, ed_tie;
                        bool temp_bool;
                        //copy self loop
                        tie(ed, temp_bool) = add_edge(progeny, progeny, graph);

                        //create edge between parent and child
                        tie(ed_tie, temp_bool) = add_edge(*it, progeny, graph);

                        //update child site
                        newSite.edges.push_back(ed);
                        newSite.edges.push_back(ed_tie);
                        graph[progeny].edgeToSite[ed] = -1;
                        graph[progeny].selfLoops[ed] = make_pair(totalSites, totalSites);
                        graph[progeny].edgeToSite[ed_tie] = totalSites;

                        //update parent site
                        graph[*it].sites[i].edges.push_back(ed_tie);
                        graph[*it].edgeToSite[ed_tie] = i;
                    }
                    //self loop on different nodes
                    else
                    {
                        assert(0==1);
                    }
                }
                else
                {
                
                    Graph::edge_descriptor ed, ed2;
                    bool temp_bool;
                    Graph::vertex_descriptor vd = edgeDest(*it, original_edge, graph);

                    tie(ed, temp_bool) = add_edge(progeny,vd,graph);

                    //update child iSite
                    newSite.edges.push_back(ed);
                    graph[progeny].edgeToSite[ed]=totalSites;

                    //update other node's iSite
                    ed2 = graph[*it].sites[i].edges[j];
                    graph[vd].sites[graph[vd].edgeToSite[ed2]].edges.push_back(ed);
                    graph[vd].edgeToSite[ed]=graph[vd].edgeToSite[ed2];

                }
            }

            graph[progeny].sites.push_back(newSite);
        }
    }

    return progeny;
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
    ++nDuplications;
    uniform_real<> real_dist(0.0, 1.0);
    variate_generator<RNGType&, uniform_real<> > rand_real(rng, real_dist);
    vector<Graph::vertex_descriptor> progenitors;
    Graph::vertex_descriptor progeny;
    progeny=duplicate(graph, indexmap, progenitors);

    //prob_asym
    pair<int, int> actualAsym(0,0); //progenitor, progeny asymmetry
    for (unsigned int pr=0, py=0; pr < progenitors.size(); ++pr)
    {
        int numSites = graph[progenitors[pr]].sites.size();
        int site;

#ifdef DEBUG
        cout<<indexmap(progenitors[pr])<<": iSites: "<<numSites<<endl;
#endif

        for (int i=0; i<numSites; ++i, ++py)
        {
            Graph::vertex_descriptor vertexLoss;

#ifdef DEBUG
            cout<<"\niSite: "<<graph[progeny].sites[py].site_name<<endl;
#endif
            if (rand_real() <= param.prob_asym) //Parent loss
            {
#ifdef DEBUG
                cout<<"Asymmetry: Progenitor"<<endl;
#endif
                vertexLoss = progenitors[pr];
                site = i;
            }
            else //Child loss
            {
#ifdef DEBUG
                cout<<"Asymmetry: Progeny"<<endl;
#endif
                vertexLoss = progeny;
                site = py;
            }

            //prob_loss
            int numEdges = graph[vertexLoss].sites[site].edges.size();
#ifdef DEBUG
            cout<< "Edges: " << numEdges << endl;
#endif
            for (int j=0; j<numEdges; ++j)
            {
#ifdef DEBUG
                cout<<"\tLoss ";
                printEdge(graph, vertexLoss, graph[vertexLoss].sites[site].edges[j], indexmap);
#endif
                Graph::edge_descriptor edgeLoss;
                edgeLoss = graph[vertexLoss].sites[site].edges[j];
                double prob_loss;
                if (source(edgeLoss, graph) != target(edgeLoss, graph))
                    prob_loss = param.prob_loss;
                else
                    prob_loss = param.prob_self;
                if (rand_real() <= prob_loss) //Edge is lost
                {
#ifdef DEBUG
                    cout<<": Yes"<<endl;
#endif
                    vertexLoss==progenitors[pr] ? ++actualAsym.first : ++actualAsym.second;

                    //Remove from vertexLoss
                    graph[vertexLoss].edgeToSite.erase(edgeLoss);
                    graph[vertexLoss].sites[site].edges.erase(
                        graph[vertexLoss].sites[site].edges.begin()+j);

                    //Check for self loop
                    if (source(edgeLoss, graph) != target(edgeLoss, graph))
                    {
                        //Update connected vertex
                        Graph::vertex_descriptor connectedVertex;
                        connectedVertex = edgeDest(vertexLoss, edgeLoss, graph);
                        int connectedSite = graph[connectedVertex].edgeToSite[edgeLoss];
                        int connectedSiteSize = graph[connectedVertex].sites[connectedSite].edges.size();
                        int k;
                        for (k=0; k<connectedSiteSize; ++k)
                        {
                            if (graph[connectedVertex].sites[connectedSite].edges[k]==edgeLoss)
                                break;
                        }
                        graph[connectedVertex].edgeToSite.erase(edgeLoss);
                        graph[connectedVertex].sites[connectedSite].edges.erase(
                            graph[connectedVertex].sites[connectedSite].edges.begin()+k);
                    }
                    else --nSelfLoops;

                    remove_edge(edgeLoss, graph);
                    --j;
                    --numEdges;
                }
#ifdef DEBUG
                else cout<<": No"<<endl;
#endif

            }//cycle through edges

        }//cycle through iSites

        //Check progenitor for empty iSites
        for (int i=0; i<numSites; ++i)
            if (graph[progenitors[pr]].sites[i].edges.size()==0)
            {
                graph[progenitors[pr]].sites.erase(
                    graph[progenitors[pr]].sites.begin() + i);
                map<edge_descriptor,int>::iterator mit;
                for (mit=graph[progenitors[pr]].edgeToSite.begin();
                     mit!=graph[progenitors[pr]].edgeToSite.end();
                     mit++)
                {
                    if (mit->second > i) mit->second--;
                }
                --i;
                --numSites;
            }
        //If progenitor has no iSites (and therefore no edges), delete it
        if (graph[progenitors[pr]].sites.empty())
        {
            remove_vertex(progenitors[pr], graph);
        }
    }//cycle through progenitors

    //Check progeny for empty iSites
    for (int i=0, size=graph[progeny].sites.size(); i<size; ++i)
        if (graph[progeny].sites[i].edges.size()==0)
        {
            graph[progeny].sites.erase(
                graph[progeny].sites.begin() + i);
            map<edge_descriptor,int>::iterator mit;
            for (mit=graph[progeny].edgeToSite.begin();
                 mit!=graph[progeny].edgeToSite.end();
                 mit++)
            {
                if (mit->second > i) mit->second--;
            }
            --i;
            --size;
        }
    //If progeny has no iSites (and therefore no edges), delete it
    if (graph[progeny].sites.empty())
    {
        remove_vertex(progeny, graph);
    }
    //Update actual asymmetry
    double thisAsym;
    if (actualAsym.first+actualAsym.second != 0)
        if (actualAsym.first > actualAsym.second)
            thisAsym = (double)actualAsym.first / ((double)actualAsym.first+(double)actualAsym.second);
        else
            thisAsym = (double)actualAsym.second / ((double)actualAsym.first+(double)actualAsym.second);
    else
        thisAsym = .5;
    asymmetry = asymmetry*(((double)nDuplications-1.0)/(double)nDuplications)+thisAsym*(1.0/(double)nDuplications);
}


/*
    Performs iSite algorithm
    1)Duplicate random node
    2)For each iSite, use prob_asym to choose parent or child
    3)For each edge on chosen site, use prob_loss to determine if lost
    4)Remove isolated vertices (vertex with degree=0)
*/
void oduplication(Graph& graph, vimap& indexmap)
{
    ++nDuplications;
    uniform_real<> real_dist(0.0, 1.0);
    variate_generator<RNGType&, uniform_real<> > rand_real(rng, real_dist);
    vector<Graph::vertex_descriptor> progenitors;
    Graph::vertex_descriptor progeny;
    progeny=duplicate(graph, indexmap, progenitors);

    //prob_asym
    pair<int,int> actualAsym(0,0); //progenitor,progeny asymmetry
    int numSites = graph[progenitors[0]].sites.size();

#ifdef DEBUG
    cout<<"iSites: "<<numSites<<endl;
#endif

    for (int i=0; i<numSites; ++i)
    {
        Graph::vertex_descriptor vertexLoss;

        //double rand_res = rand_real();
#ifdef DEBUG
        cout<<"\niSite: "<<graph[progenitors[0]].sites[i].site_name<<endl;
#endif
        if (rand_real() <= param.prob_asym) //Parent loss
        {
#ifdef DEBUG
            cout<<"Asymmetry: Progenitor"<<endl;
#endif
            vertexLoss = progenitors[0];
        }
        else //Child loss
        {
#ifdef DEBUG
            cout<<"Asymmetry: Progeny"<<endl;
#endif
            vertexLoss = progeny;
        }

        //prob_loss
        int numEdges = graph[vertexLoss].sites[i].edges.size();
#ifdef DEBUG
        cout<< "Edges: " << numEdges << endl;
#endif
        for (int j=0; j<numEdges; ++j)
        {
#ifdef DEBUG
            cout<<"\tLoss ";
            printEdge(graph, vertexLoss, graph[vertexLoss].sites[i].edges[j], indexmap);
#endif
            Graph::edge_descriptor edgeLoss;
            edgeLoss = graph[vertexLoss].sites[i].edges[j];
            double prob_loss;
            if (source(edgeLoss, graph) != target(edgeLoss, graph))
                prob_loss = param.prob_loss;
            else
                prob_loss = param.prob_self;
            //rand_res = rand_real();
            if (rand_real() <= prob_loss) //Edge is lost
            {
#ifdef DEBUG
                cout<<": Yes"<<endl;
#endif
                vertexLoss==progenitors[0] ? ++actualAsym.first : ++actualAsym.second;

                //Remove from vertexLoss
                graph[vertexLoss].edgeToSite.erase(edgeLoss);
                graph[vertexLoss].sites[i].edges.erase(
                    graph[vertexLoss].sites[i].edges.begin()+j);

                //Check for self loop
                if (source(edgeLoss, graph) != target(edgeLoss, graph))
                {
                    //Update connected vertex
                    Graph::vertex_descriptor connectedVertex;
                    connectedVertex = edgeDest(vertexLoss, edgeLoss, graph);
                    int connectedSite = graph[connectedVertex].edgeToSite[edgeLoss];
                    int connectedSiteSize = graph[connectedVertex].sites[connectedSite].edges.size();
                    int k;
                    for (k=0; k<connectedSiteSize; ++k)
                    {
                        if (graph[connectedVertex].sites[connectedSite].edges[k]==edgeLoss)
                            break;
                    }
                    graph[connectedVertex].edgeToSite.erase(edgeLoss);
                    graph[connectedVertex].sites[connectedSite].edges.erase(
                        graph[connectedVertex].sites[connectedSite].edges.begin()+k);
                }
                else --nSelfLoops;

                remove_edge(edgeLoss, graph);
                --j;
                --numEdges;
            }
#ifdef DEBUG
            else cout<<": No"<<endl;
#endif

        }//cycle through edges

    }//cycle through iSites

    //Check parent and child for empty iSites
    int j=numSites;
    for (int i=0; i<j; ++i)
        if (graph[progenitors[0]].sites[i].edges.size()==0)
        {
            graph[progenitors[0]].sites.erase(
                graph[progenitors[0]].sites.begin() + i);
            map<edge_descriptor,int>::iterator mit;
            for (mit=graph[progenitors[0]].edgeToSite.begin();
                 mit!=graph[progenitors[0]].edgeToSite.end();
                 mit++)
            {
                if (mit->second > i) mit->second--;
            }
            --i;
            --j;
        }
    for (int i=0; i<numSites; ++i)
        if (graph[progeny].sites[i].edges.size()==0)
        {
            graph[progeny].sites.erase(
                graph[progeny].sites.begin() + i);
            map<edge_descriptor,int>::iterator mit;
            for (mit=graph[progeny].edgeToSite.begin();
                 mit!=graph[progeny].edgeToSite.end();
                 mit++)
            {
                if (mit->second > i) mit->second--;
            }
            --i;
            --numSites;
        }

    //If node has no iSites (and therefore no edges), delete it
    if (graph[progenitors[0]].sites.empty())
    {
        remove_vertex(progenitors[0], graph);
    }
    else if (graph[progeny].sites.empty())
    {
        remove_vertex(progeny, graph);
    }

    //Update actual asymmetry
    double thisAsym;
    if (actualAsym.first+actualAsym.second != 0)
        if (actualAsym.first > actualAsym.second)
            thisAsym = (double)actualAsym.first / ((double)actualAsym.first+(double)actualAsym.second);
        else
            thisAsym = (double)actualAsym.second / ((double)actualAsym.first+(double)actualAsym.second);
    else
        thisAsym = .5;
    asymmetry = asymmetry*(((double)nDuplications-1.0)/(double)nDuplications)+thisAsym*(1.0/(double)nDuplications);
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
    vector<edge_descriptor> printed_self_loops;
    for (tie(vi,viend) = vertices(graph); vi!=viend; ++vi)
    {
        //Current node
        cout<<indexmap(*vi)<<":";
        out_edge_iterator oei, oeiend;
        //All connected nodes
        for (tie(oei,oeiend)=out_edges(*vi,graph); oei!=oeiend; ++oei)
        {
            Graph::vertex_descriptor vd = target(*oei,graph);
            int site = graph[*vi].edgeToSite[*oei];
            if (site == -1)
            {
                bool already_printed = false;
                for (vector<edge_descriptor>::iterator i=printed_self_loops.begin();
                     i!=printed_self_loops.end(); i++)
                {
                    if (*i == *oei)
                    {
                        already_printed = true;
                        printed_self_loops.erase(i);
                        break;
                    }
                }

                if (!already_printed)
                {
                    site = graph[*vi].selfLoops[*oei].first;
                    int site2 = graph[*vi].selfLoops[*oei].second;
                    cout<<" "<<graph[*vi].sites[site].site_name<<"->"
                        <<graph[*vi].sites[site2].site_name<<":"
                        <<indexmap(vd);
                    printed_self_loops.push_back(*oei);
                }
            }
            else
            {
                int connectedSite = graph[vd].edgeToSite[*oei];
                cout<<" "<<graph[*vi].sites[site].site_name<<"->"
                    <<graph[vd].sites[connectedSite].site_name<<":"
                    <<indexmap(vd);
            }
        }
        printed_self_loops.clear();
        cout<<endl;
    }
}

//enum output_type {print, node_summary, status};
void output_info(Graph& graph, vimap& indexmap, const string& label, output_type type)
{
#ifdef DEBUG
    Graph::vertex_iterator vi, viend;
    switch(type)
    {
        case PRINT:
            cout << endl << label << endl;
            printGraph(graph, indexmap);
            cout << endl;
        break;

        case NODE_SUMMARY:
            cout << endl << label << endl;
            for (tie(vi,viend) = vertices(graph); vi!=viend; ++vi)
            {
                Graph::vertex_descriptor vd1 = *vi; 
                cout<<"Node: " << indexmap(vd1) << endl;
                for (unsigned int i=0; i != graph[vd1].sites.size(); i++)
                    cout<<"\t"<<graph[vd1].sites[i].site_name
                        <<":: Age: "<<graph[vd1].sites[i].age
                        <<", Edges: "<<graph[vd1].sites[i].edges.size()<<endl;
                cout<<endl;
            }
        break;
        
        default:

        break;
    }
#endif

#ifdef NDEBUG
    switch(type)
    {
        case STATUS:
            if (num_vertices(graph)%100==0)
            {
                cout<<".";
                cout.flush();
            }
        break;

        default:

        break;
    }
#endif
}


int main(int argc, char* argv[])
{
    if (argc!=8)
    {
        cerr<<"Usage: ./iSite <seed-graph> <probability of subfunctionalization> "
                "<probability of assymetry> <probability of homomeric subfunctionalization> "
                "<probability of fusion> <end order> <iterations>"<<endl;
        exit(1);
    }

    param.infile=argv[1];               //Seed graph
    param.prob_loss=atof(argv[2]);      //Probability of loss of redundancy (subfunctionalization)
    param.prob_asym=atof(argv[3]);      //Probability of assymetry
    param.prob_self=atof(argv[4]);      //Probability of homomeric subfunctionalization
    param.prob_fusion=atof(argv[5]);    //Probability of fusion
    param.end_order=atoi(argv[6]);      //Order at which to stop
    param.iterations=atoi(argv[7]);     //Number of times to run algorithm
    

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
    if (param.prob_self<0.0 || param.prob_self>1.0)
    {
        cerr<<"Error: Probability of homomeric subfunctionalization must be between 0 and 1"<<endl;
        exit(1);
    }
    if (param.prob_fusion<0.0 || param.prob_self>1.0)
    {
        cerr<<"Error: Probability of fusion must be between 0 and 1"<<endl;
        exit(1);
    }
    if (param.iterations<0)
    {
        cerr<<"Error: Number of iterations must be a positive integer"<<endl;
        exit(1);
    }

    /*******Iteration***********/
#ifdef DEBUG
    int iterations=param.iterations;
#endif
    rng.seed(time(NULL)*getpid());

    //Opening output file
    ofstream outfile("result");
    if (!outfile)
    {
        cerr<<"Error opening output file: result"<<endl;
        exit(1);
    }
    outfile << "subfuncProb asymmetry selfloopLoss actualAsymmetry selfloops order size tris trips CC numComponents" << endl;

    while (param.iterations--) //needs to encompass building seed graph, too
    {

        //Opening input file
        ifstream infile(param.infile.c_str());
        if (!infile)
        {
            cerr<<"Error opening specified seed graph!"<<endl;
            exit(1);
        }

        Graph graph;                        //Protein network
        vimap indexmap = get(vertex_index, graph);
        bool newNode1 = false;
        bool newNode2 = false;
        map<string,int> nodes;
        counter=0;
        string v1,v2, s1, s2, tmpStr;
        nSelfLoops=0;
        asymmetry=0.0;
        nDuplications=0;
                
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

            //strtok might be cleaner
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


            bool temp_bool;
            tie(ed, temp_bool) = add_edge(vd1,vd2,graph);
            pair<int,int> edge_site_indices;

            //Add iSite
            map<string, int>::iterator known_site_it;
            known_site_it = graph[vd1].site_name_to_index.find(s1);
                //if site has not been added to this node before
            if (known_site_it == graph[vd1].site_name_to_index.end())
            {
                    //make new iSite with new edge in it
                graph[vd1].sites.push_back(isite(ed,0, s1));
                    //set edgeToSite for new edge
                edge_site_indices.first = graph[vd1].sites.size()-1;
                graph[vd1].edgeToSite.insert(make_pair(ed, graph[vd1].sites.size()-1));
                    //set site_name_to_index for new iSite
                graph[vd1].site_name_to_index.insert(make_pair(s1, graph[vd1].sites.size()-1));
            }
                //if site is already on the node - just add edge to it
            else
            {
                    //insert new edge in iSite
                graph[vd1].sites[known_site_it->second].edges.push_back(ed);
                    //set edgeToSite for new edge
                edge_site_indices.first = known_site_it->second;
                graph[vd1].edgeToSite.insert(make_pair(ed, known_site_it->second));
            }

            if (!((v1==v2) && (s1==s2))) //check for self loops and same iSite
            {
                known_site_it = graph[vd2].site_name_to_index.find(s2);
                    //if site has not been added to this node before
                if (known_site_it == graph[vd2].site_name_to_index.end())
                {
                        //make new iSite with new edge in it
                    graph[vd2].sites.push_back(isite(ed,0, s2));
                        //set edgeToSite for new edge
                    edge_site_indices.second = graph[vd2].sites.size()-1;
                    graph[vd2].edgeToSite.insert(make_pair(ed, graph[vd2].sites.size()-1));
                        //set site_name_to_index for new iSite
                    graph[vd2].site_name_to_index.insert(make_pair(s2, graph[vd2].sites.size()-1));
                }
                    //if site is already on the node - just add edge to it
                else
                {
                        //insert new edge in iSite
                    graph[vd2].sites[known_site_it->second].edges.push_back(ed);
                        //set edgeToSite for new edge
                    edge_site_indices.second = known_site_it->second;
                    graph[vd2].edgeToSite.insert(make_pair(ed, known_site_it->second));
                }
            }
            else
            {
                edge_site_indices.second = edge_site_indices.first;
            }

            //self loops
            if (v1==v2)
            {
                ++nSelfLoops;
                graph[vd1].edgeToSite[ed] = -1;
                graph[vd1].selfLoops.insert(make_pair(ed, edge_site_indices));
            }
            
        }

        output_info(graph, indexmap, "***Original graph***", PRINT); 

#ifdef DEBUG
        cout << "Iteration: " << iterations-param.iterations << endl;
#endif
#ifdef NDEBUG
        cout << "Working";
        cout.flush();
#endif
        while (num_vertices(graph)<param.end_order)
        {
            addAge(graph);
            duplication(graph, indexmap);

            output_info(graph, indexmap, "***Graph during algorithm***", PRINT); 
            output_info(graph, indexmap, "", STATUS); 
        }

        output_info(graph, indexmap, "***End Graph***", PRINT); 
        output_info(graph, indexmap, "***Node Summary***", NODE_SUMMARY); 

#ifdef NDEBUG
        cout<<endl;
        cout<<"Generating results"<<endl;
#endif
	
        outfile<<param.prob_loss<<" ";
        outfile<<param.prob_asym<<" ";
        outfile<<param.prob_self<<" ";
        outfile<<setprecision(3)<<asymmetry<<" ";
        outfile<<nSelfLoops<<" ";
        outfile<<num_vertices(graph)<<" ";
        outfile<<num_edges(graph)<<" ";
        int numTriangles=triangles(graph);
        int numTriples=countTriples(graph);
        outfile<<numTriangles<<" ";
        outfile<<numTriples<<" ";
        outfile<<setprecision(6)<<(double)(3*numTriangles)/numTriples<<" ";
        outfile<<components(graph, indexmap);
        outfile<<endl;

        output_info(graph, indexmap, "***Node Evolution***", EVOLUTION); 

        infile.close();
    }
    outfile.close();
}


