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
#include<sys/stat.h>
#include<sys/types.h>
#include<errno.h>
#include<iomanip>


#include<boost/graph/connected_components.hpp>
#include<boost/lexical_cast.hpp>

#include "Random.hpp"
#include "graph_properties.hpp"

using namespace std;
using namespace boost;

enum output_type {PRINT, NODE_SUMMARY, STATUS, EVOLUTION};

int nSelfLoops=0;
double asymmetry=0.0;
int nDuplications=0;


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
        cout<<" "<<site<<"->"<<site2<<":"<<indexmap(vd);
    }
    else
    {
        int connectedSite = graph[connectedNode].edgeToSite[ed];
        cout<<" "<<site<<"->"<<connectedSite<<":"<<indexmap(connectedNode);
    }
}


/*
    Calculates iSite distribution
*/
void isite_distribution(Graph& graph, vector<int>& vmap)
{
    vertex_iterator vi, viend;
    for (tie(vi, viend)=vertices(graph); vi!=viend; ++vi)
    {
        int nsites=0;
        vertexsites& node = graph[*vi];
        for (int i=0; i<node.sites.size(); ++i)
            if (!node.sites[i].edges.empty())
                ++nsites;

        if (nsites+1 > vmap.size())
            vmap.resize(nsites+1);

        ++vmap[nsites];
    }
}


/*
    Prints iSite distribution
*/
void print_distribution(vector<int>& vmap)
{
    for (int i=0; i<vmap.size(); ++i)
        cout<<i<<": "<<vmap[i]<<endl;
}


/*
    Returns the largest component of a graph
*/
Graph get_large_component(Graph& graph, vimap& indexmap)
{
    int count=0;
    vertex_iterator vi, viend;
    for (tie(vi, viend)=vertices(graph); vi!=viend; ++vi)
    {
        put(indexmap, *vi, count++);
    }
    //vector<Graph::vertex_descriptor> component(count);
    vimap component;
    connected_components(graph, component, vertex_index_map(indexmap));
    vector<int> component_counts(count, 0);
    /*
    for (int i=0; i<count; ++i)
        component_counts[component[i]]++;
    */
    vertex_iterator viter, viter_end;
    for (tie(viter, viter_end) = vertices(graph); viter != viter_end; ++viter)
        component_counts[component(*viter)]++;

    //Determine largest component
    int lsize=0;
    int whichcomponent = 0;
    for (int i=0; component_counts[i]>0; ++i)
        if (component_counts[i]>lsize)
        {
            whichcomponent = i;
            lsize = component_counts[i];
        }

    //Build largest component subgraph
    Graph lgcomponent;
    map<Graph::vertex_descriptor, Graph::vertex_descriptor> nodeToIndex;
    for (tie(viter, viter_end) = vertices(graph); viter != viter_end; ++viter)
        if (component(*viter) == whichcomponent)
        {
            Graph::vertex_descriptor node = add_vertex(lgcomponent);
            nodeToIndex[*viter] = node;
        }

    //Add edges
    edge_iterator eiter, eiter_end;
    tie(eiter,eiter_end) = edges(graph);
    for (; eiter != eiter_end; ++eiter)
    {
        Graph::vertex_descriptor src = source(*eiter, graph);
        Graph::vertex_descriptor tgt = target(*eiter, graph);
        if (nodeToIndex.count(src)==1)
            add_edge(nodeToIndex[src], nodeToIndex[tgt], lgcomponent);
    }
    return lgcomponent;
}


/*
    Removes parallel edges and self loops
*/
void simplify(Graph& graph)
{
    vertex_iterator vi, viend;
    AdjIter adji, adjiend;

    for (tie(vi, viend)=vertices(graph); vi!=viend; ++vi)
    {
        //Remove homomeric interactions
        while (edge(*vi, *vi, graph).second)
        {
            Graph::edge_descriptor ed = edge(*vi, *vi, graph).first;
            vertexsites& vs = graph[*vi];
            int site = vs.selfLoops[ed].first;
            vector<edge_descriptor>::iterator pos;
            for (pos=vs.sites[site].edges.begin(); *pos!=ed; ++pos);
            vs.sites[site].edges.erase(pos);
            remove_edge(ed, graph);
        }

        //Remove parallel interactions
        for (tie(adji, adjiend)=adjacent_vertices(*vi, graph); adji!=adjiend; ++adji)
        {
            while (edge(*vi, *adji, graph).first != edge(*adji, *vi, graph).first)
            {
                remove_edge(*vi, *adji, graph);

            }
        }
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

    for (tie(vi, viend)=vertices(graph); vi!=viend; ++vi)
        numTriples+=out_degree(*vi, graph)*(out_degree(*vi,graph)-1)/2;
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


/*
    Determines if node is isolated (only self loops)
*/
bool isolated(Graph::vertex_descriptor vd, Graph& graph)
{
    AdjIter adji, adjiend;
    for (tie(adji, adjiend) = adjacent_vertices(vd, graph);
         adji != adjiend; ++adji)
    {
        if (*adji != vd) return false;
    }
    return true;
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
static Random* rnd;

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
    progenitors.push_back(rnd->random_vertex(graph));
    //Make 'while' to do arbitrary number of nodes NOT IMPLEMENTED
    if (rnd->rand() < param.prob_fusion)
    {
        Graph::vertex_descriptor progenitor;
        while ((progenitor=rnd->random_vertex(graph)) == progenitors[0]);
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
    vertexsites& pgyref = graph[progeny];

#ifdef DEBUG
    cout << "Child: " << indexmap(progeny) << endl;
#endif

    //copy edges and isites
    int totalSites=0;
    for (vector<Graph::vertex_descriptor>::iterator it = progenitors.begin();
            it != progenitors.end(); ++it)
    {
        vertexsites& vref = graph[*it];
        int numSites = vref.sites.size();
        for (int i=0; i < numSites; ++totalSites, ++i)
        {
            isite newSite;
            newSite.age = 0;
            newSite.homomeric = false;

            int siteEdges = vref.sites[i].edges.size();
            for (int j=0; j < siteEdges; j++) //cycle through edges
            {
                Graph::edge_descriptor original_edge;
                original_edge = vref.sites[i].edges[j];
                if (edgeDest(*it, original_edge, graph) == progeny) continue;

                //self loop
                if (vref.edgeToSite[original_edge] == -1)
                {
                    newSite.homomeric = true;
                    ++nSelfLoops;
                    int site1 = vref.selfLoops[original_edge].first;
                    int site2 = vref.selfLoops[original_edge].second;
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
                        pgyref.edgeToSite[ed] = -1;
                        pgyref.selfLoops[ed] = make_pair(totalSites, totalSites);
                        pgyref.edgeToSite[ed_tie] = totalSites;

                        //update parent site
                        vref.sites[i].edges.push_back(ed_tie);
                        vref.edgeToSite[ed_tie] = i;
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
                    vertexsites& vdref = graph[vd];

                    tie(ed, temp_bool) = add_edge(progeny,vd,graph);

                    //update child iSite
                    newSite.edges.push_back(ed);
                    pgyref.edgeToSite[ed]=totalSites;

                    //update other node's iSite
                    ed2 = vref.sites[i].edges[j];
                    vdref.sites[vdref.edgeToSite[ed2]].edges.push_back(ed);
                    vdref.edgeToSite[ed]=vdref.edgeToSite[ed2];

                }
            }

            pgyref.sites.push_back(newSite);
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
    vector<Graph::vertex_descriptor> progenitors;
    Graph::vertex_descriptor progeny;
    progeny=duplicate(graph, indexmap, progenitors);
    vertexsites& pgyref = graph[progeny];

    //prob_asym
    pair<int, int> actualAsym(0,0); //progenitor, progeny asymmetry
    for (unsigned int pr=0, py=0; pr < progenitors.size(); ++pr)
    {
        vertexsites& pgrref = graph[progenitors[pr]];
        int numSites = pgrref.sites.size();
        int site;

#ifdef DEBUG
        cout<<indexmap(progenitors[pr])<<": iSites: "<<numSites<<endl;
#endif

        for (int i=0; i<numSites; ++i, ++py)
        {
            Graph::vertex_descriptor vertexLoss, notLost;

#ifdef DEBUG
            //cout<<"\niSite: "<<graph[progeny].sites[py].site_name<<endl;
#endif
            if (rnd->rand() < param.prob_asym) //Parent loss
            {
#ifdef DEBUG
                cout<<"Asymmetry: Progenitor"<<endl;
#endif
                vertexLoss = progenitors[pr];
                notLost = progeny;
                site = i;
            }
            else //Child loss
            {
#ifdef DEBUG
                cout<<"Asymmetry: Progeny"<<endl;
#endif
                vertexLoss = progeny;
                notLost = progenitors[pr];
                site = py;
            }
            vertexsites& vlref = graph[vertexLoss];

            //prob_loss
            int numEdges = vlref.sites[site].edges.size();
#ifdef DEBUG
            cout<< "Edges: " << numEdges << endl;
#endif
            for (int j=0; j<numEdges; ++j)
            {
#ifdef DEBUG
                cout<<"\tLoss ";
                printEdge(graph, vertexLoss, vlref.sites[site].edges[j], indexmap);
#endif
                Graph::edge_descriptor edgeLoss;
                edgeLoss = vlref.sites[site].edges[j];
                double prob_loss;
                /*
                if (source(edgeLoss, graph) == target(edgeLoss, graph))
                    prob_loss = param.prob_self;
                else if (edgeDest(vertexLoss, edgeLoss, graph) == notLost)
                    prob_loss = param.prob_self;
                else
                    prob_loss = param.prob_loss;
                */
                if (vlref.sites[site].homomeric)
                    prob_loss = param.prob_self;
                else
                    prob_loss = param.prob_loss;
                if (rnd->rand() < prob_loss) //Edge is lost
                {
#ifdef DEBUG
                    cout<<": Yes"<<endl;
#endif
                    vertexLoss==progenitors[pr] ? ++actualAsym.first : ++actualAsym.second;

                    //Remove from vertexLoss
                    vlref.edgeToSite.erase(edgeLoss);
                    vlref.sites[site].edges.erase(
                        vlref.sites[site].edges.begin()+j);

                    //Check for self loop
                    if (source(edgeLoss, graph) != target(edgeLoss, graph))
                    {
                        //Update connected vertex
                        Graph::vertex_descriptor connectedVertex;
                        connectedVertex = edgeDest(vertexLoss, edgeLoss, graph);
                        vertexsites& cvref = graph[connectedVertex];
                        int connectedSite = cvref.edgeToSite[edgeLoss];
                        int connectedSiteSize = cvref.sites[connectedSite].edges.size();
                        int k;
                        for (k=0; k<connectedSiteSize; ++k)
                        {
                            if (cvref.sites[connectedSite].edges[k]==edgeLoss)
                                break;
                        }
                        cvref.edgeToSite.erase(edgeLoss);
                        cvref.sites[connectedSite].edges.erase(
                            cvref.sites[connectedSite].edges.begin()+k);
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

/*
        //Check progenitor for empty iSites
        for (int i=0; i<numSites; ++i)
            if (pgrref.sites[i].edges.size()==0)
            {
                pgrref.sites.erase(
                        pgrref.sites.begin() + i);
                map<edge_descriptor,int>::iterator mit;
                for (mit=pgrref.edgeToSite.begin();
                     mit!=pgrref.edgeToSite.end();
                     mit++)
                {
                    if (mit->second > i) mit->second--;
                }
                --i;
                --numSites;
            }
*/
        //If progenitor has no iSites (and therefore no edges), delete it
        if (isolated(progenitors[pr], graph))
        {
            //clear_vertex(progenitors[pr], graph);
            remove_vertex(progenitors[pr], graph);
        }
    }//cycle through progenitors

/*
    //Check progeny for empty iSites
    for (int i=0, size=pgyref.sites.size(); i<size; ++i)
        if (pgyref.sites[i].edges.size()==0)
        {
            pgyref.sites.erase(
                pgyref.sites.begin() + i);
            map<edge_descriptor,int>::iterator mit;
            for (mit=pgyref.edgeToSite.begin();
                 mit!=pgyref.edgeToSite.end();
                 mit++)
            {
                if (mit->second > i) mit->second--;
            }
            --i;
            --size;
        }
*/
    //If progeny has no iSites (and therefore no edges), delete it
    if (isolated(progeny, graph))
    {
        //clear_vertex(progeny, graph);
        remove_vertex(progeny, graph);
        --nDuplications;
        return;
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
    Determines sites/edges to be duplicated/removed

    Returns false if duplication should be ignored (100% asym for progeny and 100% subfunc)
*/
bool determineInteractions(Graph& graph,
                           vimap& indexmap,
                           vector<Graph::vertex_descriptor>& progenitors,
                           vector<bool>& progenitorLoss,
                           vector<bool>& progenyAdd,
                           vector<pair<bool, bool> >& selfLoops)
{
    pair<int, int> actualAsym(0,0); //progenitor, progeny asymmetry

    //Choose random vertices
    progenitors.push_back(rnd->random_vertex(graph));
    int progenitorDegree = out_degree(progenitors[0], graph);
    //Make 'while' to do arbitrary number of progenitors NOT IMPLEMENTED
    if (rnd->rand() <= param.prob_fusion)
    {
        Graph::vertex_descriptor progenitor;
        //Needs to be modified to work with arbitrary number
        while ((progenitor=rnd->random_vertex(graph)) == progenitors[0]);
        progenitors.push_back(progenitor);
        progenitorDegree += out_degree(progenitor, graph);
    }

#ifdef DEBUG
    cout << "**Duplication**" << endl;
    for (unsigned int i=0; i<progenitors.size(); ++i)
        cout << "Parent: " << indexmap(progenitors[i]) << endl;
#endif

    progenitorLoss.resize(progenitorDegree);
    progenyAdd.resize(progenitorDegree);
    int progenyDegree = 0;
    bool onlySelfLoops = true;


    //Cycle through progenitors
    for (int prog=0, x=0; prog<progenitors.size(); ++prog)
    {
        //Cycle through iSites
        int numSites = graph[progenitors[prog]].sites.size();
#ifdef DEBUG
        cout<<indexmap(progenitors[prog])<<": iSites: "<<numSites<<endl;
#endif
        for (int site=0; site<numSites; ++site)
        {
            //prob_asym
            bool asym_progenitor = (rnd->rand() <= param.prob_asym);

#ifdef DEBUG
            //cout<<"\niSite: "<<graph[progenitors[prog]].sites[site].site_name<<endl;
            cout << "Asymmetry: " << (asym_progenitor ? "Progenitor" : "Progeny") << endl;
#endif

            //Cycle through edges
            int numEdges = graph[progenitors[prog]].sites[site].edges.size();
#ifdef DEBUG
            cout<< "Edges: " << numEdges << endl;
#endif
            for (int edge=0; edge<numEdges; ++edge, ++x)
            {
                double prob_loss;
                //Check for self loop
                if (graph[progenitors[prog]].edgeToSite[graph[progenitors[prog]].sites[site].edges[edge]]==-1)
                {
                    prob_loss = param.prob_self;

                    //Probability of keeping homomeric interactions
                    selfLoops.push_back(pair<bool,bool>((rnd->rand()<=prob_loss),(rnd->rand()<=prob_loss)));
                    //Paralogous interaction is determined like other heteromeric interactions
                        //but with selfLoopLoss probability
                }
                else //heteromeric
                    prob_loss = param.prob_loss;
//Need to check if edge between fused nodes
                if (rnd->rand() <= prob_loss)
                {
                    if (asym_progenitor) //loss progenitor
                    {
                        progenitorLoss[x] = true;
                        progenyAdd[x] = true;
                        ++actualAsym.first;
                        ++progenyDegree;
                        if (prob_loss == param.prob_self)
                            onlySelfLoops = false;
                    }
                    else //loss progeny
                    {
                        progenitorLoss[x] = false;
                        progenyAdd[x] = false;
                        ++actualAsym.second;
                    }
                }
                else //no loss
                {
                    progenitorLoss[x] = false;
                    progenyAdd[x] = true;
                    ++progenyDegree;
                    onlySelfLoops = false;
                } //if (loss)
            } //Cycle through edges
        } //Cycle through iSites
    } //Cycle through progenitors

    //Modify this to keep isolated nodes with self loops
    if (progenyDegree == 0 || onlySelfLoops)
        return false;

    //Duplication is happening
    ++nDuplications;

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

    return true;
}


/*
    Performs iSite algorithm
    1) Determine sites/edges to be duplicated/removed
    2) Duplicate /remove edges
*/
void preemptiveDuplication(Graph& graph, vimap& indexmap)
{
    vector<Graph::vertex_descriptor> progenitors;
    vector<bool> progenitorLoss;
    vector<bool> progenyAdd;
    vector<pair<bool,bool> > selfLoops; //progenitorLoss, progenyAdd

    //Determine sites/edges to be duplicated/removed
    if (!determineInteractions(graph, indexmap, progenitors, progenitorLoss, progenyAdd, selfLoops))
        return;

    //Add new vertex
    Graph::vertex_descriptor progeny = add_vertex(graph);
    put(indexmap, progeny, counter++);

#ifdef DEBUG
    cout << "Child: " << indexmap(progeny) << endl;
#endif

    /*
     *Add edges to progeny
     */
    //Cycle through progenitors
    for (int prog=0, x=0, progeny_site=0, self_loop_it=0; prog<progenitors.size(); ++prog)
    {
        Graph::vertex_descriptor pr = progenitors[prog];
        //Cycle through iSites
        int numSites = graph[pr].sites.size();
        for (int site=0; site<numSites; ++site, ++progeny_site)
        {
            isite newSite;
            newSite.age = 0;
            /*newSite.site_name = *lexical_cast<string>(indexmap(progenitors[prog]))
                + "."
                + graph[progenitors[prog]].sites[site].site_name;*/

            //Cycle through edges
            //int numEdges = graph[progenitors[prog]].sites[site].edges.size();
            for (int edge=0; edge<graph[pr].sites[site].edges.size(); ++edge, ++x)
            {
                Graph::edge_descriptor original_edge = graph[pr].sites[site].edges[edge];

                //Edge is new
                if (edgeDest(pr, original_edge, graph) == progeny)
                {
                    --x;
                    continue;
                }

                //Self loop
                if (graph[pr].edgeToSite[original_edge] == -1)
                {
                    //Self loop needs to be added to progeny
                    if (selfLoops[self_loop_it].second)
                    {
                        ++nSelfLoops;
                        int site1 = graph[pr].selfLoops[original_edge].first;
                        int site2 = graph[pr].selfLoops[original_edge].second;
                        
                        //Self loop on same iSite
                        if (site1 == site2)
                        {
                            //Add self loop
                            Graph::edge_descriptor ed;
                            bool temp_bool;
                            tie(ed, temp_bool) = add_edge(progeny, progeny, graph);

                            //Add to iSite
                            newSite.edges.push_back(ed);
                            graph[progeny].edgeToSite[ed] = -1;
                            graph[progeny].selfLoops[ed] = make_pair(progeny_site, progeny_site);
                        }
                        else //Self loop on different iSites NOT IMPLEMENTED
                        {
                            assert(0==1);
                        }
                    }

                    //Self loop needs to be removed from progenitor
                    if (selfLoops[self_loop_it++].first)
                    {
                        --nSelfLoops;
                        graph[pr].selfLoops.erase(original_edge);
                        graph[pr].edgeToSite.erase(original_edge);
                        graph[pr].sites[site].edges.erase(
                            graph[pr].sites[site].edges.begin()+edge);

                        remove_edge(original_edge, graph);
                    }

                    //Add paralogous interaction
                    if (progenyAdd[x] && !progenitorLoss[x])
                    {
                        Graph::edge_descriptor ed;
                        bool temp_bool;
                        tie(ed, temp_bool) = add_edge(pr, progeny, graph);

                        //Add to iSite
                        newSite.edges.push_back(ed);
                        graph[progeny].edgeToSite[ed] = progeny_site;

                        //Update progenitor iSite
                        graph[pr].sites[site].edges.push_back(ed);
                        graph[pr].edgeToSite[ed] = site;
                    }
                    continue;
                }

                //Add edge to progeny
                if (progenyAdd[x])
                {
                    Graph::edge_descriptor ed;
                    bool temp_bool;
                    Graph::vertex_descriptor vd = edgeDest(pr, original_edge, graph);

                    tie(ed, temp_bool) = add_edge(progeny, vd, graph);
                    if (!temp_bool)

                    //Add to progeny iSite
                    newSite.edges.push_back(ed);
                    graph[progeny].edgeToSite[ed] = progeny_site;

                    //Update connected node's iSite
                    graph[vd].sites[graph[vd].edgeToSite[original_edge]].edges.push_back(ed);
                    graph[vd].edgeToSite[ed] = graph[vd].edgeToSite[original_edge];
                }

                //Remove edge from progenitor
                if (progenitorLoss[x])
                {
                    //Remove from iSite
                    graph[pr].edgeToSite.erase(original_edge);
                    graph[pr].sites[site].edges.erase(
                        graph[pr].sites[site].edges.begin()+edge);

                    //Update connected node
                    Graph::vertex_descriptor vd = edgeDest(pr, original_edge, graph);
                    int connectedSite = graph[vd].edgeToSite[original_edge];
                    int connectedSiteSize = graph[vd].sites[connectedSite].edges.size();
                    int k;
                    for (k=0; k<connectedSiteSize; ++k)
                        if (graph[vd].sites[connectedSite].edges[k]==original_edge) break;
                    graph[vd].edgeToSite.erase(original_edge);
                    graph[vd].sites[connectedSite].edges.erase(
                        graph[vd].sites[connectedSite].edges.begin()+k);

                    //Remove edge
                    remove_edge(original_edge, graph);
                    --edge;
                }

            } //Cycle through edges

            graph[progeny].sites.push_back(newSite);

        } //Cycle through iSites
    } //Cycle through progenitors

    //Check for empty progenitors
    for (int i=0; i<progenitors.size(); ++i)
        if (out_degree(progenitors[i], graph)==0)
            remove_vertex(progenitors[i], graph);
        else
        {
            bool self = true;
            AdjIter adji, adjiend;
            for (tie(adji, adjiend)=adjacent_vertices(progenitors[i], graph);
                 adji!=adjiend;++adji)
                if (*adji != progenitors[i])
                    self = false;
            if (self)
                remove_vertex(progenitors[i], graph);
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
        vertexsites& vref = graph[*vi];
        int numSites = vref.sites.size();
        for( int i=0; i < numSites; i++)
        {
            vref.sites[i].age += 1;
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
            vertexsites& vref = graph[*vi];
            int site = vref.edgeToSite[*oei];
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
                    site = vref.selfLoops[*oei].first;
                    int site2 = vref.selfLoops[*oei].second;
                    cout<<" "<<site<<"->"
                        <<site2<<":"
                        <<indexmap(vd);
                    printed_self_loops.push_back(*oei);
                }
            }
            else
            {
                int connectedSite = graph[vd].edgeToSite[*oei];
                cout<<" "<<site<<"->"
                    <<connectedSite<<":"
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
                vertexsites& v1ref = graph[vd1];
                cout<<"Node: " << indexmap(vd1) << endl;
                for (unsigned int i=0; i != v1ref.sites.size(); i++)
                    cout<<"\t"<<i<<":: Age: "<<v1ref.sites[i].age
                        <<", Edges: "<<v1ref.sites[i].edges.size()<<endl;
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
    if (argc!=11)
    {
        cerr<<"Usage: ./iSite <seed-graph> <probability of subfunctionalization> "
                "<probability of assymetry> <probability of homomeric subfunctionalization> "
                "<probability of fusion> <end order> <iterations> <output dir> <output file> "
                "<PBS job-id>"<<endl;
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
    int iterations=param.iterations;
    Random random(time(NULL)*getpid());
    //Random random(10000000, "randvals");
    rnd = &random;

    //Create output directory
    string outfile_path(argv[9]);
    if (mkdir(argv[8], 0777) != -1 || errno == EEXIST)
    {
        outfile_path.insert(0, "/");
        outfile_path.insert(0, argv[8]);
    }
    else
        perror("Unable to create output directory");

    //Opening output file
    ofstream outfile(outfile_path.c_str());
    if (!outfile)
    {
        cerr<<"Error opening output file: "<<outfile_path<<endl;
        exit(1);
    }
    outfile << "JOBID subfuncProb asymmetry selfloopLoss fusionProb actualAsymmetry selfloops order size tris trips CC numComponents lgComponentOrder lgComponentSize lgComponentTris lgComponentTips lgComponentCC" << endl;

    while (param.iterations--)
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
            vertexsites& v1ref = graph[vd1];
            vertexsites& v2ref = graph[vd2];
            Graph::edge_descriptor ed;


            bool temp_bool;
            tie(ed, temp_bool) = add_edge(vd1,vd2,graph);
            pair<int,int> edge_site_indices;

            //Add iSite
            map<string, int>::iterator known_site_it;
            known_site_it = v1ref.site_name_to_index.find(s1);
                //if site has not been added to this node before
            if (known_site_it == v1ref.site_name_to_index.end())
            {
                    //make new iSite with new edge in it
                v1ref.sites.push_back(isite(ed,0));
                    //set edgeToSite for new edge
                edge_site_indices.first = v1ref.sites.size()-1;
                v1ref.edgeToSite.insert(make_pair(ed, v1ref.sites.size()-1));
                    //set site_name_to_index for new iSite
                v1ref.site_name_to_index.insert(make_pair(s1, v1ref.sites.size()-1));
            }
                //if site is already on the node - just add edge to it
            else
            {
                    //insert new edge in iSite
                v1ref.sites[known_site_it->second].edges.push_back(ed);
                    //set edgeToSite for new edge
                edge_site_indices.first = known_site_it->second;
                v1ref.edgeToSite.insert(make_pair(ed, known_site_it->second));
            }

            if (!((v1==v2) && (s1==s2))) //check for self loops and same iSite
            {
                known_site_it = v2ref.site_name_to_index.find(s2);
                    //if site has not been added to this node before
                if (known_site_it == v2ref.site_name_to_index.end())
                {
                        //make new iSite with new edge in it
                    v2ref.sites.push_back(isite(ed,0));
                        //set edgeToSite for new edge
                    edge_site_indices.second = v2ref.sites.size()-1;
                    v2ref.edgeToSite.insert(make_pair(ed, v2ref.sites.size()-1));
                        //set site_name_to_index for new iSite
                    v2ref.site_name_to_index.insert(make_pair(s2, v2ref.sites.size()-1));
                }
                    //if site is already on the node - just add edge to it
                else
                {
                        //insert new edge in iSite
                    v2ref.sites[known_site_it->second].edges.push_back(ed);
                        //set edgeToSite for new edge
                    edge_site_indices.second = known_site_it->second;
                    v2ref.edgeToSite.insert(make_pair(ed, known_site_it->second));
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
                v1ref.edgeToSite[ed] = -1;
                v1ref.selfLoops.insert(make_pair(ed, edge_site_indices));
                //Mark site as homomeric
                v1ref.sites[edge_site_indices.first].homomeric = true;
            }
            
        }

        output_info(graph, indexmap, "***Original graph***", PRINT); 

#ifdef DEBUG
        cout << "Iteration: " << iterations-param.iterations << endl;
#endif
#ifdef NDEBUG
        //cout << "Working";
        cout << "Iteration: " << iterations-param.iterations << endl;
#endif
        while (num_vertices(graph)<param.end_order)
        {
            addAge(graph);
            duplication(graph, indexmap);

            output_info(graph, indexmap, "***Graph during algorithm***", PRINT); 
            //output_info(graph, indexmap, "", STATUS); 

            //Uncomment to print graph after each duplication
            //printGraph(graph, indexmap);
            //cout<<endl;
        }

        output_info(graph, indexmap, "***End Graph***", PRINT); 
        output_info(graph, indexmap, "***Node Summary***", NODE_SUMMARY); 

        //Uncomment to print resulting graph
        //printGraph(graph, indexmap);

#ifdef NDEBUG
        cout<<"Generating results"<<endl;
#endif
        //vector<int> dist;

        outfile<<argv[10]<<" ";
        outfile<<param.prob_loss<<" ";
        outfile<<param.prob_asym<<" ";
        outfile<<param.prob_self<<" ";
        outfile<<param.prob_fusion<<" ";
        outfile<<setprecision(3)<<asymmetry<<" ";
        outfile<<nSelfLoops<<" ";
        outfile<<num_vertices(graph)<<" ";
        /*
        cout<<"Before simplification"<<endl;
        isite_distribution(graph, dist);
        print_distribution(dist);
        cout<<endl;
        */
        simplify(graph);
        /*
        cout<<"After simplification"<<endl;
        dist.clear();
        isite_distribution(graph, dist);
        print_distribution(dist);
        */
        outfile<<num_edges(graph)<<" ";
        int numTriangles=triangles(graph);
        int numTriples=countTriples(graph);
        outfile<<numTriangles<<" ";
        outfile<<numTriples<<" ";
        outfile<<setprecision(6)<<(double)(3*numTriangles)/numTriples<<" ";
        outfile<<components(graph, indexmap)<<" ";
        //Largest component
        Graph lgcomponent = get_large_component(graph, indexmap);
        outfile<<num_vertices(lgcomponent)<<" ";
        outfile<<num_edges(lgcomponent)<<" ";
        numTriangles = triangles(lgcomponent);
        numTriples = countTriples(lgcomponent);
        outfile<<numTriangles<<" ";
        outfile<<numTriples<<" ";
        outfile<<setprecision(6)<<(double)(3*numTriangles)/numTriples;
        outfile<<endl;

        output_info(graph, indexmap, "***Node Evolution***", EVOLUTION); 

        infile.close();
    }
    outfile.close();
}


