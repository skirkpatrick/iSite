iSite
-------------
The [iSite algorithm][1] is a proposed model for the evolution of protein
interaction networks developed by Todd A. Gibson and Debra S. Goldberg.

This project is a C++ implementation of the iSite algorithm.

[1]: http://bioinformatics.oxfordjournals.org/content/27/3/376.full

Compiling
-------------
To compile the standard version with basic, end-state output, run "make".

To compile the developer version with debug statements, run "make debug".

Debug Statements
------------------
Debug statements are formatted as follows:

###Building graph from file
####Outputs the rules for the seed graph as it is read from the file.
<em>NodeName</em> refers to a string identifier for a node.
<em>SiteName</em> refers to a string identifier for an iSite.
    NodeName:SiteName<->SiteName:NodeName
    NodeName:SiteName<->SiteName:NodeName

###Graph output
####Outputs the current graph
<em>NodeIndex</em> refers to the mapped index of a node.
<em>SiteName</em> refers to a string identifier for an iSite of <em>NodeIndex</em>.
<em>ConnectedNode</em> refers to the mapped index of a node connected to <em>NodeIndex</em> at <em>SiteName</em>.
    NodeIndex: SiteName->ConnectedNode SiteName->ConnectedNode...
    NodeIndex: SiteName->ConnectedNode SiteName->ConnectedNode...

###Duplication
####Outputs data regarding node duplication
<em>ParentIndex</em> refers to the mapped index of the node chosen to be duplicated.
<em>ChildIndex</em> refers to the mapped index of the newly created node.
<em>NumberOfiSites</em> refers to the number of iSites that were duplicated.
The node listed under <em>Asymetry</em> is the node chosen by the probability of asymmetry.
The value under <em>Loss</em> indicates whether a particular edge is lost as determined by the probability of loss.
    Parent: ParentIndex
    Child: ChildIndex
    iSites: NumberOfiSites
    (blankline)
    Asymetry: \{ParentIndex | ChildIndex\}
        Loss: \{Yes | No\}
    (blankline)
    .
    .
    .

###Node summary
####Outputs a summary of each node and the corresponding iSites therein
<em>NodeIndex</em> refers to the mapped index of a node.
<em>SiteName</em> refers to a string identifier for an iSite of <em>NodeIndex</em>.
<em>SiteAge</em> refers to the age of <em>SiteName</em> in number of duplications that have passed since the iSite's creation.
<em>NumberOfEdges</em> refers to the number of edges connected to <em>SiteName</em>.
    Node: NodeIndex
        SiteName:: Age: SiteAge, Edges: NumberOfEdges
    (blankline)
    Node: NodeIndex
        SiteName:: Age: SiteAge, Edges: NumberOfEdges

###Node evolution
####Outputs a trace of the nodes from which each site was duplicated
<em>NodeIndex</em> refers to the mapped index of a node.
<em>OriginNode</em> refers to the original node that was duplicated.
<em>DuplicationN</em> refers to the Nth duplication, duplicated from the N-1th duplication.
<em>ParentNode</em> refers to the node from which <em>NodeIndex</em> was directly duplicated.
    NodeIndex: OriginNode->Duplication1->...->ParentNode->NodeIndex
    NodeIndex: OriginNode->Duplication1->...->ParentNode->NodeIndex
