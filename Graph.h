/* -------------------------------------------------
      _       _     ___                            
 __ _| |_ _ _(_)___/ __| __ __ _ _ _  _ _  ___ _ _ 
/ _` |  _| '_| / -_)__ \/ _/ _` | ' \| ' \/ -_) '_|
\__, |\__|_| |_\___|___/\__\__,_|_||_|_||_\___|_|  
|___/                                          
    
gtrieScanner: quick discovery of network motifs
Released under Artistic License 2.0
(see README and LICENSE)

Pedro Ribeiro - CRACS & INESC-TEC, DCC/FCUP

----------------------------------------------------
Partially Abstract Base Graph Class

Last Update: 11/02/2012
---------------------------------------------------- */

#ifndef _GRAPH_
#define _GRAPH_

#include "Common.h"

typedef enum{DIRECTED, UNDIRECTED} GraphType;

class Graph {
 public:
  virtual ~Graph() {};

  virtual void createGraph(int n, GraphType t) = 0; // create graph with n nodes
  virtual GraphType type() = 0;           // Graph Type
  virtual void zero() = 0;                // remove all connections
  bool isDirected(){ return !type();}
  virtual int addEdge(int a, int b) = 0; // add edge from a to b
  virtual int rmEdge(int a, int b)  = 0; // remove edge from a to b



  virtual int numNodes()  = 0; // Number of nodes in graph 
  virtual int numEdges()  = 0; // Number of edges in graph 
 
  virtual bool hasEdge(int a, int b) = 0; // Edge between a and b?
  virtual bool isConnected(int a, int b) = 0;  // Edge (a,b) or (b,a)?

  virtual int nodeOutEdges(int a) = 0;    // nr Edges from node a
  virtual int nodeInEdges(int a) = 0;     // nr Edges to   node a
  virtual int numNeighbours(int a) = 0;   // Nr Neighbours of node a

  virtual void sortNeighbours() = 0;       // All neighbours sorted in increasing order (sort vectors)
  virtual void sortNeighboursArray() = 0;  // All neighbours sorted in increasing order (sort arrays)

  virtual void makeArrayNeighbours() = 0;  // Create arrays of neighbours and discard vectors
  virtual void makeVectorNeighbours() = 0; // Create vectors of neighbours and discard arrays

  virtual void prepareGraph() = 0;

  virtual vector<int> *neighbours(int a) = 0; // Neighbours of node a
  virtual int **matrixNeighbours() = 0;            // Neighbours of node a in array form
  virtual int *arrayNeighbours(int a) = 0;         // Neighbours of node a in array form
  virtual int *arrayNumNeighbours() = 0;           // Numbers of neighbours in array form
  virtual vector<int> *outEdges(int a) = 0;   // Outgoing edges of node a
  virtual vector<int> *inEdges(int a) = 0;    // Ingoing edges of node a
//  virtual int getTag(int a) = 0;               // Tag of node a
  //STREAM
  virtual Graph * getEdgeNeighbourhood(int a, int b, int depth) = 0;
  //virtual void printNeighbourhood();
  virtual void printGraph() = 0;
  virtual void reverseEdge(int a, int b) = 0;
};

#endif
