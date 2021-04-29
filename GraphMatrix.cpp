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
Graphs Implementation with Adj. Matrix and Adj. List

Last Update: 11/02/2012
---------------------------------------------------- */

#include "GraphMatrix.h"
#include "GraphUtils.h"
#include <stdio.h>
#include <algorithm>
#include <queue>

GraphMatrix::GraphMatrix() {
  _init();
}

GraphMatrix::~GraphMatrix() {
  _delete();
}

// ------------------------------
// Graph Creation
// ------------------------------

void GraphMatrix::_init() {
  _num_nodes = _num_edges = 0;

  _adjM             = NULL;
  _adjOut           = NULL;
  _adjIn            = NULL;
  _neighbours       = NULL;  
  _in               = NULL;
  _out              = NULL;
  _num_neighbours   = NULL;
  _array_neighbours = NULL;
  //_inNeighbourhood = NULL;
 
}

void GraphMatrix::_delete() {
  int i;

  if (_adjM!=NULL) {
    for (i=0; i<_num_nodes; i++)
      if (_adjM[i]!=NULL) delete[] _adjM[i];
    delete[] _adjM;
  }
  if (_adjIn!=NULL) delete[] _adjIn;
  if (_adjOut!=NULL) delete[] _adjOut;
  if (_neighbours!=NULL) delete[] _neighbours;

  if (_in!=NULL) delete[] _in;
  if (_out!=NULL) delete[] _out;
  if (_out!=NULL) delete[] _num_neighbours;

  if (_array_neighbours!=NULL) {
    for (i=0; i<_num_nodes; i++)
      if (_array_neighbours[i]!=NULL) free(_array_neighbours[i]);
    delete[] _array_neighbours;
  }
}

void GraphMatrix::zero() {
  int i,j;
  _num_edges = 0;

  for (i=0; i<_num_nodes;i++) {
    _in[i] = 0;
    _out[i] = 0;
    _num_neighbours[i] = 0;
    _adjIn[i].clear();
    _adjOut[i].clear();
    _neighbours[i].clear();
    for (j=0; j<_num_nodes;j++)
      _adjM[i][j]=false;
  }
}

void GraphMatrix::createGraph(int n, GraphType t) {
  int i;

  _num_nodes = n;
  _type = t;

  _delete();

  _adjM = new bool*[n];  
  for (i=0; i<n; i++) _adjM[i] = new bool[n];
  _adjIn      = new vector<int>[n];
  _adjOut     = new vector<int>[n];
  _neighbours = new vector<int>[n];

  _in             = new int[n]; 
  _out            = new int[n];
  _num_neighbours = new int[n];

  zero();
}

int GraphMatrix::addEdge(int a, int b) {
  //printf("ADDING %d and %d\n", a,b);
  if (_adjM[a][b]){
    //printf("Edge already exists between %d and %d\n",a,b);
    return 0;
  } 

  _adjM[a][b] = true;

  _adjOut[a].push_back(b);
  _out[a]++;

  _adjIn[b].push_back(a);
  _in[b]++;

  _num_edges++;

  if (!_adjM[b][a]) {
    _num_neighbours[a]++;
    _num_neighbours[b]++;
    if(_array_neighbours){
      _insertArrayNeighbours(a,b);
    }
    else{
      _neighbours[a].push_back(b);
      _neighbours[b].push_back(a);
      
    }
      
    
  }
  return 1;
}



int GraphMatrix::rmEdge(int a, int b) {
  //printf("Removing %d and %d\n",a,b);
  if (!_adjM[a][b]){
    printf("Edge does not exist\n");
    return 0;
  }

  _adjM[a][b] = false;

  _removeVector(_adjOut[a], b);
  _out[a]--;

  _removeVector(_adjIn[b], a);
  _in[b]--;

  _num_edges--;

  if (!_adjM[b][a]) {
    _num_neighbours[a]--;
    _num_neighbours[b]--;
    if(_array_neighbours){
      _removeArrayNeighbours(a,b);
    }
    else{

      _removeVector(_neighbours[a],b);
      _removeVector(_neighbours[b],a);
    }
    
  }
  return 1;
}

void GraphMatrix::_removeVector(vector<int> &v, int b) {
  int i, s = v.size();
  for (i=0; i<s; i++)
    if (v[i] == b) break;
  if (i<s) v.erase(v.begin()+i);
}

void GraphMatrix::sortNeighbours() {
  int i;
  for (i=0; i<_num_nodes; i++)
    sort(_neighbours[i].begin(), _neighbours[i].begin()+_neighbours[i].size());
}

void GraphMatrix::sortNeighboursArray() {
  int i;
  for (i=0; i<_num_nodes; i++)
    qsort(_array_neighbours[i], _num_neighbours[i], sizeof(int), GraphUtils::int_compare);
}

void GraphMatrix::makeArrayNeighbours() {
  int i,j;
  vector<int>:: iterator ii;
  _array_neighbours = new int*[_num_nodes];  
  for (i=0; i<_num_nodes; i++) {
    _array_neighbours[i] = (int *) malloc(sizeof(int) * _neighbours[i].size());
    //_array_neighbours[i] = new int[_neighbours[i].size()];
    for (ii=_neighbours[i].begin(), j=0; ii!=_neighbours[i].end(); ++ii, j++)
      _array_neighbours[i][j] = *ii;
    _neighbours[i].clear();
  }
}

void GraphMatrix::makeVectorNeighbours() {
  int i,j;
  vector<int>:: iterator ii;

  for (i=0; i<_num_nodes; i++)
    for (j=0; j<_num_neighbours[i]; j++)
      _neighbours[i].push_back(_array_neighbours[i][j]);

  if (_array_neighbours!=NULL) {
    for (i=0; i<_num_nodes; i++)
      if (_array_neighbours[i]!=NULL) delete[] _array_neighbours[i];
    delete[] _array_neighbours;
  }
}

void GraphMatrix::_printArrayNeighbours(int a){
  printf("Neighbours[%d] of %d: ",_num_neighbours[a],a);
  for(int i = 0; i < _num_neighbours[a]; i++){
    printf("%d , ",_array_neighbours[a][i]);
  }
  printf("\n");

}

void GraphMatrix::_insertArrayNeighbours(int a, int b){
  //printf("inserting %d and %d in ARRAY_NEIGHBOURS\n",a,b);
    _array_neighbours[a] = (int*) realloc(_array_neighbours[a],_num_neighbours[a] * sizeof(int));
    _array_neighbours[b] = (int*) realloc(_array_neighbours[b],_num_neighbours[b] * sizeof(int));
    _array_neighbours[a][_num_neighbours[a]-1] = b;
    _array_neighbours[b][_num_neighbours[b]-1] = a;
    return;

}


void GraphMatrix::_removeArrayNeighbours(int a, int b){
  int * _aux_neighbours_a = (int *) malloc(sizeof(int) * _num_neighbours[a]);
  int * _aux_neighbours_b = (int *) malloc(sizeof(int) * _num_neighbours[b]);
  int i = 0, j = 0;
  while(j < _num_neighbours[a]){
    if(_array_neighbours[a][i] == b){
      i++;
      continue;
    }
    else{
      _aux_neighbours_a[j] = _array_neighbours[a][i];
      i++;
      j++;
    }
  }
  i = 0;
  j = 0;
  while(j < _num_neighbours[b]){
    if(_array_neighbours[b][i] == a){
      i++;
      continue;
    }
    else{
      _aux_neighbours_b[j] = _array_neighbours[b][i];
      i++;\
      j++;
    }
  }
  free(_array_neighbours[a]);
  free(_array_neighbours[b]);
  _array_neighbours[a] = _aux_neighbours_a;
  _array_neighbours[b] = _aux_neighbours_b;
}





Graph * GraphMatrix::getEdgeNeighbourhood(int a, int b, int depth){
  //printf("\n*********************\nGETTING EDGE NEIGHBOURHOOD OF %d and %d\n*********************\n",a,b);
  //nodeQueue<int> nodeQueue;
  //printf("checkpoint1");
  Graph * newGraph = new GraphMatrix();
  
  int * correspIndex = new int[_num_nodes];
  memset(correspIndex,-1,sizeof(int) * _num_nodes);
  correspIndex[a] = 0;
  correspIndex[b] = 1;
  int currIndex = 2;
  vector<pair<int,int>> newEdges;
  //_neighbourhood.clear();
  //_exclusiveNeighbourhood.clear();
  queue<int> nodeQueue;
  queue<int> depthQueue;
  bool * isNeighbour = new bool[_num_nodes];
 // bool * isExclusive = new bool[_num_nodes];
  //_inNeighbourhood = new bool[_num_nodes];
  //memset(_inNeighbourhood,false,sizeof(bool) * _num_nodes);
  memset(isNeighbour,false,sizeof(bool) * _num_nodes);
  //memset(isExclusive,false,sizeof(bool) * _num_nodes);
  nodeQueue.push(a);
  //printf("pushing a into nodeQueue: %d\n",a);
  //printf("checkpoint1\n");
  depthQueue.push(depth-2);
  isNeighbour[b] = true;
  while(!nodeQueue.empty()){
    //printf("nodeQueue in 'a' not empty!\n");
    int currNode = nodeQueue.front();
    nodeQueue.pop();
    int currDepth = depthQueue.front();
    //printf("In currNode %d at depth %d\n",currNode,currDepth);
    depthQueue.pop();
    if(currNode == b) continue;
    if(correspIndex[currNode] == -1){
      correspIndex[currNode] = currIndex;
      currIndex++;
    }
    //printf("Currently at node %d at depth %d\n",currNode,currDepth);
    isNeighbour[currNode] = true;
    //isExclusive[currNode] = true;
   // _neighbourhood.push_back(currNode);
    //_inNeighbourhood[currNode] = true;
    if(currDepth == 0) continue;
    for(int i = 0; i < _num_neighbours[currNode]; i++){
      int currNeighbour = _array_neighbours[currNode][i];
      //printf("ANALYZING NEIGHBOURS OF %d - NUM: %d - CURRENT: %d\n",currNode,_num_neighbours[currNode],currNeighbour);
     // if(_adjM[currNode][currNeighbour] && !newGraph->hasEdge(currNode,currNeighbour)) newGraph->addEdge(currNode,currNeighbour);
     // if(_adjM[currNeighbour][currNode] && !newGraph->hasEdge(currNeighbour,currNode)) newGraph->addEdge(currNeighbour,currNode);
      if(_adjM[currNode][currNeighbour]) newEdges.push_back(make_pair(currNode,currNeighbour));
      if(_adjM[currNeighbour][currNode]) newEdges.push_back(make_pair(currNeighbour,currNode));
      if(isNeighbour[currNeighbour]) continue;
      nodeQueue.push(currNeighbour);
      depthQueue.push(currDepth-1);
    }
    
  }
  //printf("checkpoint2\n");
  memset(isNeighbour,false,sizeof(bool) * _num_nodes);
  isNeighbour[a] = true;
  nodeQueue.push(b);
  depthQueue.push(depth-2);
  //printf("Moving on to the second part!\n");
   while(!nodeQueue.empty()){
    //printf("nodeQueue in 'b' not empty!\n");
    int currNode = nodeQueue.front();

    nodeQueue.pop();
    int currDepth = depthQueue.front();
    //printf("In currNode %d at depth %d\n",currNode,currDepth);
    depthQueue.pop();
    if (currNode == a) continue;
    isNeighbour[currNode] = true;
    //printf("Currently at node %d at depth %d\n",currNode,currDepth);
    /*if(isExclusive[currNode]){
      //_exclusiveNeighbourhood.push_back(currNode);
      //printf("This node is exclusive!\n");
    } */
    if(correspIndex[currNode] == -1){
      correspIndex[currNode] = currIndex;
      currIndex++;
    }
    if(currDepth == 0) continue;

    for(int i = 0; i < _num_neighbours[currNode]; i++){
      int currNeighbour = _array_neighbours[currNode][i];
      //printf("ANALYZING NEIGHBOURS OF %d - NUM: %d - CURRENT: %d\n",currNode,_num_neighbours[currNode],currNeighbour);
      //if(_adjM[currNode][currNeighbour]) newGraph->addEdge(currNode,currNeighbour);
      //if(_adjM[currNeighbour][currNode]) newGraph->addEdge(currNeighbour,currNode);
     // if(_adjM[currNode][currNeighbour] && !newGraph->hasEdge(currNode,currNeighbour)) newGraph->addEdge(currNode,currNeighbour);
     // if(_adjM[currNeighbour][currNode] && !newGraph->hasEdge(currNeighbour,currNode)) newGraph->addEdge(currNeighbour,currNode);
      if(_adjM[currNode][currNeighbour]) newEdges.push_back(make_pair(currNode,currNeighbour));
      if(_adjM[currNeighbour][currNode]) newEdges.push_back(make_pair(currNeighbour,currNode));
      if(isNeighbour[currNeighbour]) continue;
      nodeQueue.push(currNeighbour);
      depthQueue.push(currDepth-1);
    }
    
  }

  //_exclusiveNeighbourhood.push_back(a);
  //_exclusiveNeighbourhood.push_back(b);
  if(_adjM[a][b] ){
    //printf("Adding %d and %d\n",a,b);
    newEdges.push_back(make_pair(a,b));
  }
  if(_adjM[b][a] ){
    //printf("Adding %d and %d\n",b,a);
    newEdges.push_back(make_pair(b,a));
  }


  newGraph->createGraph(currIndex,_type);
  for(int i = 0; i < newEdges.size(); i++){
    pair<int,int> currEdge = newEdges[i];
    int correspA , correspB;
    correspA = correspIndex[currEdge.first];
    correspB = correspIndex[currEdge.second];
    if(!newGraph->hasEdge(correspA,correspB)){
      newGraph->addEdge(correspA,correspB);
    }
  }

  
  delete [] isNeighbour;
  delete [] correspIndex;
  //delete [] isExclusive;
  //printf("\n*********************\nFINISHING EDGE NEIGHBOURHOOD OF %d and %d\n*********************\n",a,b);
  return newGraph;
}

/*
void GraphMatrix::printNeighbourhood(){
  printf("///////////////\n");
  printf("NEIGHBOURHOOD[%d]:\n",_neighbourhood.size());
  for(int i = 0; i < _neighbourhood.size(); i++){
    printf("%d ",_neighbourhood[i]);
  }
  printf("\n");
  printf("///////////////\n\n");
  printf("///////////////\n");
  printf("EXCLUSIVE NEIGHBOURHOOD[%d]:\n",_exclusiveNeighbourhood.size());
  for(int i = 0; i < _exclusiveNeighbourhood.size(); i++){
    printf("%d ",_exclusiveNeighbourhood[i]);
  }
  printf("\n");
  printf("///////////////\n");
}

void GraphMatrix::getSubgraphNeighbourhood(int a, int b, int depth){
  int i,j;
  vector<int>:: iterator ii;
  getEdgeNeighbourhood(a,b,depth);
  _array_subgraph_neighbours = new int*[_neighbourhood.size()];  
  for (i=0; i<_neighbourhood.size(); i++) {
    _array_subgraph_neighbours[i] = new int[_neighbours[i].size()];
    for (ii=_neighbours[i].begin(), j=0; ii!=_neighbours[i].end(); ++ii, j++)
      _array_neighbours[i][j] = *ii;
    //_neighbours[i].clear();
  }
}
*/


void GraphMatrix::printGraph(){
  //printf("hello!\n");
  for(int i = 0; i < _num_nodes; i++){
    if(_num_neighbours[i] > 0){
      printf("[%d] -> ",i);
      for(int j = 0; j < _num_neighbours[i]; j++){
        printf(" %d ,",_array_neighbours[i][j]);
      }
      printf("\n");
    }
  }
}

//REFAZER FUNÇAO MAS PARA CRIAR LOGO ARRAY_NEIGHBOURS A MEDIDA QUE SE FAZ BFS