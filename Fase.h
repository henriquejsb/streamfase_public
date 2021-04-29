#ifndef _FASE_
#define _FASE_

#include "Common.h"
#include "Graph.h"
#include "Random.h"
#include "Label.h"
#include "IGtrie.h"
#include "Isomorphism.h"

class Fase
{
 private:
  bool directed;
  bool sampling;
  Graph * graph;
  int K;
  int motifCount;
  IGtrie igtrie;
  map<string, int> canonicalTypes;
  int A,B;

  unsigned int Kmask;

  int* ALT_COMPS;
  int* OGN_COMPS;
  int ** ALT_Context;
  int ** OGN_Context;
  pair<int,int> * ALT_COMPS_Changes;
  pair<int,int> * OGN_COMPS_Changes;

  int* vsub;
  int* vextIndex;
  int* vextList;
  bool * usedNeighbour;
  int * usedNeighbourChanges;
  //int nUsedNeighbourChanges;
  int nExpanded;
  double* sampProb;
  char sadjM[MAXMOTIF * MAXMOTIF + 1];
  char nauty_s[MAXMOTIF * MAXMOTIF + 1];

  
  void updateExistingSubgraph(int currentVertex, bool decrement);

  void reduceCanonicalTypes();
  void expandEnumeration(int index, int depth, int labelNode, long long int label);
  void expandEnumerationStream(int index, int depth, int vextSize, int nUsedNeighbourChanges, int nOriginChanges, int labelNode, long long int label, bool decrement, bool existingSubgraph, int altLabelNode, long long int altlabel);
  void getSubgraphFrequency(pair<long long int, int> element, Isomorphism* iso);
  void checkConnectedEndpoint(int currentVertex, bool * connectedA, bool * connectedB);


  void printVSUB();
  //int ** OGN_COMPS;
  //int ** ALT_COMPS;
  //int nOGNComps;
  //int nALTComps;
  //bool * mergeOGNComp;
  //bool * mergeALTComp;

  bool * forbidden;
  int * forbiddenChanges;


  int ** batchOperationsMatrix;

 
  void DESU(int depth, int index, int vextSize, int nForbiddenChanges, int nUsedNeighbourChanges, int nOGNChanges, int nALTChanges, int labelNode, long long int label,int altLabelNode, long long int altlabel);
  void updateToExpand(int currentVertex, int depth, int * size , int * nForbiddenChanges, int * nUsedNeighbourChanges, int OGNmask, int * nOGNChanges, int ALTmask, int * nALTChanges);
  bool checkIn(vector<pair<int,int>> * operations, int v1, int v2);
  bool checkForbidden(int a, int b, int ** batchOperationsMatrix);
  bool checkOGNConnected(int a, int b, int ** batchOperationsMatrix);
  bool checkALTConnected(int a, int b, int ** batchOperationsMatrix);
  //void updateComponents(int depth, int v2, int ** batchOperationsMatrix);
 
 public:
  Fase(Graph* _g, int _K, bool _directed, bool debugprints);
  ~Fase();

  int getTypes();
  void runCensus();
  void initSampling(int sz, double* _sampProb);
  int getMotifCount() {return motifCount;}
  vector<pair<int, string> > subgraphCount();
  void runCensusSubgraph(int _K, int a, int b, bool decrement);

  void runCensusBatch(int K, pair<int,int> * additions, int nAdditions, pair<int,int> * deletions, int nDeletions, int ** batchOperationsMatrix);
};

#endif
