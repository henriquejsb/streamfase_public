#ifndef _LABEL_
#define _LABEL_

#include "Common.h"
#include "Graph.h"

class Label
{
 private:
  static Graph* graph;
  static bool directed;


 public:
  static void init(Graph* _g, bool _directed);
  static long long int updateLabel(int* vsub, int currentVertex, int depth);
  static long long int updateLabelDESUBefore(int* vsub, int currentVertex, int depth, int ** batchOperationsMatrix);
  static long long int updateLabelDESUAfter(int* vsub, int currentVertex, int depth, int ** batchOperationsMatrix);
  static int repDigits(int depth);
  static void fillNautyMatrix(char* sadjM, int K, long long int mask);
};

#endif
