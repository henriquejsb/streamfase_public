#ifndef _IGTRIE_
#define _IGTRIE_

#include "Common.h"
#include "Isomorphism.h"
#include "Label.h"
#include "Timer.h"
#define LB_WORD_LEN 4
#define LB_WORD_SIZE 16 // 1 << LB_WORD_LEN

class IGtrie
{
 private:
  int maxLabels;
  int numLabels;
  int** labelPaths;
  int* labelLeaf;
  int* labelCount;
  int K;
  bool directed;
  vector<pair<long long int, int> > enumeration;
  map<string, int> canonicalCounts;
  map<int,string> canonicalTypes;
  char sadjM[MAXMOTIF * MAXMOTIF + 1];
  char nauty_s[MAXMOTIF * MAXMOTIF + 1];
  Isomorphism *iso;
  void expand();
  void enumerateFrom(int currentNode, long long int label, long long int parLabel, int parSize, int remaining);
  string calculateIsomorphism(long long int label);
 
 public:
  IGtrie();
  ~IGtrie();

  void init(int K, bool dir);
  void destroy();
  void checkIsomorphism(int labelNode, long long int label);
  void incrementLabel(int labelNode, int value);
  int insertLabel(int labelNode, long long int label, int digits);
  vector<pair<long long int, int> > enumerate(int K);
  void subgraphCount(vector<pair<int,string>> *subgraphVector);
  int getTypes();
};

#endif
