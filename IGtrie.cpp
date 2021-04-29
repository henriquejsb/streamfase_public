#include "IGtrie.h"

IGtrie::IGtrie()
{
}

IGtrie::~IGtrie()
{
  destroy();
}

void IGtrie::init(int _K, bool dir)
{
  K = _K;
  directed = dir;
  maxLabels = 10 * _K * _K;
  numLabels = 1;
  labelPaths = (int**) malloc(sizeof(int*) * maxLabels);
  labelLeaf = (int*) malloc(sizeof(int) * maxLabels);
  labelCount = (int*) malloc(sizeof(int) * maxLabels);

  labelPaths[0] = new int[LB_WORD_SIZE];
  labelLeaf[0] = 1;
  labelCount[0] = 0;
  memset(labelPaths[0], -1, sizeof(int) * LB_WORD_SIZE);

  iso = new Isomorphism();
  iso->initNauty(K, directed);
}

void IGtrie::destroy()
{
  if (!numLabels)
    return;

  int i;
  for (i = 0; i < numLabels; i++)
    delete[] labelPaths[i];

  free(labelPaths);
  free(labelLeaf);
  free(labelCount);
  numLabels = 0;
  maxLabels = 0;
  iso->finishNauty();
  delete iso;
}

void IGtrie::expand()
{
  maxLabels *= 2;
  labelPaths = (int**) realloc(labelPaths, sizeof(int*) * maxLabels);
  labelLeaf = (int*) realloc(labelLeaf, sizeof(int) * maxLabels);
  labelCount = (int*) realloc(labelCount, sizeof(int) * maxLabels);  
}

string IGtrie::calculateIsomorphism(long long int label){
  Label::fillNautyMatrix(sadjM, K, label);
  nauty_s[0] = '\0';
  iso->canonicalStrNauty(sadjM, nauty_s);
  string str = string(nauty_s);
  return str;
}


int IGtrie::getTypes(){
  return (int)canonicalCounts.size();
}


void IGtrie::subgraphCount(vector<pair<int,string>> *subgraphVector){
  for (auto element : canonicalCounts)
    subgraphVector->push_back(make_pair(element.second, element.first));
}

void IGtrie::checkIsomorphism(int labelNode, long long int label){
  if(!canonicalTypes.count(labelNode)){
     Timer * iso_timer = new Timer();
     iso_timer->start();
     // printf("No canonical for -> %d\n",labelNode);
      string canon = calculateIsomorphism(label);
     // printf("%s\n",canon.c_str());
      canonicalTypes[labelNode] = canon;
     if(!canonicalCounts.count(canon)){
      canonicalCounts[canon] += 0;
     }
     iso_timer->stop();
     //printf("Iso Computation Time (ms): %0.6lf\n", iso_timer->elapsed());
     delete iso_timer;
    }
}


void IGtrie::incrementLabel(int labelNode, int value)
{
  //printf("INCREMENTING %d label by %d\n",labelNode,value);
  
  labelCount[labelNode] += value;
  canonicalCounts[canonicalTypes[labelNode]] += value;
  
}

int IGtrie::insertLabel(int labelNode, long long int label, int digits)
{
  if (labelPaths[labelNode][label & (LB_WORD_SIZE - 1)] == -1)
  {
    //If the label hasn't been inserted, creates new label node
    if (numLabels == maxLabels)
      expand();

    int newNode = numLabels++;
    //printf("newNode - %d\n",newNode);
    labelPaths[newNode] = new int[LB_WORD_SIZE];
    memset(labelPaths[newNode], -1, sizeof(int) * LB_WORD_SIZE);
    labelPaths[labelNode][label & (LB_WORD_SIZE - 1)] = newNode;
    labelLeaf[newNode] = ((digits <= LB_WORD_LEN) ? digits : 0);
    labelCount[newNode] = 0;
  }

  int nextNode = labelPaths[labelNode][label & (LB_WORD_SIZE - 1)];

  if (!labelLeaf[nextNode])
    return insertLabel(nextNode, label >> LB_WORD_LEN, digits - LB_WORD_LEN);
  return nextNode;
}

vector<pair<long long int, int> > IGtrie::enumerate(int K)
{
  enumeration.clear();
  enumerateFrom(0, 0, 0, 0, K - 1);

  return enumeration;
}

void IGtrie::enumerateFrom(int currentNode, long long int label, long long int parLabel, int parSize, int remaining)
{
  
  if (remaining == 0)
  {
    enumeration.push_back(make_pair(label, labelCount[currentNode]));
    return;
  }

  int i;
  for (i = 0; i < LB_WORD_SIZE; i++)
    if (labelPaths[currentNode][i] != -1)
    {
      long long int tmpLabel = parLabel;
      int tmpSize = parSize;

      if (labelLeaf[labelPaths[currentNode][i]])
      {
        int digits = labelLeaf[labelPaths[currentNode][i]];
        tmpLabel |= (i << tmpSize);
        tmpSize += digits;
        enumerateFrom(labelPaths[currentNode][i], ((label << tmpSize) | tmpLabel), 0, 0, remaining - 1);
      }
      else
      {
        int digits = LB_WORD_LEN;
        tmpLabel |= (i << tmpSize);
        tmpSize += digits;
        enumerateFrom(labelPaths[currentNode][i], label, tmpLabel, tmpSize, remaining);
      }
    }
}
