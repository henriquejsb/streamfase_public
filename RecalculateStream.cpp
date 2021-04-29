#include "RecalculateStream.h"
#include <stdio.h>
#include <algorithm>


RecalculateStream::RecalculateStream() : Stream(){
  streamType = "RecalculateStream";

}


int RecalculateStream::parseUpdate(vector<tuple<char,int,int>> * updateList, int index, int batchSize){
  int nError = 0;
  int nNodes = graph->numNodes();
  for(int i = index; i < index + batchSize; i++){
    char c = get<0>((*updateList)[i]);
    int a = get<1>((*updateList)[i]);
    int b = get<2>((*updateList)[i]); 
    a = a - base ;
    b = b - base;

    if(a == b){
        printf("No self edges allowed\n");
        nError++;
        continue;
      }
    
     if(c != 'A' && c != 'R'){
      printf("Unknown command\n");
      nError += 1;
      continue;
    }
    
    if(a >= nNodes || b >= nNodes){
        printf("Edge exceedes graph size!\n");
        nError++;
        continue;
      }
    if((c == 'A' && graph->hasEdge(a,b)) || (c == 'R' && !graph->hasEdge(a,b))){
      printf("Invalid edge operation\n");
      nError += 1;
      continue;
    }


    if(c == 'A'){
      addEdge(a,b);
    }
    else if(c == 'R'){
      removeEdge(a,b);
    }
  }
  singleCensus();
  //delete fase;
  if(nError == batchSize) return 1;
  return 0;

}



