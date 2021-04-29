#include "SubgraphStream.h"
#include <stdio.h>
#include <algorithm>


SubgraphStream::SubgraphStream(){
  streamType = "SubgraphStream";
}
	

int SubgraphStream::parseUpdate(vector<tuple<char,int,int>> * updateList, int i, int batchSize){
  int nError = 0;
  int nNodes = graph->numNodes();
  for(int j = i; j < i + batchSize; j++){
      char c = get<0>((*updateList)[j]);
      int a = get<1>((*updateList)[j]);
      int b = get<2>((*updateList)[j]);
      
      //printf("BASE - %d\n",base);
      if(a == b){
        printf("No self edges allowed\n");
      	nError++;
      	continue;
      }

    	a = a  - base;
      b = b  - base;
      if(a >= nNodes || b >= nNodes){
        printf("Edge exceedes graph size!\n");
        nError++;
        continue;
      }
      if(c != 'A' && c != 'R'){
        printf("Unknown command\n");
        nError++;
        continue;
      }
      if((c == 'A' && graph->hasEdge(a,b)) || (c == 'R' && !graph->hasEdge(a,b))){
        printf("\n\nInvalid edge operation\n\n");
        nError++;
        continue;
      }
      bool decrement = false;

      if(c == 'A'){

    		addEdge(a,b);
        fase->runCensusSubgraph(K,a,b,decrement);
      
    	}
    	else if(c == 'R'){

    		decrement = true;
        fase->runCensusSubgraph(K,a,b,decrement);
        removeEdge(a,b);
       
    	}
  }
  if(nError == batchSize) return 1;
  return 0;
}








