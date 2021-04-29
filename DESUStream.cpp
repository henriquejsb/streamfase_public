#include "DESUStream.h"
#include <stdio.h>
#include <algorithm>


DESUStream::DESUStream() {
  streamType = "DESUStream";
  

}



int DESUStream::parseUpdate(vector<tuple<char,int,int>> * updateList, int index, int batchSize){
  preprocess_timer->start();
  if(batchOperationsMatrix == NULL){
    //Inicializar a matrix batchOperationsMatrix NxN (N -> tamanho do grafo)
    //Cada elemento da matriz pode ser -2,-1,0,1 ou2
    // Se for 1 -> Essa aresta é uma aresta que irá ser adicionada (está em additions)
    // Se for 2 -> Essa aresta é uma aresta de additions que já foi expandida pelo DESU (aresta proibida)
    // Se for -1 -> Essa aresta é uma aresta que irá ser removida (está em deletions)
    // Se for -2 -> Essa aresta é uma aresta que já foi removida (aresta proibida)

    //Todos os valores são inicializados a 0
  	
  }

  //vector<pair<int,int>> additions;
  //vector<pair<int,int>> deletions;


  int nError = 0;
  int nAdditions = 0;
  int nDeletions = 0;
  int nNodes = graph->numNodes();
  for(int i = index; i < index + batchSize; i++){
      char c = get<0>((*updateList)[i]);
      int a = get<1>((*updateList)[i]);
      int b = get<2>((*updateList)[i]);
      //printf("BASE - %d\n",base);
      if(a == b){
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
      /*
      if((c == 'A' && graph->hasEdge(a,b)) || (c == 'R' && !graph->hasEdge(a,b))){
        printf("Invalid edge operation\n");
        nError++;
        continue;
      }
		*/


      if(c == 'A'){
        //Tratar as operações de adições 
      		if(batchOperationsMatrix[a][b] == -1){
      			//Quer dizer que neste batch já recebemos uma remoção da mesma aresta antes por isso ficam ambas sem efeito
      			//printf("FIX HERE\n");
            //printf("A after R\n");
            batchOperationsMatrix[a][b] = 0;
      			if(!directed) batchOperationsMatrix[b][a] = 0;

      			int j = 0;
      			while(j < nDeletions){
              //Removemos a operação de remoção que ficou sem efeito da lista de deletions
      				if(deletions[j].first == a && deletions[j].second == b){
      					//printf("FIX HERE\n");
                //deletions.erase(deletions.begin() + j);
                deletions[j].first = -1;
      					break;
      				}
      				if(!directed && deletions[j].first == b && deletions[j].second == a){
      					//printf("FIX HERE\n");
                deletions[j].first = -1;
                //deletions.erase(deletions.begin() + j);
      					break;
      				}
              j++;
      			}

            continue;
      		}else if(batchOperationsMatrix[a][b] == 1){
      			//Temos uma adição duplicada, não fazemos nada e não adicionados a operação ao batch
      			
            nError++;
      			continue;
      		}else if(graph->hasEdge(a,b)){
            //O grafo já tem a aresta que estamos a tentar adicionar, não fazemos nada
      			nError++;
        		continue;
      		}
          //Nenhuma das condições de erro se verificou, colocamos batchOperationsMatrix a 1 e adicionamos a operação
          //à lista additions
      		batchOperationsMatrix[a][b] = 1;
      		if(!directed) batchOperationsMatrix[b][a] = 1;
    		  //additions.push_back(make_pair(a,b));
          additions[nAdditions++] = make_pair(a,b);
          addEdge(a,b);
        	
    	}
    	else if(c == 'R'){
        //Análogo ao exemplo anterior mas para remoções
    		if(batchOperationsMatrix[a][b] == 1){
    			//Removing an edge that is going to be added before
            //printf("FIX HERE\n");
            //printf("R after A\n");
      			batchOperationsMatrix[a][b] = 0;
      			if(!directed) batchOperationsMatrix[b][a] = 0;
            removeEdge(a,b);
      			int j = 0;
      			while(j < nAdditions){
      				if(additions[j].first == a && additions[j].second == b){
      					//additions.erase(additions.begin() + j);
                additions[j].first = -1;
      					break;
      				}
      				if(!directed && additions[j].first == b && additions[j].second == a){
      					//additions.erase(additions.begin() + j);
                additions[j].first = -1;
      					break;
      				}
              j++;
      			}
            continue;
    		}else if(batchOperationsMatrix[a][b] == -1){
    			//Duplicated udpate
    			nError++;
    			continue;
    		}else if(!graph->hasEdge(a,b)){
    			nError++;
    			continue;
    		}
        //Nenhuma das condições de erro se verificou, colocamos batchOperationsMatrix a 1 e adicionamos a operação
          //à lista additions
          batchOperationsMatrix[a][b] = -1;
          if(!directed) batchOperationsMatrix[b][a] = -1;
          //deletions.push_back(make_pair(a,b));
          deletions[nDeletions++] = make_pair(a,b);

    		
        
       
    	}
    
  }

  //Já processámos as arestas do batch que estamos a tratar, agora vamos adicionar as novas arestas ao grafo
/*
  for(int i = 0; i < nAdditions; i++){
  	int a = additions[i].first;
  	int b = additions[i].second;
  	addEdge(a,b);
  }
*/

  //preprocess_timer->stop();
  //double elapsed = preprocess_timer->elapsed();

  //Passamos para a parte do DESU
  fase->runCensusBatch(K, additions, nAdditions, deletions, nDeletions, batchOperationsMatrix);

  //preprocess_timer->start();

  //printf("\t\t----------> %d - %d\n",additions.size(),deletions.size());

  //Damos reset à matrix batchOperationsMatrix para o próximo batch
  //e aplicamos as remoções de arestas
  for(int i = 0; i < nAdditions; i++){
    int a = additions[i].first;
    int b = additions[i].second;
    if(a == -1) continue;
    batchOperationsMatrix[a][b] = 0;
    if(!directed) batchOperationsMatrix[b][a] = 0;
  }

  for(int i = 0; i < nDeletions; i++){
    
  	int a = deletions[i].first;
  	int b = deletions[i].second;
    if(a == -1) continue;
    batchOperationsMatrix[a][b] = 0;
    if(!directed) batchOperationsMatrix[b][a] = 0;
  	removeEdge(a,b);
  }
  
  
  if(nError == batchSize) return 1;
  return 0;
}