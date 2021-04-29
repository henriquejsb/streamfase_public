#include "Fase.h"

Fase::Fase(Graph* _g, int _K, bool _directed,bool _debugprints)
{
  directed = _directed;
  graph = _g;
  K = _K;
  Kmask = (1U << K) - 1U;
  sampling = false;
  vextList = new int[graph->numNodes()];

  vextIndex = new int[graph->numNodes()];
  memset(vextIndex,0,sizeof(int)* graph->numNodes());


  ALT_COMPS = new int[graph->numNodes()];
  OGN_COMPS = new int[graph->numNodes()];
  memset(ALT_COMPS,0,sizeof(int) * graph->numNodes());
  memset(OGN_COMPS,0,sizeof(int) * graph->numNodes());
  OGN_COMPS_Changes = new pair<int,int>[graph->numNodes() * K];
  ALT_COMPS_Changes = new pair<int,int>[graph->numNodes() * K];

  OGN_Context = new int*[K];
  ALT_Context = new int*[K];
  for(int i = 0; i < K; i++){
    OGN_Context[i] = new int[i+2];
    ALT_Context[i] = new int[i+2];
    memset(OGN_Context[i],0,sizeof(int) * (i+2));
    memset(ALT_Context[i],0,sizeof(int) * (i+2));
    //OGN_Context[i][0] = 0;
    //ALT_Context[i][0] = 0;
  }


  vsub = new int[MAXMOTIF];
  sampProb = new double[MAXMOTIF];
  A = -1;
  B = -1;
  Label::init(_g, _directed);


  usedNeighbour = new bool[graph->numNodes()];
  memset(usedNeighbour,false,sizeof(bool) * graph->numNodes());
  usedNeighbourChanges = new int[graph->numNodes()];
  memset(usedNeighbourChanges,0,sizeof(int) * graph->numNodes());
  

  forbidden = new bool[graph->numNodes()];
  memset(forbidden,false,sizeof(bool) * graph->numNodes());
  forbiddenChanges = new int[graph->numNodes()];
  memset(forbiddenChanges,0,sizeof(int) * graph->numNodes());
  


  
}

Fase::~Fase()
{

  //delete[] vextOrigin;
  //delete[] originChanges;
  delete[] vsub;
  delete[] sampProb;
  delete[] vextList;
  delete[] usedNeighbour;

  delete[] vextIndex;
 
  delete[] forbidden;
  
  /*delete[] mergeOGNComp;
  delete[] mergeALTComp;
  for(int i = 1; i < K; i++){
    delete[] OGN_COMPS[i];
    delete[] ALT_COMPS[i];
  }*/
  delete[] OGN_COMPS;
  delete[] ALT_COMPS;
  delete[] OGN_COMPS_Changes;
  delete[] ALT_COMPS_Changes;

  for(int i = 0; i < K; i++){
    delete[] ALT_Context[i];
    delete[] OGN_Context[i];
   }
   delete[] ALT_Context;
   delete[] OGN_Context;


  delete[] forbiddenChanges;
  delete[] usedNeighbourChanges;


  igtrie.destroy();
}

void Fase::initSampling(int sz, double* _sampProb)
{
  int i;
  for (i = 0; i < sz; i++)
    sampProb[i] = _sampProb[i];

  sampling = true;
}

/*-------------------------------------------------------------------------------------------
|
|
|
|                       FASE normal para o censo inicial
|
|
---------------------------------------------------------------------------------------------*/




void Fase::runCensus()
{
  //K = _K;
  motifCount = 0;
  //printf("Init\n");
  igtrie.init(K,directed);
  for (int i = 0; i < graph->numNodes() + 1 - K; i++){
    vsub[0] = i;
    int *nei = graph->arrayNeighbours(i);
    int neiNum = graph->numNeighbours(i);
    nExpanded = 0;
    //vextSz[1] = 0;
    for (int j = 0; j < neiNum; j++)
      if (nei[j] > i){
        vextList[nExpanded++] = nei[j];
        //vext[1][vextSz[1]++] = nei[j];
      }
  
    expandEnumeration(0, 1, 0, 0LL);
  }
}



void Fase::expandEnumeration(int index, int depth, int labelNode, long long int label)
{
  if (depth == K - 1)
  {
    //while (vextSz[depth])
    for(int I = index; I < nExpanded; I++)
    {
      //int currentVertex = vext[depth][--vextSz[depth]];
      int currentVertex = vextList[I];

      //if (!sampling || Random::testProb(sampProb[depth]))
      //{
     
      long long int clabel = Label::updateLabel(vsub, currentVertex, depth);
      long long int auxlabel = label;
      auxlabel = (label << Label::repDigits(depth)) | clabel;
      int gtrieNode = igtrie.insertLabel(labelNode, clabel, Label::repDigits(depth));
      igtrie.checkIsomorphism(gtrieNode,auxlabel);
      igtrie.incrementLabel(gtrieNode, 1);
      motifCount++;
      
      //}
    }

    return;
  }

  int i, j;
  long long int clabel;
  int clabelNode = labelNode;
  


  int auxNExpanded = nExpanded;
  for(int I = index; I < auxNExpanded; I++)
  {
  
    int currentVertex = vextList[I];
    long long int auxlabel = label;

   
    vsub[depth] = currentVertex;

    int *eExcl = graph->arrayNeighbours(currentVertex);
    int eExclNum = graph->numNeighbours(currentVertex);
    
    for (i = 0; i < eExclNum; i++)
    {
      if (eExcl[i] <= vsub[0])
        continue;

      for (j = 0; j < depth; j++)
        if (eExcl[i] == vsub[j] || graph->isConnected(eExcl[i], vsub[j]))
          break;

      if (j == depth)
        
        vextList[nExpanded++] = eExcl[i];
    }

    if (depth >= 1)
    {
      clabel = Label::updateLabel(vsub, currentVertex, depth);
      clabelNode = igtrie.insertLabel(labelNode, clabel, Label::repDigits(depth));
      auxlabel = (label << Label::repDigits(depth)) | clabel;
    }

    expandEnumeration(I+1, depth + 1, clabelNode, auxlabel);
  
    nExpanded = auxNExpanded;
  }
}

/*-------------------------------------------------------------------------------------------
|
|
|
|                       Operações de Stream (versão estável)
|
|
---------------------------------------------------------------------------------------------*/



void Fase::runCensusSubgraph(int _K, int a, int b, bool decrement){

  A = a;
  B = b;
  
  vsub[0] = A;
  vsub[1] = B;
  int *neiA = graph->arrayNeighbours(a);
  int neiNumA = graph->numNeighbours(a);
  int *neiB = graph->arrayNeighbours(b);
  int neiNumB = graph->numNeighbours(b);
  

 
  //memset(usedNeighbour,false,sizeof(bool) * graph->numNodes());
  
  //memset(OGN_COMPS,-1,sizeof(int) * graph->numNodes());

  long long int clabel = Label::updateLabel(vsub, B, 1);
  int clabelNode = igtrie.insertLabel(0, clabel, Label::repDigits(1));
 

  graph->reverseEdge(A,B);
  if(!directed) graph->reverseEdge(B,A);

  long long int altclabel = Label::updateLabel(vsub,B,1);
  int altclabelNode = igtrie.insertLabel(0,altclabel,Label::repDigits(1));

  graph->reverseEdge(A,B);
  if(!directed) graph->reverseEdge(B,A);


  bool existingSubgraph = false;
  

  if(directed && graph->hasEdge(B,A)){
    //printf("EVERYTHING GOES MY WAY!\n");
    existingSubgraph = true;
  }


  //nExpanded = 0;

  int vextSize = 0;
  int nUsedNeighbourChanges = 0;

  for (int j = 0; j < neiNumA; j++){
      if(neiA[j] == B) continue;
      OGN_COMPS[neiA[j]] = A;
      vextIndex[neiA[j]] = vextSize;
      vextList[vextSize++] = neiA[j];
      usedNeighbour[neiA[j]] = true;
      usedNeighbourChanges[nUsedNeighbourChanges++] = neiA[j];

   
  }
  
  for (int j = 0; j < neiNumB; j++){
      if(neiB[j] == A) continue;
      if(usedNeighbour[neiB[j]]){
        OGN_COMPS[neiB[j]] = A + B + 1;
        continue;
      } 
      vextIndex[neiB[j]] = vextSize;
      vextList[vextSize++] = neiB[j];
      usedNeighbour[neiB[j]] = true;
      usedNeighbourChanges[nUsedNeighbourChanges++] = neiB[j];
      OGN_COMPS[neiB[j]] = B;
      
  }

  forbidden[A] = true;
  forbidden[B] = true;
  
  
  expandEnumerationStream(0,2,vextSize,nUsedNeighbourChanges,0,clabelNode, clabel,decrement, existingSubgraph, altclabelNode, altclabel);
  
  forbidden[A] = false;
  forbidden[B] = false;

  for(int j = 0; j < nUsedNeighbourChanges; j++){
    int vertex = usedNeighbourChanges[j];
    usedNeighbour[vertex] = false;
  }


}

void Fase::printVSUB(){
  printf("VSUB: ");
  for(int i = 0; i < K; i++){
    printf("%d ",vsub[i]);
  }
  printf("\n");
  for(int i = 0; i < K; i++){
    int v1 = vsub[i];
    for(int j = 0; j < i; j++){
      int v2 = vsub[j];
      if(graph->hasEdge(v1,v2)){
        printf("%d->%d\n",v1,v2);
      }
    }
  }
}

void Fase::expandEnumerationStream(int index, int depth, int vextSize, int nUsedNeighbourChanges, int nOriginChanges, int labelNode, long long int label, bool decrement, bool existingSubgraph, int altLabelNode, long long int altlabel)
{

  if (depth == K - 1)
  {
   
    for(int I = index; I < vextSize; I++)
    {
     
      int currentVertex = vextList[I];
      bool auxExistingSubgraph = existingSubgraph || (OGN_COMPS[currentVertex] == A + B + 1);
      
      long long int auxlabel = label;
      long long int auxaltlabel = altlabel;


      long long int clabel = Label::updateLabel(vsub, currentVertex, depth);
      //printf("LABEL -> %d | LABELNODE -> %d\n",clabel,labelNode);
      
      auxlabel = (auxlabel << Label::repDigits(depth)) | clabel;
      int gtrieNode = igtrie.insertLabel(labelNode, clabel, Label::repDigits(depth));
      igtrie.checkIsomorphism(gtrieNode,auxlabel);
      igtrie.incrementLabel(gtrieNode, 1 -2 * decrement);

      motifCount += 1 -2*decrement;

      if(auxExistingSubgraph){
        //vsub[depth] = currentVertex;
        //printVSUB();
        motifCount += 1 - 2 * !decrement;
        auxaltlabel = (auxaltlabel << Label::repDigits(depth)) | clabel;
       
        int gtrieNode = igtrie.insertLabel(altLabelNode,clabel,Label::repDigits(depth));
        igtrie.checkIsomorphism(gtrieNode,auxaltlabel);
        igtrie.incrementLabel(gtrieNode, 1-2*!decrement);
      }
      
    }
    return;
  }

  int i, j;
  long long int clabel;
  int clabelNode = labelNode;

  
  int altclabelNode = altLabelNode;

 
  for(int I = index; I < vextSize; I++)
  {
    
    int auxVextSize = vextSize;
    int auxNUsedNeighbourChanges = nUsedNeighbourChanges;
    int auxNOriginChanges = nOriginChanges;
    int currentVertex = vextList[I];
    int currentOrigin = OGN_COMPS[currentVertex];

    bool auxExistingSubgraph = existingSubgraph || (currentOrigin == A + B + 1);

    long long int auxlabel = label;
    long long int auxaltlabel = altlabel;

   
      
   
    vsub[depth] = currentVertex;
    forbidden[currentVertex] = true;
   
    
    int *eExcl = graph->arrayNeighbours(currentVertex);
    int eExclNum = graph->numNeighbours(currentVertex);
    
    //vector<pair<int,int>> changesOriginVext;

    for (i = 0; i < eExclNum; i++)
    {
    
      int currNeighbour = eExcl[i];

      if(forbidden[currNeighbour]){
        //Neighbour already in vsub
        continue;
      }
      if(usedNeighbour[currNeighbour]){
        //Neighbour already expanded in vext
        //If this neighbour is connected to a vertex in the subgraph we are considering,
        //We check if they have different origins. This neighbour will already be in vExt,
        //But we notify future recursion calls that when they add this neighbour to the subgraph,
        //the subgraph considered will be an already existing subgraph
        if(vextIndex[currNeighbour] > I && OGN_COMPS[currentVertex] != A+B+1 && OGN_COMPS[currNeighbour] != A+B+1 && OGN_COMPS[currentVertex] != OGN_COMPS[currNeighbour]){
            //changesOriginVext.push_back(make_pair(currNeighbour,OGN_COMPS[currNeighbour]));
            //OGN_COMPS_Changes[auxNOriginChanges++] = make_pair(currNeighbour,OGN_COMPS[currNeighbour]);
            OGN_COMPS_Changes[auxNOriginChanges].first = currNeighbour;
            OGN_COMPS_Changes[auxNOriginChanges++].second = OGN_COMPS[currNeighbour];
            OGN_COMPS[currNeighbour] = A + B + 1;
          }
        continue;
      }
      //If none of the above happens, we add the neighbour to vExt and we 
      //define its origin as the same origin as currentVertex. We push this change
      //to undo after the recursion
      //vext[depth + 1][vextSz[depth + 1]++] = eExcl[i];
      vextIndex[currNeighbour] = auxVextSize;
      vextList[auxVextSize++] = currNeighbour;
      //changesOriginVext.push_back(make_pair(eExcl[i],OGN_COMPS[eExcl[i]]));
      OGN_COMPS[currNeighbour] = OGN_COMPS[currentVertex];
      usedNeighbour[currNeighbour] = true;
      usedNeighbourChanges[auxNUsedNeighbourChanges++] = currNeighbour;



  
      }

  
      clabel = Label::updateLabel(vsub, currentVertex, depth);
      clabelNode = igtrie.insertLabel(labelNode, clabel, Label::repDigits(depth));

      auxlabel = (auxlabel << Label::repDigits(depth)) | clabel;
      auxaltlabel = (auxaltlabel << Label::repDigits(depth)) | clabel;
      
      altclabelNode = igtrie.insertLabel(altLabelNode,clabel,Label::repDigits(depth));
     

      expandEnumerationStream(I+1, depth + 1, auxVextSize, auxNUsedNeighbourChanges, auxNOriginChanges, clabelNode, auxlabel, decrement, auxExistingSubgraph, altclabelNode, auxaltlabel);

      for(int i = nOriginChanges; i < auxNOriginChanges; i++){
        pair <int,int> change = OGN_COMPS_Changes[i];
        int vertex = change.first;
        int origin = change.second;
        OGN_COMPS[vertex] = origin;
      }
      for(int i = nUsedNeighbourChanges; i < auxNUsedNeighbourChanges; i++){
        int vertex = usedNeighbourChanges[i];
        usedNeighbour[vertex] = false;
      }

      forbidden[currentVertex] = false;


    

    
  }
}



/*-------------------------------------------------------------------------------------------
|
|
|
|                       DESU e operações em batch
|
|
---------------------------------------------------------------------------------------------*/




void Fase::runCensusBatch(int K, pair<int,int> * additions, int nAdditions, pair<int,int> * deletions, int nDeletions, int ** auxBatchOperationsMatrix){
  //Função que irá chamar o DESU para cada aresta nova ou removida
  batchOperationsMatrix = auxBatchOperationsMatrix;
  for(int i = 0; i < nAdditions + nDeletions; i++){

    //printf("In operation %d / %d\n",i+1,nAdditions + nDeletions);


    int A,B;
    int j = i - nAdditions;
    
    //Itera primeiro sobre as adições e depois sobre as remoções
    if(j < 0){
      A = additions[i].first;
      B = additions[i].second;
    }
    else{
      A = deletions[j].first;
      B = deletions[j].second;
    }
    if(A == -1) continue;
    if(checkForbidden(A,B,batchOperationsMatrix)){
      continue;
    }
    //Colocamos A e B em vsub e expandimos os seus vizinhos
    vsub[0] = A;
    vsub[1] = B;
    
    int * neiA = graph->arrayNeighbours(A);
    int neiNumA = graph->numNeighbours(A);
    int * neiB = graph->arrayNeighbours(B);
    int neiNumB = graph->numNeighbours(B);

    //TODO replace memsets with undo operations of array *Changes
    //memset(usedNeighbour,false,sizeof(bool) * graph->numNodes());
    //memset(forbidden,false,sizeof(bool) * graph->numNodes());
    
    //memset(ALT_COMPS,0,sizeof(int) * graph->numNodes());

    int nUsedNeighbourChanges = 0;
    int nForbiddenChanges = 0;
    int nALTChanges = 0;
    int nOGNChanges = 0;

    int ALTmaskA = 1;
    int ALTmaskB = 2;

    int OGNmaskA = 1;
    int OGNmaskB = 2;

    
    usedNeighbour[A] = true;
    usedNeighbour[B] = true;

    forbidden[A]  = true;
    forbidden[B]  = true;

    //memset(OGN_COMPS,0,sizeof(int) * graph->numNodes());

    int sizeVext = 0;

    if(j < 0){
    //Estamos a tratar uma adição
    //OGN_COMPS[1][0] = 0;
      if(directed && graph->hasEdge(B,A) && batchOperationsMatrix[B][A] <= 0){
        
        OGNmaskA |= OGNmaskB;
        OGNmaskB |= OGNmaskA;
        //OGN_COMPS[1][1] = 0;
        //nOGNComps = 1;
      }/*else{
        
        OGN_COMPS[1][1] = 1;
        nOGNComps = 2;
      }*/
      //ALT_COMPS[1][0] = 0;
      //ALT_COMPS[1][1] = 0;
      ALTmaskA |= ALTmaskB;
      ALTmaskB |= ALTmaskA;
    //nALTComps = 1;
   
    }
    else{
      //Estamos a tratar uma remoção
      //ALT_COMPS[1][0] = 0;
      if(directed && graph->hasEdge(B,A) && batchOperationsMatrix[B][A] >= 0){
        ALTmaskA |= ALTmaskB;
        ALTmaskB |= ALTmaskA;
        //ALT_COMPS[1][1] = 0;
        //nALTComps = 1;
      }/*else{
        ALT_COMPS[1][1] = 1;
        nALTComps = 2;
        
      }*/
      //OGN_COMPS[1][0] = 0;
      //OGN_COMPS[1][1] = 0;
      //nOGNComps = 1;
      OGNmaskA |= OGNmaskB;
      OGNmaskB |= OGNmaskA;
    
    }

    for (int j = 0; j < neiNumA; j++){
      if(neiA[j] == B) continue;
      if(checkForbidden(neiA[j],A,batchOperationsMatrix) || checkForbidden(neiA[j],B,batchOperationsMatrix)){
        forbiddenChanges[nForbiddenChanges++] = neiA[j];
        forbidden[neiA[j]] = true;
        continue;
      }
      
      vextIndex[neiA[j]] = sizeVext;
      vextList[sizeVext++] = neiA[j];
      usedNeighbour[neiA[j]] = true;
      usedNeighbourChanges[nUsedNeighbourChanges++] = neiA[j];
      if(checkALTConnected(A,neiA[j],batchOperationsMatrix)){
        ALT_COMPS[neiA[j]] = ALTmaskA;
      }else{
        ALT_COMPS[neiA[j]] = 0;
      }
      if(checkOGNConnected(A,neiA[j],batchOperationsMatrix)){
        OGN_COMPS[neiA[j]] = OGNmaskA;
      }else{
        OGN_COMPS[neiA[j]] = 0;
      }
   
    }

    for (int j = 0; j < neiNumB; j++){
      if(neiB[j] == A) continue;
      if(forbidden[neiB[j]]){
        continue;
      }
      if(usedNeighbour[neiB[j]]){
        if(checkALTConnected(B,neiB[j],batchOperationsMatrix)){
          ALT_COMPS[neiB[j]] |= ALTmaskB;
        }
        if(checkOGNConnected(B,neiB[j],batchOperationsMatrix)){
          OGN_COMPS[neiB[j]] |= OGNmaskB;
        }
        continue;
      } 
      if(checkForbidden(neiB[j],A,batchOperationsMatrix) || checkForbidden(neiB[j],B,batchOperationsMatrix)){
        forbiddenChanges[nForbiddenChanges++] = neiB[j];
        //TODO is this really necessary? cant we just avoid expansion?
        forbidden[neiB[j]] = true;
        continue;
      }
      vextIndex[neiB[j]] = sizeVext;
      vextList[sizeVext++] = neiB[j];
      usedNeighbour[neiB[j]] = true;
      usedNeighbourChanges[nUsedNeighbourChanges++] = neiB[j];
      if(checkALTConnected(B,neiB[j],batchOperationsMatrix)){
        ALT_COMPS[neiB[j]] = ALTmaskB;
      }else{
        ALT_COMPS[neiB[j]] = 0;
      }
      if(checkOGNConnected(B,neiB[j],batchOperationsMatrix)){
        OGN_COMPS[neiB[j]] = OGNmaskB;
      }else{
        OGN_COMPS[neiB[j]] = 0;
      }
   
    }

  

    //Para começar a construir as labels e o caminho da g-trie para subgrafos novos

    long long int clabel = Label::updateLabelDESUAfter(vsub, B, 1,batchOperationsMatrix);
    int clabelNode = igtrie.insertLabel(0, clabel, Label::repDigits(1));
 

    //Para começar a construir as labels e o caminho da g-trie para subgrafos que vão deixar de existir

    long long int altclabel = Label::updateLabelDESUBefore(vsub,B,1,batchOperationsMatrix);
    int altclabelNode = igtrie.insertLabel(0,altclabel,Label::repDigits(1));
    

  
  

    

    //printf("Starting DESU\n");
    //Chamada ao DESU para a aresta que estamos a tratar
    

    DESU(2,0,sizeVext,nForbiddenChanges, nUsedNeighbourChanges, 0, 0, clabelNode,clabel,altclabelNode,altclabel);

    
        
   
    //Multiplicamos por 2 para indicar que esta aresta agora é proibida
    //Se for adição passa de 1 a 2, se for remoção passa de -1 a -2
    batchOperationsMatrix[A][B] *= 2;
    if(!directed) batchOperationsMatrix[B][A] *= 2;

    for(int j = 0; j < nForbiddenChanges; j++){
      int vertex = forbiddenChanges[j];
      forbidden[vertex] = false;
    }
    for(int j = 0; j < nUsedNeighbourChanges; j++){
      int vertex = usedNeighbourChanges[j];
      usedNeighbour[vertex] = false;
    }
    forbidden[A] = false;
    forbidden[B] = false;
    usedNeighbour[A] = false;
    usedNeighbour[B] = false;
  }
  

  


}

void Fase::DESU(int depth, int index, int sizeVext, int nForbiddenChanges, int nUsedNeighbourChanges, int nOGNChanges, int nALTChanges, int labelNode, long long int label,int altLabelNode, long long int altlabel){

  
  if(depth == K){

    //Chegámos ao fim da recursão

    //printf("\tFINISHED A RECUSION!\n");
    //int res = 0;
    //if(currALTmask == Kmask){
    if(ALT_Context[depth-1][0] == 1 && ALT_Context[depth-1][1] == Kmask){
    //if(nALTComps == 1){
      //ALT OCCURENCE
    //  res++;
      igtrie.checkIsomorphism(labelNode,label);
      igtrie.incrementLabel(labelNode,1);
      motifCount++;
      


    }
    //if(nOGNComps == 1){
    //if(currOGNmask == Kmask){
    if(OGN_Context[depth-1][0] == 1 && OGN_Context[depth-1][1] == Kmask){
      //OGN OCCURENCE
     // res++;
      igtrie.checkIsomorphism(altLabelNode,altlabel);
      igtrie.incrementLabel(altLabelNode,-1);
      motifCount--;
      //printf("ALTETED OCCURENCE!\n");
    }
    /*if(res == 2){
      printVSUB();
    }*/


 
    return;
  }


  for(int I = index; I < sizeVext; I++){
 
      int currentVertex = vextList[I];

      if(forbidden[currentVertex]){
        continue;
      }


      int auxNExpanded = sizeVext;
      int auxNUsedNeighbourChanges = nUsedNeighbourChanges;
      int auxNForbiddenChanges = nForbiddenChanges;

      int auxALTmask;
      int vertexALTmask = ALT_COMPS[currentVertex];
      int auxNALTChanges = nALTChanges;

      int auxOGNmask;
      int vertexOGNmask = OGN_COMPS[currentVertex];
      int auxNOGNChanges = nOGNChanges;


      vertexALTmask |= (1U << depth);
      vertexOGNmask |= (1U << depth);


      /*
      if(auxALTmask & vertexALTmask) auxALTmask |= vertexALTmask;
      else auxALTmask = vertexALTmask;

      if(auxOGNmask & vertexOGNmask) auxOGNmask |= vertexOGNmask;
      else auxOGNmask = vertexOGNmask;

      //auxALTmask |= (1U << depth);

      auxOGNmask |= (1U << depth);

    */
      vsub[depth] = currentVertex;



      OGN_Context[depth][0] = 0;
      ALT_Context[depth][0] = 0;

      int auxKmask = (1U << (depth+1)) - 1;

      for(int i = 1; i < OGN_Context[depth-1][0] + 1; i++){
        //printf("hey\n");
        auxOGNmask = OGN_Context[depth-1][i];
        if(auxOGNmask & vertexOGNmask){

          vertexOGNmask |= auxOGNmask;
          if(vertexOGNmask == auxKmask) break;
        }
        else{
          OGN_Context[depth][1 + OGN_Context[depth][0]++] = auxOGNmask;
        }

      }
      OGN_Context[depth][1 + OGN_Context[depth][0]++] = vertexOGNmask;

      //printf("OGN CONTEXT DEPTH 2 SIZE = %d\n",OGN_Context[depth][0]);

      for(int i = 1; i < ALT_Context[depth-1][0] + 1; i++){
        auxALTmask = ALT_Context[depth-1][i];
        if(auxALTmask & vertexALTmask){

          vertexALTmask |= auxALTmask;
          if(vertexALTmask == auxKmask) break;
        }
        else{
          ALT_Context[depth][1 + ALT_Context[depth][0]++] = auxALTmask;
        }

      }
      ALT_Context[depth][1 + ALT_Context[depth][0]++] = vertexALTmask;

      //printf("ALT CONTEXT DEPTH 2 SIZE = %d\n",ALT_Context[depth][0]);

      
      
      //Setting current vertex to forbidden to distinguish between already expanded vertices
      //and already considered vertices
      forbidden[currentVertex] = true;

      //Tratamento das labels e caminhos da G-Trie
      long long int auxlabel = label;
      long long int auxaltlabel = altlabel;
      long long int clabel;
      int clabelNode;
      long long int altclabel;
      int altclabelNode;


      if(depth < K-1 || ALT_Context[depth][1] == Kmask){

        clabel = Label::updateLabelDESUAfter(vsub, currentVertex, depth, batchOperationsMatrix);
        clabelNode = igtrie.insertLabel(labelNode, clabel, Label::repDigits(depth));
        auxlabel = (auxlabel << Label::repDigits(depth)) | clabel;
      }

      if(depth < K-1 || OGN_Context[depth][1] == Kmask){
        altclabel = Label::updateLabelDESUBefore(vsub, currentVertex, depth, batchOperationsMatrix);
        altclabelNode = igtrie.insertLabel(altLabelNode, altclabel, Label::repDigits(depth));
        auxaltlabel = (auxaltlabel << Label::repDigits(depth)) | altclabel;
      }

      
      
      //Não vale a pena expandir os vizinhos se já estivermos no último vértice do subgrafo porque não vamos expandir mais

     

      
     
      if(depth < K-1){
       

      updateToExpand(I, depth, &auxNExpanded, &auxNForbiddenChanges,&auxNUsedNeighbourChanges, vertexOGNmask, &auxNOGNChanges, vertexALTmask, &auxNALTChanges);

        
        
      }

      //int prevNOGNComps = nOGNComps;
      //int prevNALTComps = nALTComps;

      

      //updateComponents(depth,currentVertex,batchOperationsMatrix);

      


      DESU(depth+1, I+1, auxNExpanded, auxNForbiddenChanges, auxNUsedNeighbourChanges, auxNOGNChanges, auxNALTChanges,clabelNode, auxlabel, altclabelNode, auxaltlabel);

      

      //nOGNComps = prevNOGNComps;
      //nALTComps = prevNALTComps;
      //sizeVext = auxNExpanded;


      for(int i = nForbiddenChanges; i < auxNForbiddenChanges; i++){
        int vertex = forbiddenChanges[i];
        forbidden[vertex] = false;
      }
      for(int i = nUsedNeighbourChanges; i < auxNUsedNeighbourChanges; i++){
        int vertex = usedNeighbourChanges[i];
        usedNeighbour[vertex] = false;
      }
      for(int i = nALTChanges; i < auxNALTChanges; i++){
        int vertex = ALT_COMPS_Changes[i].first;
        int mask = ALT_COMPS_Changes[i].second;
        ALT_COMPS[vertex] = mask;
      }
      for(int i = nOGNChanges; i < auxNOGNChanges; i++){
        int vertex = OGN_COMPS_Changes[i].first;
        int mask = OGN_COMPS_Changes[i].second;
        //printf("[%d] %d %d\n",Kmask,vertex,mask);
        OGN_COMPS[vertex] = mask;
      }

      forbidden[currentVertex] = false;

      //nForbiddenChanges = auxNForbiddenChanges;
      //nUsedNeighbourChanges = auxNUsedNeighbourChanges;
     
  
  }
 
  
}

bool Fase::checkForbidden(int a, int b, int ** batchOperationsMatrix){
  if(batchOperationsMatrix[a][b] == 2 || batchOperationsMatrix[a][b] == -2 || batchOperationsMatrix[b][a] == 2 || batchOperationsMatrix[b][a] == -2){
    //printf("\t\tFORBIDDEN!\n");
    return true;  
  } 
  return false;
}



void Fase::updateToExpand(int currentIndex, int depth, int * sizeVext, int * nForbiddenChanges, int * nUsedNeighbourChanges, int OGNmask, int * nOGNChanges, int ALTmask, int * nALTChanges){
  //Dou reset ao número de vértices expandidos nesta profundidade para "limpar" expansões de outras recursões

  int currentVertex = vextList[currentIndex];
  int *eExcl = graph->arrayNeighbours(currentVertex);
  int eExclNum = graph->numNeighbours(currentVertex);
  
  

  for(int i = 0; i < eExclNum; i++){
    bool alreadyConsidered = false;
    //ITERATING OVER NEIGHBOURS
    int currNeighbour = eExcl[i];

    
    if(forbidden[currNeighbour]){
      continue;
    }
    if(checkForbidden(currNeighbour,currentVertex,batchOperationsMatrix)){
        //Verifico logo se o nó atual e o vizinho já são arestas proibidas
        
        forbidden[currNeighbour] = true;
        forbiddenChanges[(*nForbiddenChanges)++] = currNeighbour;
        //if((*nForbiddenChanges) > 29) printf("[%d/%d]\n",(*nForbiddenChanges),graph->numNodes());
        continue;
    }

    if(usedNeighbour[currNeighbour]){
      if(vextIndex[currNeighbour] > currentIndex){
        if(checkALTConnected(currentVertex,currNeighbour,batchOperationsMatrix)){
          //Maybe test if zero first so you just need to = and not add to changes(?)
          if(ALTmask | ALT_COMPS[currNeighbour] != ALT_COMPS[currNeighbour]){
            //ALT_COMPS_Changes[(*nALTChanges)].first = currNeighbour;
            //ALT_COMPS_Changes[(*nALTChanges)++].second = ALT_COMPS[currNeighbour];
            ALT_COMPS_Changes[(*nALTChanges)++] = make_pair(currNeighbour,ALT_COMPS[currNeighbour]);
            ALT_COMPS[currNeighbour] |= ALTmask;
          }else{
            printf("!!\n");
          }
        }
        if(checkOGNConnected(currentVertex,currNeighbour,batchOperationsMatrix)){
          //Maybe test if zero first so you just need to = and not add to changes(?)
          if(OGNmask | OGN_COMPS[currNeighbour] != OGN_COMPS[currNeighbour]){
            //OGN_COMPS_Changes[(*nOGNChanges)].first = currNeighbour;
            //OGN_COMPS_Changes[(*nOGNChanges)++].second = OGN_COMPS[currNeighbour];
            OGN_COMPS_Changes[(*nOGNChanges)++] = make_pair(currNeighbour,OGN_COMPS[currNeighbour]);
            OGN_COMPS[currNeighbour] |= OGNmask;
          }else{
            printf("??\n");
          }
        }
      }
      continue;
    }
    vextIndex[currNeighbour] = (*sizeVext);
    vextList[(*sizeVext)++] = currNeighbour;
    usedNeighbour[currNeighbour] = true;
    usedNeighbourChanges[(*nUsedNeighbourChanges)++] = currNeighbour;
    if(checkALTConnected(currentVertex,currNeighbour,batchOperationsMatrix)){
      ALT_COMPS[currNeighbour] = ALTmask;
    }else{
      ALT_COMPS[currNeighbour] = 0;
    }
    if(checkOGNConnected(currentVertex,currNeighbour,batchOperationsMatrix)){
      OGN_COMPS[currNeighbour] = OGNmask;
    }else{
      OGN_COMPS[currNeighbour] = 0;
    }

  }
 
  

}

bool Fase::checkIn(vector<pair<int,int>> * operations, int v1, int v2){
  for(int i = 0; i < (*operations).size(); i++){
    pair<int,int> operation = (*operations)[i];
    if(operation.first == v1 && operation.second == v2) return true;
    if(!directed){
      if(operation.second == v1 && operation.first == v2) return true;
    }
  }
  return false;
}

bool Fase::checkOGNConnected(int a, int b, int ** batchOperationsMatrix){
  return (graph->hasEdge(a,b) && batchOperationsMatrix[a][b] <= 0) || (graph->hasEdge(b,a) && batchOperationsMatrix[b][a] <= 0);
  
}

bool Fase::checkALTConnected(int a, int b, int ** batchOperationsMatrix){
  return (graph->hasEdge(a,b) && batchOperationsMatrix[a][b] >= 0) || (graph->hasEdge(b,a) && batchOperationsMatrix[b][a] >= 0);
  
}

/*
void Fase::updateComponents(int depth, int vertex, int ** batchOperationsMatrix){
  

  int newComponent = depth;

  memset(mergeOGNComp,false,sizeof(bool) * depth);
  memset(mergeALTComp,false,sizeof(bool) * depth);

  nOGNComps++;
  nALTComps++;


  for(int i = 0; i < depth; i++){
    int v = vsub[i];
    int compOGN = OGN_COMPS[depth-1][i];
    int compALT = ALT_COMPS[depth-1][i];
   
    if(!mergeOGNComp[compOGN]){
      if(checkOGNConnected(v,vertex,batchOperationsMatrix)){
        mergeOGNComp[compOGN] = true;
        nOGNComps--;
      }
      
    }

    //ALT_COMPS
    if(!mergeALTComp[compALT]){
      if(checkALTConnected(v,vertex,batchOperationsMatrix)){
        mergeALTComp[compALT] = true;
        nALTComps--;
        
      }
    }
  }


  for(int i = 0; i < depth; i++){
    int compOGN = OGN_COMPS[depth-1][i];
    int compALT = ALT_COMPS[depth-1][i];
    if(mergeOGNComp[compOGN]) OGN_COMPS[depth][i] = newComponent;
    else OGN_COMPS[depth][i] = compOGN;
    if(mergeALTComp[compALT]) ALT_COMPS[depth][i] = newComponent;
    else ALT_COMPS[depth][i] = compALT;
  }
  
  OGN_COMPS[depth][depth] = newComponent;
  ALT_COMPS[depth][depth] = newComponent;

}
*/


/*-------------------------------------------------------------------------------------------
|
|
|
|                       Código para obter contagens e subgrafos para dar output
|
|
---------------------------------------------------------------------------------------------*/






void Fase::getSubgraphFrequency(pair<long long int, int> element, Isomorphism* iso)
{
  Label::fillNautyMatrix(sadjM, K, element.first);

  nauty_s[0] = '\0';
  iso->canonicalStrNauty(sadjM, nauty_s);
  string str = string(nauty_s);
  canonicalTypes[str] += element.second;
}

void Fase::reduceCanonicalTypes()
{
  if (!canonicalTypes.empty())
    return;

  Isomorphism *iso = new Isomorphism();
  iso->initNauty(K, directed);
  for (auto element : igtrie.enumerate(K))
    getSubgraphFrequency(element, iso);
  iso->finishNauty();
  delete iso;
}

int Fase::getTypes()
{
 // canonicalTypes.clear();
 // reduceCanonicalTypes();
  return igtrie.getTypes();
}

vector<pair<int, string> > Fase::subgraphCount()
{
  //reduceCanonicalTypes();

  vector<pair<int, string> > subgraphVector;
  igtrie.subgraphCount(&subgraphVector);
  /*for (auto element : canonicalTypes)
    subgraphVector.push_back(make_pair(element.second, element.first));
*/
  sort(subgraphVector.begin(), subgraphVector.end());
  reverse(subgraphVector.begin(), subgraphVector.end());

  return subgraphVector;
}
