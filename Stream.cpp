#include "Stream.h"
#include <stdio.h>
#include <algorithm>

Stream::Stream(){

}

Stream::~Stream(){
	
}


void Stream::createStream(Graph * _g, bool _samp,  double * _sampProb, bool _dir, char * _streamF, char * _outFile, bool _streaming, int _base, bool _detailed){
  
  graph = _g;
  directed = _dir;
  if (_streamF[0] == '0' && _streamF[1] == '\0')
      streamFile = NULL;
    else
      streamFile = fopen(_streamF, "r");  
  sampProb = _sampProb;
  samp = _samp;
  base = _base;
  streaming = _streaming;
  detailed = _detailed;
  fase = NULL;
  K = 0;
  if (_outFile[0] == '0' && _outFile[1] == '\0')
      outFile = stdout;
    else
      outFile = fopen(_outFile, "w");
  global_timer = new Timer();
  local_timer = new Timer();
  iso_timer = new Timer();
  preprocess_timer = new Timer();
  batchOperationsMatrix = NULL;
  debugprints = false;
}


void Stream::initSamplingProbabilities()
{
  int i;
  prob = 1.0;
  for (i = 0; i < K; i++)
    prob *= sampProb[i];

  if (samp && fabs(prob) > 10e-7)
    fase->initSampling(K, sampProb);
  else
    prob = 1.0;
}


void Stream::finish()
{
  delete global_timer;
  delete local_timer;
  delete iso_timer;
  delete preprocess_timer;
  if(streamType != "RecalculateStream") delete fase;
  if(batchOperationsMatrix != NULL){

    for(int i = 0; i < graph->numNodes(); i++){
      delete [] batchOperationsMatrix[i];
    }
    delete [] batchOperationsMatrix;
    delete [] additions;
    delete [] deletions;
  }
  delete graph;
  fclose(outFile);
}

void Stream::_readFile(){
	if (!streamFile)
    	exit(EXIT_FAILURE);
  char c;
  int a,b;
	while(fscanf(streamFile,"%c %d %d\n",&c,&a,&b) != EOF){
    tuple<char,int,int> newUpdate = make_tuple(c,a,b);
		updateList.push_back(newUpdate);
	}
	fclose(streamFile);

}

void Stream::addEdge(int a, int b){
  graph->addEdge(a,b);
  if (!directed) graph->addEdge(b,a);
}

void Stream::removeEdge(int a, int b){
  graph->rmEdge(a,b);
  if (!directed) graph->rmEdge(b,a);
}
  
void Stream::singleCensus(){
  fase = new Fase(graph, K, directed,debugprints);
  fase->runCensus();
}

void Stream::output()
{
  FILE * f = outFile;
  //fprintf(f,"Finished Calculating\n");
  //fprintf(f, "\tOutput:\n");
  //fprintf(f, "Network: %s\n", ifilename);
  //fprintf(f, "Directed: %s\n", directed ? "Yes" : "No");
  //fprintf(f, "Nodes: %d\n",graph->numNodes());
  //fprintf(f, "Edges: %d\n", graph->numEdges() / (directed ? 1 : 2));
 // fprintf(f, "Subgraph Size: %d\n", K);
 /* if (largeScale)
    fprintf(f, "Graph Representation: Large Scale\n");
  else
    fprintf(f, "Graph Representation: Adjacency Matrix\n");
  */
  /*
  t_end = time(0);
  struct tm *tm_start = localtime(&t_start);
  fprintf(f, "Start of Computation: %02dh%02dm%02ds %02d/%02d/%02d\n\
", tm_start->tm_hour, tm_start->tm_min, tm_start->tm_sec, tm_start->tm_mday, tm_start->tm_mon + 1, 1900 + tm_start->tm_year);
  struct tm *tm_end   = localtime(&t_end);
  fprintf(f, "End of Computation: %02dh%02dm%02ds %02d/%02d/%02d\n", tm_end->tm_hour, tm_end->tm_min, tm_end->tm_sec, tm_end->tm_mday, tm_end->tm_mon + 1, 1900 + tm_end->tm_year);
  */
  //fprintf(f, "\n\n\tResults:\n");
  fprintf(f, "Subgraph Occurrences: %lld\n",(long long int)(fase->getMotifCount() ));
  //fprintf(f, "Subgraph Types: %d\n", fase->getTypes());
  fprintf(f, "Global Computation time (ms): %0.4lf\n", global_timer->elapsed());
  fprintf(f, "Computation Time (ms): %0.6lf\n", local_timer->elapsed());

  if(detailed){
    
  
    fprintf(f, "\n\tDetailed Output:\n");
    for (auto element : fase->subgraphCount())
        fprintf(f, "%s: %d occurrences\n", element.second.c_str(), element.first);
 }
}


void Stream::runCensus(int _K, int batchSize, bool _debugprints){
 
  K = _K;
  global_timer->start();
  local_timer->start();
  debugprints = _debugprints;
  //printf("Hello!\n");
  singleCensus();
  
  local_timer->stop();
  global_timer->stop();
  output();
  //TODO ????? 
  batchOperationsMatrix = new int*[graph->numNodes()];
  for(int i = 0; i < graph->numNodes(); i++){
    batchOperationsMatrix[i] = new int[graph->numNodes()];
    memset(batchOperationsMatrix[i],0, sizeof(int) * graph->numNodes());
  }
  additions = new pair<int,int>[batchSize];
  deletions = new pair<int,int>[batchSize];

  if (streaming){
    if (streamFile){
     _readFile();

      for(int i = 0; i < updateList.size(); i++){
        //fprintf(outFile,"\n\n---------------\n[%d/%d] %c %d %d\n---------------\n",i+1,updateList.size(),get<0>(updateList[i]),get<1>(updateList[i]),get<2>(updateList[i]));
        local_timer->start();
        int fail = parseUpdate(&updateList,i,min(batchSize,(int) updateList.size() - i));
        i += batchSize - 1;
        local_timer->stop();
        global_timer->stop();
        //printf("Computation Time (ms): %0.6lf\n", local_timer->elapsed());
        if(fail) continue;
        output();
        if(streamType == "RecalculateStream") delete fase;
      }
    }
    else{
      char c;
      int a,b;
      fflush(stdin);
      while (scanf("%c %d %d\n",&c,&a,&b)){
        //parseUpdate(c,a,b);
        fflush(stdin);
        local_timer->stop();
        global_timer->stop();
        output();
       
      }
    }
  } 
  global_timer->stop();
  //printf("Finished!\n");
  //output();
  finish();
}
