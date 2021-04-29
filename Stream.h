#ifndef _STREAM_
#define _STREAM_

#include "Graph.h"
#include "Fase.h"
#include "Fase.h"
#include "DynamicGraph.h"
#include "GraphMatrix.h"
#include "Timer.h"
#include "Random.h"

class Stream {
	protected:
		//char * _streamFile;
		Graph * graph;
		bool samp;
		FILE * streamFile;
		bool directed, streaming;
		double * sampProb;
		char * streamType;
		int ** batchOperationsMatrix; // -1 -> Edge to be deleted -2 -> Edge was deleted 1 -> Edge to be added 2 -> Edge was added
		pair<int,int> * additions;
		pair<int,int> * deletions;
		double prob;
		int base;
		vector<tuple<char,int,int>> updateList;
		int K;
		Fase * fase;
		FILE * outFile;
		Timer * global_timer;
		Timer * local_timer;
		Timer * iso_timer;
		Timer * preprocess_timer;
		bool debugprints;
		char * sfilename;
		char * ofilename;
		void _readFile();
		void initSamplingProbabilities();
		void finish();
		void singleCensus();
		void output();
		void addEdge(int a, int b);
		void removeEdge(int a, int b);
		virtual int parseUpdate(vector<tuple<char,int,int>> * updateList, int i, int batchSize) = 0;
		//time_t t_start, t_end;
		//void output();

	public:
		Stream();
		~Stream();
		void runCensus(int K, int batchSize,bool debugprints);
		//void printFilename() {printf("Well N nodes = %d and sfilename = %c\n",graph->numNodes(),sfilename[0]);}
		void createStream(Graph * _g, bool samp, double * _sampProb, bool _dir, char * _streamF, char * _outFile, bool _streaming, int _base);
		

};


#endif