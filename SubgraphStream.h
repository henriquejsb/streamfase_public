#ifndef _SUBGRAPHSTREAM_
#define _SUBGRAPHSTREAM_


#include "Stream.h"

class SubgraphStream : public Stream {
	private:
		int parseUpdate(vector<tuple<char,int,int>> * updateList, int i, int batchSize);
	
		

	public:
		SubgraphStream();
		
		
};


#endif