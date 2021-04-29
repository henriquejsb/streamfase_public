#ifndef _DESUSTREAM_
#define _DESUSTREAM_

#include "Stream.h"

class DESUStream : public Stream {
	private:
			
		int parseUpdate(vector<tuple<char,int,int>> * updateList, int i, int batchSize);
		
		
	

	public:
		DESUStream();
		
		
};


#endif