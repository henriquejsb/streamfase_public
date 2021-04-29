#ifndef _RECALCULATESTREAM_
#define _RECALCULATESTREAM_

#include "Stream.h"

class RecalculateStream : public Stream {
	private:
			
		int parseUpdate(vector<tuple<char,int,int>> * updateList, int i, int batchSize);
		
		
	

	public:
		RecalculateStream();
		
		
};


#endif