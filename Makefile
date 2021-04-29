EXEC_NAME=FASE
CC=g++

CFLAGS= -w -Wno-write-strings -O2 -std=c++11 -g  
#CFLAGS = -Wall -Wno-write-strings -O1 -g -std=c++11 -g
#CFLAGS = -Wall -Wno-write-strings  -g -std=c++11 -g


SRC =                   \
	Fase.cpp 	\
	Label.cpp	\
	IGtrie.cpp	\
	Isomorphism.cpp	\
	Timer.cpp	\
	DynamicGraph.cpp\
	GraphMatrix.cpp	\
	GraphUtils.cpp	\
	Random.cpp	\
	nauty/nauty.c	\
	nauty/nautil.c	\
	nauty/naugraph.c \
	Stream.cpp \
	RecalculateStream.cpp \
	SubgraphStream.cpp \
	DESUStream.cpp \
	FaseMain.cpp 

OBJ =  ${SRC:.cpp=.o}

#------------------------------------------------------------

all: ${EXEC_NAME}

${EXEC_NAME}: ${OBJ}
	${CC} ${CFLAGS} ${CLIBS} -o ${EXEC_NAME} ${OBJ}

%.o: %.cpp
	${CC} ${CFLAGS} -c -o $@ $+

clean:
	rm ${EXEC_NAME} *.o *~ *# -rf
