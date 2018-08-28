
OBJECTS=main.o Extract.o
PROGNAME=Extractor
CXX=g++

CXXFLAGS= -std=c++0x -ggdb -Wall -pedantic -pipe -O3 -fopenmp

.PHONY: all clean

all: ${PROGNAME}

%.o: %.cpp %.h
	${CXX} ${CXXFLAGS} -c $<

%.o: %.cpp
	${CXX} ${CXXFLAGS} -c $<

${PROGNAME}: ${OBJECTS}
	${CXX} ${CXXFLAGS} ${OBJECTS} -o $@ 

clean:
	rm -f ${OBJECTS} ${PROGNAME}
