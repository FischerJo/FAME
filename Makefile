
OBJECTS=gzstream.o RefReader_istr.o RefGenome.o DnaBitStr.o main.o\
		ReadQueue.o Read.o CONST.o ShiftAnd.o LevenshtDP.o
PROGNAME=Metal
CXX=g++

CXXFLAGS= -std=c++0x -ggdb -Wall -pedantic -pipe -O3 -fopenmp
GZFLAGS= -lz

.PHONY: all clean profile

all: ${PROGNAME}

profile: ${PROGNAME}Profile

gzstream.o: gzstream/gzstream.C gzstream/gzstream.h
	${CXX} ${CXXFLAGS} -c $<

ReadQueue.o: ReadQueue.cpp ReadQueue.h gzstream/gzstream.C
	${CXX} ${CXXFLAGS} -c $<

%.o: %.cpp %.h
	${CXX} ${CXXFLAGS} -c $<

%.o: %.cpp
	${CXX} ${CXXFLAGS} -c $<

${PROGNAME}: ${OBJECTS}
	${CXX} ${CXXFLAGS} ${OBJECTS} ${GZFLAGS} -o $@ 

${PROGNAME}Profile: ${OBJECTS}
	${CXX} ${CXXFLAGS} -pg -rdynamic ${OBJECTS} ${GZFLAGS} -o $@

clean:
	rm -f ${OBJECTS} ${PROGNAME}
