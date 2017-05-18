
OBJECTS=RefReader_istr.o RefGenome.o DnaBitStr.o main.o ReadQueue.o Read.o CONST.o
PROGNAME=Metal
CXX=g++

CXXFLAGS= -std=c++0x -ggdb -Wall -pedantic -pipe -O3 -fopenmp

.PHONY: all clean profile

all: ${PROGNAME}

profile: ${PROGNAME}Profile

%.o: %.cpp %.h
	${CXX} ${CXXFLAGS} -c $<

%.o: %.cpp
	${CXX} ${CXXFLAGS} -c $<

${PROGNAME}: ${OBJECTS}
	${CXX} ${CXXFLAGS} ${OBJECTS} -o $@ 

${PROGNAME}Profile: ${OBJECTS}
	${CXX} ${CXXFLAGS} -pg -rdynamic ${OBJECTS} -o $@

clean:
	rm -f ${OBJECTS} ${PROGNAME}
