
OBJECTS=RefReader.o RefGenome.o DnaBitStr.o main.o
PROGNAME=Metal
CXX=g++

CXXFLAGS= -std=c++0x -ggdb -Wall -pedantic -pipe -fopenmp

.PHONY: all clean

all: ${PROGNAME}

# test: testsuite

%.o: %.cpp %.h
	${CXX} ${CXXFLAGS} -c $<

%.o: %.cpp
	${CXX} ${CXXFLAGS} -c $<

# testsuite: ${OBJECTS} main_test.o
# 	${CXX} ${CXXFLAGS} ${OBJECTS} main_test.o -o $@ 

${PROGNAME}: ${OBJECTS}
	${CXX} ${CXXFLAGS} ${OBJECTS} -o $@ 

clean:
	rm -f ${OBJECTS} ${PROGNAME} main_test.o
