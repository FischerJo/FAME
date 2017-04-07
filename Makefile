
OBJECTS=RefReader_istr.o RefGenome.o DnaBitStr.o main.o
PROGNAME=Metal
CXX=g++

CXXFLAGS= -std=c++0x -ggdb -Wall -pedantic -pipe -O3

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
	${CXX} ${CXXFLAGS} -pg ${OBJECTS} -o $@

clean:
	rm -f ${OBJECTS} ${PROGNAME}
