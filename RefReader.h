#ifndef REFREADER_H
#define REFREADER_H

#include <string>
#include <vector>


#include "structs.h"

// read reference FASTA file from ifs
// produces a vector of all CpGs present in the genome written to cpgTab
// produces sequence strings seperated by chromosome saved to genSeq
//      underlying vector should be empty on calling
//      convention: table index 0-21  autosome 1-22, 22-23 allosome X,Y
void readReference(const char* const filename, std::vector<struct CpG>& cpgTab, std::vector<std::string>& genSeq);


// internal function to read buffer slice corresponding to primary assembly
// appends read letters to chromSeq
// starts reading buffer at start ending at start+offset
// returns number of characters appended to chromSeq
unsigned int readBufferSlice(char* const start, const unsigned int offset, std::string& chromSeq);


// internal function to construct CpG structs from given char array
// fills cpgTab with CpG structs found in array, using chrIndex and seqLength (number of bps read before this slice in this chromosome) as position information
// return true iff the last read char was a C
bool constructCpgs(char* const  start, const unsigned int offset, const uint8_t chrIndex, const unsigned int SeqLength, std::vector<struct CpG>& cpgTab);

#endif /* REFREADER_H */
