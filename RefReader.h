#ifndef REFREADER_H
#define REFREADER_H

#include <string>
#include <vector>


#include "structs.h"

// read reference FASTA file from ifs
// produces a vector of all CpGs present in the genome written to cpgTab
// produces sequence strings seperated by chromosome saved to genSeq, their length to genSeqLen
//      underlying vectors should be empty on calling
//      convention: table index 0-21  autosome 1-22, 22-23 allosome X,Y
void readReference(const char* const filename, std::vector<struct CpG>& cpgTab, std::vector<const char*>& genSeq, std::vector<std::size_t> genSeqLen);


// internal function to read buffer slice corresponding to primary assembly
// appends read letters to chromSeq
// starts reading buffer at start ending at start+offset
// returns number of characters appended to chromSeq
unsigned int readBufferSlice(char* const start, const unsigned int offset, std::string& chromSeq);


// internal function to construct CpG structs from given char array
// fills cpgTab with CpG structs found in array, using chrIndex and seqLength (number of bps read before this slice in this chromosome) as position information
// cEndFlag states if the last char of previous buffer was a C
void constructCpgs(const char* const  start, const unsigned int offset, const uint8_t chrIndex, const unsigned int SeqLength, std::vector<struct CpG>& cpgTab, bool& cEndFlag );


// extracts (the position offsets of) all lines that are id tags (starting with '>')
// and saves them in idQueue
// Arguments:
//              start:  where to start the search
//              offset: how many chars to read
//              idQueue:where to save the result
inline void extractIdLines(const char* const start, const unsigned int offset, std::list<struct idPos>& idQueue)
{
    // extract id lines from buffer
    for(char* idLine = start; (idLine = (char*) memchr(idLine, '>', (start + offset) - idLine)); )
    {

        char* startSec = ((char*) memchr(idLine, '\n', (start + offset) - idLine)) + 1;
        // check if id line corresponds to primary assembly of a chromosome
        // (i.e. following letter is a C)
        if (*(idLine + 1) == 'C')
        {

            idQueue.emplace_back(startSec, idLine, true);


        } else {

            idQueue.emplace_back(startSec, idLine, false);
        }
        idLine = startSec;
    }
}


#endif /* REFREADER_H */
