#ifndef REFREADER_H
#define REFREADER_H

#include <string>


#include "structs.h"

// read reference FASTA file from ifs
// produces a vector of all CpGs present in the genome written to cpgtab
// produces sequence strings seperated by chromosome saved to genSeq
//      underlying vector should be empty on calling
//      convention: table index 0-21  autosome 1-22, 22-23 allosome X,Y
void readReference(const char * filename, std::vector<const struct CpG>* cpgtab, std::vector<const std::string>* genSeq);

#endif /* REFREADER_H */
