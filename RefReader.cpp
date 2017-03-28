#include <stdio.h>
#include <cstring>      // memchr()
#include <list>
#include <fcntl.h>      // posix_fadvise()
#include <iostream>
#include <unistd.h>     // read()

#include "RefReader.h"
#include "CONST.h"

using namespace std;

void readReference(const char* const filename, vector<struct CpG>& cpgTab, vector<const char*>& genSeq, vector<size_t> genSeqLen)
{

    // reserve CPGMAX many entries (should be roughly more than CpGs in human genome)
    cpgTab.reserve(MyConst::CPGMAX);
    genSeq.reserve(MyConst::CHROMNUM);

    // stores the length of the sequence read so far
    unsigned int seqLength = 0;

    // index which chromosome is read to access genSeq
    uint8_t chrIndex = 0;

    // chromosome string
    string chromSeq;
    // reserve length of CHROMMAX many entries (should be roughly more than number of bps in biggest chromosome)
    // done to avoid constant reallocations during input read
    chromSeq.reserve(MyConst::CHROMMAX);


    // buffer for IO
    char fileBuf[MyConst::BUFSIZE + 1];

    // open file descriptor
    int file = open(filename, O_RDONLY);

    // unsuccessfull
    if (file == -1)
    {

        cerr << "Error opening file " << filename << "\n\n";
        exit(EXIT_FAILURE);
    }

    // tell the kernel that data is read sequentially to allow optimizations
    if (posix_fadvise(file, 0, 0, POSIX_FADV_SEQUENTIAL) != 0)
    {

        cerr << "Error emitting access pattern for " << filename << " with file descriptor " << file << " to kernel\n\n";
        exit(EXIT_FAILURE);
    }

    // flag which states if we are currently reading part of a chromosome sequence (true) or unmapped stuff (false)
    bool contFlag = false;

    // flag which states if the last buffer ended with a C
    bool cEndFlag = false;

    size_t charsRead = read(file, fileBuf, MyConst::BUFSIZE);

    // read to buffer BUFSIZE many chars
    for ( ; charsRead > (size_t) 0; charsRead = read(file, fileBuf, MyConst::BUFSIZE) )
    {

        // contains positions of all '>' in buffer
        list<struct idPos> idQueue;
        extractIdLines(fileBuf, charsRead, idQueue);


        // if no new sequence id tags, read through
        if (idQueue.empty())
        {

            if (contFlag)
            {
                // TODO: USE INFORMATION OF CpGs ACROSS DIFFERENT BUFFERS!
                constructCpgs(fileBuf, charsRead, chrIndex, seqLength, cpgTab, cEndFlag);
                seqLength += readBufferSlice(fileBuf, charsRead, chromSeq);
            }
            // read new buffer
            continue;

        } else {

            // if we are in primary assembly, read slice until next id tag
            if (contFlag)
            {

                constructCpgs(fileBuf, idQueue.front().id - fileBuf, chrIndex, seqLength, cpgTab, cEndFlag);
                seqLength += readBufferSlice(fileBuf, idQueue.front().id - fileBuf, chromSeq);

                chromSeq.shrink_to_fit();
                // save chromosome sequence
                genSeq.push_back(chromSeq.c_str());
                genSeqLen.push_back(chromSeq.size());
                // make new sequence buffer
                chromSeq.clear();
                chromSeq.reserve(MyConst::CHROMMAX);
                ++chrIndex;

            }
        }

        //  where to start reading buffer
        char * primStart;
        unsigned int offset;

        while (!idQueue.empty()) {

            // if the next id tag is still no primary assembly skip this blog
            if (!idQueue.front().imp)
            {
                // if we read primary assembly previously, save it
                if (contFlag)
                {
                    chromSeq.shrink_to_fit();
                    // save chromosome sequence
                    genSeq.push_back(chromSeq.c_str());
                    genSeqLen.push_back(chromSeq.size());
                    // make new sequence buffer
                    chromSeq.clear();
                    chromSeq.reserve(MyConst::CHROMMAX);
                    ++chrIndex;
                }

                idQueue.pop_front();
                contFlag = false;
                continue;

            // if we are in primary assembly
            } else {

                // set start to start of primary assembly section
                primStart = idQueue.front().sec;
                idQueue.pop_front();

                if (idQueue.empty())
                {
                    // set offset to buffer end
                    offset = fileBuf + charsRead - primStart;

                } else {

                    // set offset to next id tag
                    offset = idQueue.front().id - primStart;
                }
                constructCpgs(primStart, offset, chrIndex, seqLength, cpgTab, cEndFlag);
                seqLength += readBufferSlice(primStart, offset, chromSeq);
                contFlag = true;
            }

        }

    }

    if (charsRead < 0)
    {

        cerr << "Error reading file " << filename << "\n\n";
        exit(EXIT_FAILURE);
    }
}


inline unsigned int readBufferSlice(char* const start, const unsigned int offset, std::string& chromSeq)
{

    unsigned int savedChars = 0;
    // stores previous newline + 1 (i.e. beginning of next line)
    char* prevNL;
    // find the first newline
    if ( (prevNL =(char*) memchr(start, '\n', offset)) )
    {

        // append chars until newline
        chromSeq.append(start, prevNL - start);
        savedChars += prevNL - start;

    } else {

        // append all chars
        chromSeq.append(start, offset);
        // no further newline to find continue reading to buffer
        return offset;

    }
    ++prevNL;

    // read buffer until next newline, append this part to sequence
    for(char * line = start; (line = (char*) memchr(line, '\n', (start + offset) - line)); )
    {

        // append the line to sequence
        chromSeq.append(prevNL, line - prevNL);
        savedChars += line - prevNL;

        ++line;
        prevNL = line;
    }

    // append rest of buffer
    chromSeq.append(prevNL, start + offset - prevNL);
    savedChars += start + offset - prevNL;
    return savedChars;
}


bool constructCpgs(char* const start, const unsigned int offset, const uint8_t chrIndex, const unsigned int SeqLength, std::vector<struct CpG>& cpgTab, bool& cEndFlag)
{

    // TODO: skip newline from count...
    //
    // where to set the start of the CpG
    constexpr unsigned int winStart = MyConst::READLEN + 2;

    // if previous char was c, check if first char is g
    if (cEndFlag)
    {
        if (*start == 'G' || *start == 'g')
        {

            // another off by one because we are already at the 'G'
            cpgTab.emplace_back(chrIndex, SeqLength - winStart - 1);

        }
    }

    // find Cs
    for (char* cBase = start; (cBase = (char*) memchr(cBase, 'C', (start + offset - 1) - cBase)); ++cBase)
    {

        // if next character is G, construct CpG
        if (*(cBase + 1) == 'G' || *(cBase + 1) == 'g')
        {

            cpgTab.emplace_back(chrIndex, SeqLength + (cBase - start) - winStart);
        }

    }
    // find cs (lowercase mask for repeats)
    for (char* cBase = start; (cBase = (char*) memchr(cBase, 'c', (start + offset - 1) - cBase)); ++cBase)
    {

        // if next character is G, construct CpG
        if (*(cBase + 1) == 'G' || *(cBase + 1) == 'g')
        {

            cpgTab.emplace_back(chrIndex, SeqLength + (cBase - start) - winStart);
        }

    }

    // check last char
    if (*(start + offset) == 'C' || *(start + offset) == 'c')
    {
        cEndFlag = true;
    }

    cEndFlag = false;
}
