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
    bool cont_flag = false;

    size_t charsRead = read(file, fileBuf, MyConst::BUFSIZE);

    // read to buffer BUFSIZE many chars
    for ( ; charsRead > (size_t) 0; charsRead = read(file, fileBuf, MyConst::BUFSIZE) )
    {

        // contains positions of all '>' in buffer
        list<struct idPos> idQueue;

        // read buffer until next id line
        for(char* idLine = fileBuf; (idLine = (char*) memchr(idLine, '>', (fileBuf + charsRead) - idLine)); )
        {

            char* startSec = ((char*) memchr(idLine, '\n', (fileBuf + charsRead) - idLine)) + 1;
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

        // if no new sequence id tags, read through
        if (idQueue.empty())
        {

            if (cont_flag)
            {
                constructCpgs(fileBuf, charsRead, chrIndex, seqLength, cpgTab);
                seqLength += readBufferSlice(fileBuf, charsRead, chromSeq);
            }
            // read new buffer
            continue;

        } else {

            // if we are in primary assembly, read slice until next id tag
            if (cont_flag)
            {

                constructCpgs(fileBuf, idQueue.front().id - fileBuf, chrIndex, seqLength, cpgTab);
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
                if (cont_flag)
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
                cont_flag = false;
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
                constructCpgs(primStart, offset, chrIndex, seqLength, cpgTab);
                seqLength += readBufferSlice(primStart, offset, chromSeq);
                cont_flag = true;
            }

        }

    }

    if (charsRead < 0)
    {

        cerr << "Error reading file " << filename << "\n\n";
        exit(EXIT_FAILURE);
    }





    // TODO
    // while (getline(ifs, line)) {
    //
    //     // Throw out unlocalized contigs
    //     if (line.begin() != line.end() && *(line.begin()) == '>' && ++(line.begin()) != line.end() && *(++(line.begin())) != 'C')
    //     {
    //
    //         cont_flag = false;
    //         continue;
    //     }
    //     // if we have description tag of real chromosome assembly
    //     if (line.begin() != line.end() && *(line.begin()) == '>' && ++(line.begin()) != line.end() && *(++(line.begin())) == 'C')
    //     {
    //
    //         // if we have read something already read, dump it to the output reset chromSeq variable
    //         if (chrIndex > -1)
    //         {
    //
    //             // shrink size to actual sequence length
    //             chromSeq.shrink_to_fit();
    //             // dump chromosome sequence to output table
    //             genSeq[chrIndex] = chromSeq;
    //             // keeps capacity untouched and erases content
    //             chromSeq.clear();
    //             // extend capacity again
    //             chromSeq.reserve(260000000);
    //         }
    //         ++chrIndex;
    //         cont_flag = true;
    //         continue;
    //     }
    //     // ignore descriptions and comments
    //     if (!cont_flag || line.begin() == line.end() || *(line.begin()) == ';' )
    //     {
    //
    //         continue;
    //     }
    //
    //     for (string::iterator it = line.begin(); it != line.end(); ++it)
    //     {
    //         ++seq_length;
    //
    //         // pop Cs that are out of the window
    //         if (!c_queue.empty())
    //         {
    //
    //             if  (c_queue.back() < (seq_length - 2*READL))
    //             {
    //                 --c_count;
    //                 c_queue.pop_back();
    //             }
    //         }
    //         // pop cpgs that are out of the window, update statistics
    //         if (!cpg_queue.empty())
    //         {
    //             if  (cpg_queue.back() < (seq_length - READL + 1))
    //             {
    //
    //                 --cpg_count;
    //                 ++cpg_adj_c[c_count];
    //                 ++cpg_adj_cpg[cpg_count];
    //                 cpg_queue.pop_back();
    //             }
    //         }
    //
    //
    //         if (*it == 'C' || *it == 'c')
    //         {
    //
    //             ++c_count;
    //             c_queue.push_front(seq_length);
    //
    //             if ( (++it) != line.end() )
    //             {
    //                 ++seq_length;
    //
    //                 if (*it == 'G' || *it == 'g')
    //                 {
    //
    //                     ++cpg_num;
    //                     ++cpg_count;
    //                     cpg_queue.push_front(seq_length - 1);
    //                 }
    //
    //             } else {
    //
    //                 break;
    //             }
    //         }
    //     }
    // }

    // discard unused memory
    cpgTab.shrink_to_fit();

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


bool constructCpgs(char* const start, const unsigned int offset, const uint8_t chrIndex, const unsigned int SeqLength, std::vector<struct CpG>& cpgTab)
{

    // find Cs
    // last char cannot be a c of cpg in this buffer
    for (char* cBase = start; (cBase = (char*) memchr(cBase, 'C', (start + offset - 1) - cBase)); ++cBase)
    {

        // if next character is G, construct CpG
        if (*(cBase + 1) == 'G' || *(cBase + 1) == 'g')
        {

            cpgTab.emplace_back(chrIndex, SeqLength + (cBase - start + 1));
        }

    }
    // find cs (lowercase mask for repeats)
    for (char* cBase = start; (cBase = (char*) memchr(cBase, 'c', (start + offset - 1) - cBase)); ++cBase)
    {

        // if next character is G, construct CpG
        if (*(cBase + 1) == 'G' || *(cBase + 1) == 'g')
        {

            cpgTab.emplace_back(chrIndex, SeqLength + (cBase - start + 1));
        }

    }

    // check last char
    if (*(start + offset) == 'C' || *(start + offset) == 'c')
    {
        return true;
    }

    return false;
}
