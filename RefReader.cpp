#include <stdio.h>
#include <cstring>

#include "RefReader.h"

using namespace std;

void readReference(char const * filename, vector<const struct CpG>* cpgtab, vector<const string>* genSeq)
{

    // reserve 24 Mil. entries (roughly more then CpGs in human genome)
    cpgtab.reserve(24000000);
    genSeq.reserve(24);

    // stores the length of the sequence read so far (i.e. invariant holds:
    // read characters so far is seq_length for the while loop)
    int seq_length = 0;

    // number of CpGs in whole sequence
    int cpg_num = 0;

    // chromosome index
    // used as index for genSeq
    int chr_index = -1;

    // chromosome string
    string chromSeq;
    // reserve length of 260 Mil. (~240Mil is number of bps in largest human chr.)
    // done to avoid constant reallocations during input read
    chromSeq.reserve(260000000);


    // buffer for IO
    char fileBuf[BUFSIZE + 1];

    // open file descriptor
    int file = open(filename, O_RDONLY);

    // unsuccessfull
    if (file == -1)
    {

        cerr << "Error opening file " << filename << "\n\n";
        exit(EXIT_FAILURE);
    }

    // tell the kernel that data is read sequentially to allow optimizations
    if (posix_fadvice(file, 0, 0, POSIX_FADV_SEQUENTIAL) != 0)
    {

        cerr << "Error emitting access pattern for " << filename << " with file descriptor " << file << " to kernel\n\n";
        exit(EXIT_FAILURE);
    }

    // flag which states if we are currently reading part of a chromosome sequence (true) or unmapped stuff (false)
    bool cont_flag = false;

    size_t charsRead = -1;

    // read to buffer BUFSIZE many chars
    for (charsRead = read(file, fileBuf, BUFSIZE); charsRead > (size_t) 0; charsRead = read(file, fileBuf, BUFSIZE))
    {

        // contains positions of all '>' in buffer
        list<const struct id> idQueue;

        // read buffer until next id line
        for(idLine = fileBuf; (idLine = (char*) memchr(idLine, '>', (fileBuf + charsRead) - idLine)); ++idLine)
        {

            char* startSec = ++((char*) memchr(idLine, '\n', (fileBuf + charsRead) - idLine));
            // check if id line corresponds to primary assembly of a chromosome
            // (i.e. following letter is a C)
            if (*(idLine + 1) == 'C')
            {

                idQueue.emplace_back({startSec, idLine, true});


            } else {

                idQueue.emplace_back({startSec, idLine, false});
            }
        }
        if (!cont_flag) {

            while (!idQueue.isEmpty()) {

                // TODO Wrap the whole line reading into this.
                //      use offsets to read only allowed part of the buffer (i.e. assembly sequence id tag
            }
        }

        // find the first newline
        char* prevNL;
        if ( (prevNL =(char*) memchr(fileBuf, '\n', charsRead)) )
        {

            // append chars until newline
            chromSeq.append(fileBuf, prevNL - fileBuf);

        } else {

            // append all chars
            chromSeq.append(fileBuf, charsRead);
            // no further newline to find continue reading to buffer
            continue;

        }

        // read buffer until next newline, append this part to sequence
        for(char * line = fileBuf; (line = (char*) memchr(line, '\n', (fileBuf + charsRead) - line)); ++line)
        {

            // append the line to sequence
            chromSeq.append(fileBuf, line - prevNL);

            prevNL = line;
        }

        // append rest of buffer
        chromSeq.append(fileBuf, fileBuf + charsRead - prevNL);
    }

    if (charsRead < 0)
    {

        cerr << "Error reading file " << filename << "\n\n";
        exit(EXIT_FAILURE);
    }





    // TODO
    while (getline(ifs, line)) {

        // Throw out unlocalized contigs
        if (line.begin() != line.end() && *(line.begin()) == '>' && ++(line.begin()) != line.end() && *(++(line.begin())) != 'C')
        {

            cont_flag = false;
            continue;
        }
        // if we have description tag of real chromosome assembly
        if (line.begin() != line.end() && *(line.begin()) == '>' && ++(line.begin()) != line.end() && *(++(line.begin())) == 'C')
        {

            // if we have read something already read, dump it to the output reset chromSeq variable
            if (chr_index > -1)
            {

                // shrink size to actual sequence length
                chromSeq.shrink_to_fit();
                // dump chromosome sequence to output table
                genSeq[chr_index] = chromSeq;
                // keeps capacity untouched and erases content
                chromSeq.clear();
                // extend capacity again
                chromSeq.reserve(260000000);
            }
            ++chr_index;
            cont_flag = true;
            continue;
        }
        // ignore descriptions and comments
        if (!cont_flag || line.begin() == line.end() || *(line.begin()) == ';' )
        {

            continue;
        }

        for (string::iterator it = line.begin(); it != line.end(); ++it)
        {
            ++seq_length;

            // pop Cs that are out of the window
            if (!c_queue.empty())
            {

                if  (c_queue.back() < (seq_length - 2*READL))
                {
                    --c_count;
                    c_queue.pop_back();
                }
            }
            // pop cpgs that are out of the window, update statistics
            if (!cpg_queue.empty())
            {
                if  (cpg_queue.back() < (seq_length - READL + 1))
                {

                    --cpg_count;
                    ++cpg_adj_c[c_count];
                    ++cpg_adj_cpg[cpg_count];
                    cpg_queue.pop_back();
                }
            }


            if (*it == 'C' || *it == 'c')
            {

                ++c_count;
                c_queue.push_front(seq_length);

                if ( (++it) != line.end() )
                {
                    ++seq_length;

                    if (*it == 'G' || *it == 'g')
                    {

                        ++cpg_num;
                        ++cpg_count;
                        cpg_queue.push_front(seq_length - 1);
                    }

                } else {

                    break;
                }
            }
        }
    }

    // discard unused memory
    cpgtab.shrink_to_fit();

}
