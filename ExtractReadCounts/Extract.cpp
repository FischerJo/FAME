
#include <list>
#include <iostream>


#include "Extract.h"

unsigned int READLEN = 100;

Extractor::Extractor(char* filepath) :
        ifs(filepath)
    ,   wrongCounter(0)
{

}

void Extractor::parseLines()
{

    std::string line;
    int prevChr = -1;
    // used as queue
    // front is oldest item (with lowest pos)
    std::list<ReadCount> fwdCounts;
    std::list<ReadCount> revCounts;
    while (std::getline(ifs, line))
    {

        size_t linePos = 0;
        size_t newPos = line.find_first_of('\t', linePos);
        std::string buffer(line, linePos, newPos - linePos);
        linePos = newPos + 1;
        int chrom = std::stoi(buffer);
        if (chrom == 21)
        {
            prevChr = 21;
            continue;
        }
        if (chrom != prevChr)
        {
            // clear queues
            wrongCounter += fwdCounts.size();
            fwdCounts.clear();
            wrongCounter += revCounts.size();
            revCounts.clear();
            prevChr = chrom;
        }
        // get the new position
        newPos = line.find_first_of('\t', linePos);
        buffer = std::string(line, linePos, newPos - linePos);
        linePos = newPos + 1;
        unsigned long pos = std::stoul(buffer);

        // get the 4 counts
        newPos = line.find_first_of('\t', linePos);
        buffer = std::string(line, linePos, newPos - linePos);
        linePos = newPos + 1;
        int methFwdCount = std::stoi(buffer);

        newPos = line.find_first_of('\t', linePos);
        buffer = std::string(line, linePos, newPos - linePos);
        linePos = newPos + 1;
        int unmethFwdCount = std::stoi(buffer);
        parseCountType(fwdCounts, pos, unmethFwdCount + methFwdCount);

        newPos = line.find_first_of('\t', linePos);
        buffer = std::string(line, linePos, newPos - linePos);
        linePos = newPos + 1;
        int methRevCount = std::stoi(buffer);

        buffer = std::string(line, linePos, std::string::npos);
        int unmethRevCount = std::stoi(buffer);
        parseCountType(revCounts, pos, unmethRevCount + methRevCount);


    }

}



void Extractor::parseCountType(std::list<ReadCount>& l, unsigned long pos, int count)
{
    // test for reads in queue if they overlap
    for (auto lIt = l.begin(); lIt != l.end(); )
    {
        // if read does not overlap, throw it out and count one up
        if (pos - lIt->start > READLEN)
        {
            ++lIt;
            l.pop_front();
            ++wrongCounter;

        } else {

            ++lIt;
        }
    }
    // check if we have less reads than count, if so put in new reads
    if (count > l.size())
    {
        for (unsigned int i = 0; i < count - l.size(); ++i)
        {
            l.push_back({pos});
        }
    }
}
