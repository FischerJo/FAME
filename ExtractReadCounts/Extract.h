#ifndef EXTRACT_H
#define EXTRACT_H

#include <string>
#include <fstream>
#include <list>


class Extractor
{
    public:

        Extractor(char* filepath);

        inline unsigned long getCount()
        {
            return wrongCounter;
        }

        void parseLines();

    private:

        struct ReadCount
        {
            // (upper bound on) start position
            unsigned long start;

            // chromosome
            // unsigned int chr;
        };
        void parseCountType(std::list<ReadCount>& l, unsigned long pos, int count);

        std::ifstream ifs;

        // counts how many reads are mapped to the wrong chromosome (not 21)
        unsigned long wrongCounter;

};

#endif /* EXTRACT_H */
