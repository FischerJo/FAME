//	Metal - A fast methylation alignment and calling tool for WGBS data.
//	Copyright (C) 2017  Jonas Fischer
//
//	This program is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//
//	This program is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//	Jonas Fischer	jonaspost@web.de

#ifndef READQUEUE_H
#define READQUEUE_H

#include <string>
#include <fstream>
#include <array>
#include <algorithm> // reverse, sort
#include <numeric> // iota

#ifdef _OPENMP
#include <omp.h>
#endif

// #include <sparsehash/dense_hash_map>
#include <hopscotch_map.h>

#include "gzstream/gzstream.h"
#include "CONST.h"
#include "RefGenome.h"
#include "Read.h"
#include "ShiftAnd.h"
#include "LevenshtDP.h"


class ReadQueue
{

    public:

        // Ctor -----


        ReadQueue() = delete;

        // ARGUMENTS:
        //          filePath    path to file containing reads in fastq format
        //          ref         internal representation of reference genome
        //          isGZ        flag - true iff file is gzipped
		//          bsFlag		flag - true iff there is no orientation of the read (i.e. could be C->T or G->A converted)
        ReadQueue(const char* filePath, RefGenome& ref, const bool isGZ, const bool bsFlag);

        // for paired end
        // ARGUMENTS:
        //          filePath    path to file containing read1 of paired reads in fastq format
        //          filePath2   path to file containing read2 of paired reads in fastq format
        //          ref         internal representation of reference genome
        //          isGZ        flag - true iff file is gzipped
		//          bsFlag		flag - true iff there is no orientation of read 1 of the read pair
        //
        // NOTE:
        //          provided files are ASSUMED to have equal number of reads in correct (paired) order!
        ReadQueue(const char* filePath, const char* filePath2, RefGenome& reference, const bool isGZ, const bool bsFlag);
		// for single cell paired end
		// ARGUMENTS:
		// 			...
		// 			scOutputPath	path to file which will contain single cell counts after processing
		// 			isP				flag if reads are paired
		ReadQueue(const char* scOutputPath, RefGenome& reference, const bool isGZ, const bool bsFlag, const bool isP);

        // -----------

        // Parses a chunk of the ifstream file/ igzstream file, respectively
        // Reads up to MyConst::CHUNKSIZE many reads and saves them
        // ARGUMENT:
        //          procReads   will contain number of reads that have been read into buffer
        // returns true if neither read error nor EOF occured, false otherwise
        bool parseChunk(unsigned int& procReads);
        bool parseChunkGZ(unsigned int& procReads);

		// Decides to which strand r1 should always be matched against
		void decideStrand();

        // Match all reads in readBuffer to reference genome
        // ARGUMENT:
        //          procReads   number of reads to match
		//          *Match		counter for matching statistics
		//          getStranded	flag for initial stranding counters
        // Make full match using bit shift trick
        // Align match using banded levenshtein alignment to update methylation counts
        bool matchReads(const unsigned int& procReads, uint64_t& succMatch, uint64_t& nonUniqueMatch, uint64_t& unSuccMatch, const bool getStranded);
        bool matchPairedReads(const unsigned int& procReads, uint64_t& succMatch, uint64_t& nonUniqueMatch, uint64_t& unSuccMatch, uint64_t& succPairedMatch, uint64_t& tooShortCountMatch, const bool getStranded);

		// Match batch of single cell reads given by file
		bool matchSCBatch(const char* scFile, const std::string scId, const bool isGZ);
		bool matchSCBatchPaired(const char* scFile1, const char* scFile2, const std::string scId, const bool isGZ);

        // Print the CpG methylation levels to the given filename
        // Two files are generated, one called filename_cpg.tsv
        // the other one called filename_cpg_start.tsv
        // The latter file contaings the CpGs very close to the border of the chromosomes.
        // Such CpGs are rarely found in large genome assemblys
        //
        // Both contain methylation counts for each CpG in the reference genome
        // using the following (tab separated) format:
        //
        // Chromosome   Position    Fwd_methylated  Fwd_unmethylated    Rev_methylated  Rev_unmethylated
        //
        //, where
        // Chromosome is the chromosome ID
        //
        // Position is the offset of the (C on the forward strand of the) CpG inside the
        // chromosome, measured from the start of the chromosomal sequence and is zero based
        //
        // Fwd_methylated is the number of aligned methylated CpGs on the forward strand CpG
        // Fwd_unmethylated  -------" " ---------- unmeethylated --------- " " ------------
        // Rev_methylated is the number of aligned methylated CpGs on the reverse strand CpG
        // Rev_unmethylated  -------" " ---------- unmeethylated --------- " " ------------
        //
        // ARGUMENT:
        //          filename    desired basename for the output files
        void printMethylationLevels(std::string& filename);

		void printSCMethylationLevels(const std::string scID);


    private:

        // hash function for std::unordered_set of meta ids
        struct MetaHash
        {
                // just use meta IDs as they are unique
               size_t operator() (const uint32_t& metaID) const {

                   return static_cast<size_t>(metaID);
               }
        };

		// query the k-mers in internal data structures (countsFwdStart/countsRevStart, paired_counts...) to given modified shift-and automaton
		// only queries to MetaCpGs with enough k-mers (# >= qThreshold)
		//
		// ARGUMENTS:
		// 			sa			Shift-And automaton for reads sequence
		// 			mat			empty struct/DS that will hold best match/matches found
		// 			qThreshold	minimum number of k-mers needed for Meta CpG to be queried to shift-and
		//
		// MODIFICATION:
		// 			Adds the best found match/ matches to mat/mats
		inline int saQuerySeedSetRef(ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa, MATCH::match& mat, uint16_t& qThreshold);
		inline void saQuerySeedSetRefFirst(ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa, std::vector<MATCH::match>& mats, const uint16_t& qThreshold);
		inline void saQuerySeedSetRefSecond(ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa, std::vector<MATCH::match>& mats, const uint16_t& qThreshold);

        // count all metaCpG occurences of k-mers appearing in seq
        //
        // ARGUMENTS:
        //          seq             sequence of the read to query to hash table
        //
        // MODIFICATION:
        //          The threadCount* fields are modified such that they have the count of metaCpGs after
        //          a call to this function
		inline void getSeedRefs(const std::string& seq, const size_t& readSize, const uint16_t qThreshold);
		inline uint16_t getSeedRefsFirstRead(const std::string& seq, const size_t& readSize, const uint16_t qThreshold);
		inline bool getSeedRefsSecondRead(const std::string& seq, const size_t& readSize, const uint16_t qThreshold);

		inline bool matchFwdFirst(std::pair<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool> >& meta, uint8_t& prevChr, uint64_t& prevOff, std::vector<MATCH::match>& mats, int32_t& bmCount, uint16_t qThreshold, ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa);
		inline bool matchRevFirst(std::pair<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool> >& meta, uint8_t& prevChr, uint64_t& prevOff, std::vector<MATCH::match>& mats, int32_t& bmCount, uint16_t qThreshold, ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa);
		inline bool matchFwdSecond(std::pair<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool> >& meta, uint8_t& prevChr, uint64_t& prevOff, std::vector<MATCH::match>& mats, int32_t& bmCount, uint16_t qThreshold, ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa);
		inline bool matchRevSecond(std::pair<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool> >& meta, uint8_t& prevChr, uint64_t& prevOff, std::vector<MATCH::match>& mats, int32_t& bmCount, uint16_t qThreshold, ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa);

		inline void sort_by_n(std::vector<unsigned int>::iterator it_start, std::vector<unsigned int>::iterator it_n, std::vector<unsigned int>::iterator it_end, std::vector<uint64_t>& sliceOff, std::vector<bool>& sliceIsDone)
		{
			// order descending on window and reverse strand > fwd strand
			// Order:
			//		return false (i.e. id1 > id2) if
			//			has no windows left to process in 1
			//		return true (i.e. id1 < id2) if
			//			window id 1 > window id 2 and
			//				strand is rev for 1 and is fwd for 2
			//			else return false
			std::nth_element(it_start, it_n, it_end,
					[&](unsigned int id1, unsigned int id2){
						if (sliceIsDone[id1])
						{
							return false;
						} else if (sliceIsDone[id2])
						{
							return true;
						}

						const auto id1meta = KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff[id1]]);
						const auto id2meta = KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff[id2]]);
						if (id1meta > id2meta)
						{
							return true;

						} else {
							// test for strand if windows are equal
							if (id1meta == id2meta &&
									(ref.strandTable[sliceOff[id1]] > ref.strandTable[sliceOff[id2]]))
							{
								// Case id2 window is reverse strand, id1 window is fwd strand
								return true;
							} else {
								return false;
							}
						}
					});
		}

        // Extract single match for given lists of matches of fwd and reverse complement of a single read
        // Internally updates methylation counts
        //
        // ARGUMENTS:
        //          fwdMatches  array of matches retrieved for original sequence
        //          revMatches  array of matches retrieved for reverse complement of sequence
        //          r           read representation
        //          revSeq      reverse complement sequence of read
        //
        // RETURN:  true iff successfully extracted match
        //
        inline bool extractSingleMatch(std::vector<MATCH::match>& fwdMatches, std::vector<MATCH::match>& revMatches, Read& r, std::string& revSeq);

        // Examines if two matchings build a pair
        //
        // ARGUMENTS:
        //          mat1    matching of first read
        //          mat2    matching of second read
        //
        // RETURN:
        //          -1      iff no pairing
        //          n       iff pairing, where n is the number of errors summed over both matchings
        //
        inline int extractPairedMatch(MATCH::match& mat1, MATCH::match& mat2);


        // compute the methylation levels for the given read by traversing the CpGs of the matched meta CpG
        //
        // ARGUMENTS:
        //          mat     match to process
        //          seq     sequence that was matched (i.e. r.seq or revSeq in main query routine)
        //
        // MODIFICATIONS:
        //              will modify internal methLevel counters
        inline void computeMethLvl(MATCH::match& mat, std::string& seq);

        // input stream of file given as path to Ctor
        std::ifstream file;
        igzstream igz;
        // second file if paired
        std::ifstream file2;
        igzstream igz2;


        // representation of the reference genome
        RefGenome& ref;

        // buffer holding MyConst::CHUNKSIZE many reads
        std::vector<Read> readBuffer;
        // second buffer for paired reads
        std::vector<Read> readBuffer2;

        // mapping of letters to array indices for shift and algorithm
        // 'A' -> 0
        // 'C' -> 1
        // 'G' -> 2
        // 'T' -> 3
        std::array<uint8_t, 16> lmap;

        std::array<tsl::hopscotch_map<uint32_t, uint16_t, MetaHash>, CORENUM> fwdMetaIDs;
        std::array<tsl::hopscotch_map<uint32_t, uint16_t, MetaHash>, CORENUM> revMetaIDs;
        // Holds counts for each thread for counting heuristic
        // KEY: Meta CpG ID
        // VALUE:
        //      1) K-mer count of first read
        //      2) K-mer count of second read
		//      3) Boolean flag that is true iff second read has enough kmers in this or adjacent MetaCpGs
        //      4) Boolean flag that is true iff first read is matched to this or adjacent Meta CpGs
        //
        // std::array<google::dense_hash_map<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>, MetaHash>, CORENUM> paired_fwdMetaIDs;
        // std::array<google::dense_hash_map<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>, MetaHash>, CORENUM> paired_revMetaIDs;
        std::array<tsl::hopscotch_map<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>, MetaHash>, CORENUM> paired_fwdMetaIDs;
        std::array<tsl::hopscotch_map<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>, MetaHash>, CORENUM> paired_revMetaIDs;

        bool isPaired;
		bool isSC;

        // comparison for match
        struct CompiFwd {

            inline uint16_t operator() (const char& cRead, const char& cRef)
            {
                // bisulfite antisymmetry
                if (cRead == 'T')
                {
                    if (cRef == 'C')
                        return 0;
                }
                return cRead == cRef ? 0 : 1;
            }

        } cmpFwd;
        // reverse complement matching using the original strand
        struct CompiRev {

            inline uint16_t operator() (const char& cRead, const char& cRef)
            {
                switch (cRead)
                {
                    case ('T') :
                        // bisulfite antisymmetry
                        if (cRef == 'A' || cRef == 'G')
                            return 0;
                        break;
                    case ('G') :
                        if (cRef == 'C')
                            return 0;
                        break;
                    case ('A'):
                        if (cRef == 'T')
                            return 0;
                        break;
                    case ('C'):
                        if (cRef == 'G')
                            return 0;
                        break;
                }

                return 1;
            }

        } cmpRev;

		struct CompiMetaFirst {
			bool operator()(std::pair<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>>& a, std::pair<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>>& b) const
			{
				return std::get<0>(a.second) > std::get<0>(b.second);
			}
		} cmpMetaFirst;
		struct CompiMetaSecond {
			bool operator()(std::pair<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>>& a, std::pair<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>>& b) const
			{
				return std::get<1>(a.second) > std::get<1>(b.second);
			}
		} cmpMetaSecond;

        // holding the counts for a single CpG
        struct methLvl
        {
            // for forward strand
            uint16_t methFwd;
            uint16_t unmethFwd;

            // for reverse strand
            uint16_t methRev;
            uint16_t unmethRev;
        };
        // holds the counts for each CpG
        // indexed by the same indices as of the cpgTable vector in RefGenome class
        std::vector<struct methLvl> methLevels;
        std::vector<struct methLvl> methLevelsStart;

		// Information for single cell data
		// Holds for the most recent single cell the individual counts
		std::vector<struct methLvl> methLevelsSc;
		// Mapping of internal ids to external identifier tags
		std::unordered_map<size_t, std::string> scIds;
		// output file for sc data
		std::ofstream scOutput;

        // holds matching info
        std::array<uint64_t, CORENUM> matchStats;
        std::array<uint64_t, CORENUM> nonUniqueStats;
        std::array<uint64_t, CORENUM> noMatchStats;
        std::array<uint64_t, CORENUM> matchPairedStats;
        std::array<uint64_t, CORENUM> tooShortCounts;

		bool bothStrandsFlag;
		// counter for read 1 matches to fwd strand
		uint64_t r1FwdMatches;
		// counter for read 1 matches to rev strand
		uint64_t r1RevMatches;
		bool matchR1Fwd;

        // TODO
        std::ofstream of;
        inline void printMatch(std::ostream& o, MATCH::match& mat)
        {
            o << "Match at offset " << MATCH::getOffset(mat) + ref.metaWindows[MATCH::getMetaID(mat)].startPos << \
                " on chromosome " << static_cast<uint64_t>(ref.metaWindows[MATCH::getMetaID(mat)].chrom) << \
                " on strand " << (MATCH::isFwd(mat) ? "fwd" : "rev") << " with " << MATCH::getErrNum(mat) << " many errors. Meta CpG " << MATCH::getMetaID(mat);
        }


};

#endif /* READQUEUE_H */
