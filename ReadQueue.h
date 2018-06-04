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

#include <sparsehash/dense_hash_map>

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
		inline void saQuerySeedSetRefFirst(ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa, std::vector<MATCH::match>& mats, uint16_t& qThreshold);
		inline void saQuerySeedSetRefSecond(ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa, std::vector<MATCH::match>& mats, uint16_t& qThreshold);

        // count all metaCpG occurences of k-mers appearing in seq
        //
        // ARGUMENTS:
        //          seq             sequence of the read to query to hash table
        //
        // MODIFICATION:
        //          The threadCount* fields are modified such that they have the count of metaCpGs after
        //          a call to this function
		inline void getSeedRefs(const std::string& seq, const size_t& readSize, const uint16_t qThreshold);
        // TODO start Metas
		inline void getSeedRefsFirstRead(const std::string& seq, const size_t& readSize, const uint16_t qThreshold);
        // TODO start Metas
		inline void getSeedRefsSecondRead(const std::string& seq, const size_t& readSize, const uint16_t qThreshold);


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

		// Karl's matching
		//
		// ARGUMENTS:
		//			seq			read sequence to match
		//			qThreshold	minimum number of k-mers required in Meta CpG to test for matching
		//			sa			ShiftAnd automaton
		//			mat			Best match found (or dummy if none)
		//			threadnum	Number of the thread in which this function is called
		//
		//	RETURN:
		//			-1			iff no successfull match (e.g. nonunique)
		//			0			no match at all
		//			1			successfull match
		//
		//	MODIFICATION:
		//			Adapts qThreshold if match is found.
		//
		inline int matchSingle(const std::string& seq, uint16_t& qThreshold, ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa, MATCH::match& mat, const int threadnum)
		{

			// std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
			const size_t kmerNum = seq.size() - MyConst::KMERLEN + 1;

			// slices of hash table currently looking at
			std::vector<uint64_t>& sliceOff = sliceOffThreads[threadnum];
			std::vector<uint64_t>& sliceEnd = sliceEndThreads[threadnum];
			// referencing indices correspond to indices of sliceHashes
			std::vector<unsigned int>& sliceSortedIds = sliceSortedIdsThreads[threadnum];
			sliceSortedIds.resize(kmerNum);
			// flags if whole slice is already processed
			std::vector<bool>& sliceIsDone = sliceIsDoneThreads[threadnum];
			std::fill(sliceIsDone.begin(), sliceIsDone.end(), false);


			std::iota(sliceSortedIds.begin(), sliceSortedIds.end(), 0);

            // retrieve start and end of reference kmer list for initial read kmer
			// uint64_t fhVal = ntHash::NTP64(seq.data());
			// TODO: spaced
			uint64_t fhVal;
			// std::chrono::high_resolution_clock::time_point startTime2 = std::chrono::high_resolution_clock::now();
			uint64_t sfVal = ntHash::NTPS64(seq.data(), MyConst::SEED, MyConst::KMERLEN, fhVal);
			// std::chrono::high_resolution_clock::time_point endTime2 = std::chrono::high_resolution_clock::now();
			// auto runtime2 = std::chrono::duration_cast<std::chrono::microseconds>(endTime2 - startTime2).count();

			sliceOff[0] = ref.tabIndex[sfVal % MyConst::HTABSIZE];
			sliceEnd[0] = ref.tabIndex[(sfVal % MyConst::HTABSIZE) + 1];

			for (size_t i = 1; i < kmerNum; ++i)
			{
				// ntHash::NTP64(fhVal, seq[i-1], seq[i-1+MyConst::KMERLEN]);
				// TODO: spaced
				// std::chrono::high_resolution_clock::time_point startTime3 = std::chrono::high_resolution_clock::now();
				sfVal = ntHash::NTPS64(seq.data()+i, MyConst::SEED, seq[i-1], seq[i-1+MyConst::KMERLEN], MyConst::KMERLEN, fhVal);
				// std::chrono::high_resolution_clock::time_point endTime3 = std::chrono::high_resolution_clock::now();
				// auto runtime3 = std::chrono::duration_cast<std::chrono::microseconds>(endTime3 - startTime3).count();
				// of << "Rolling hash " << i << " runtime: " << runtime3 << "\n";
				sliceOff[i] = ref.tabIndex[sfVal % MyConst::HTABSIZE];
				sliceEnd[i] = ref.tabIndex[(sfVal % MyConst::HTABSIZE) + 1];
			}
			for (size_t i = 0; i < kmerNum; ++i)
			{
				if (sliceOff[i] >= sliceEnd[i])
				{
					sliceIsDone[i] = true;
					// of << "Happening! End for " << i << "\n";
				}

			}
			// std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
			// auto runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
            //
			// of <<"Runtime init: " << runtime << "\n";
			// of <<"Initital hash runtime: " << runtime2 <<"\n";
            //
			// // TODO
			// for (size_t j = 0; j < kmerNum; ++j)
			// {
			// 	of << "Cell " << j << ": " << sliceEnd[j] - sliceOff[j] << "\n";
			// 	// for (unsigned int l = sliceOff[j]; l < sliceEnd[j]; ++l)
			// 	// {
			// 	// 	of << KMER_S::getMetaCpG(ref.kmerTableSmall[l]) << " / " << ref.cpgTable[ref.metaCpGs[ref.kmerTableSmall[l]].start].pos << "\n";
			// 	// }
			// 	of << "\n---\n\n";
			// }
			// counter for how often we had a match
			std::array<uint8_t, MyConst::ADDMIS + MyConst::MISCOUNT + 1> multiMatch;
			multiMatch.fill(0);

			// will contain matches iff match is found for number of errors specified by index
			std::array<MATCH::match, MyConst::ADDMIS + MyConst::MISCOUNT + 1> uniqueMatches;
			// store the last match found in current MetaCpG
			uint8_t prevChr = 0;
			uint64_t prevOff = 0xffffffffffffffffULL;

			sort_by_n(sliceSortedIds.begin(), sliceSortedIds.begin() + qThreshold - 1, sliceSortedIds.end(), sliceOff, sliceIsDone);
			// of << "SortedIds:\n\t";
			// for (size_t i = 0; i < kmerNum; ++i)
			// {
			// 	of << sliceSortedIds[i] << "\t";
			// }
			// of << "\n";

			while (!sliceIsDone[sliceSortedIds[qThreshold-1]])
			{

				bool unchanged = true;
				uint32_t qWindow = KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff[sliceSortedIds[qThreshold-1]]]);
				//TODO
				// of << "qWindow ID: " << qWindow << "\t";
				// of << ref.cpgTable[ref.metaCpGs[qWindow].start].pos << "\n";
                //
				// startTime = std::chrono::high_resolution_clock::now();
				for (size_t i = 0; i < qThreshold-1; ++i)
				{
					// of << sliceSortedIds[i] << "\t";
					// advance pointers while window id is larger
					while (KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff[sliceSortedIds[i]]]) > qWindow)
					{
						//TODO
						// of << KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff[sliceSortedIds[i]]]) << "++";
						++sliceOff[sliceSortedIds[i]];
						unchanged = false;
						// reached end of vector slice
						if (sliceOff[sliceSortedIds[i]] >= sliceEnd[sliceSortedIds[i]])
						{
							// TODO
							// of << "Setting " << sliceSortedIds[i] << " to End\n";
							sliceIsDone[sliceSortedIds[i]] = true;
							break;
						}
					}
					// of << "\n";

					// advance pointer if strand is fwd but qStrand is rev
					if (KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff[sliceSortedIds[i]]]) == qWindow &&
							!sliceIsDone[sliceSortedIds[i]] &&
							ref.strandTable[sliceOff[sliceSortedIds[i]]] > ref.strandTable[sliceOff[sliceSortedIds[qThreshold-1]]])
					{
						//TODO
						// of << "WrongStrand++" << KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff[sliceSortedIds[i]]]) << "\n";
						++sliceOff[sliceSortedIds[i]];
						unchanged = false;
						// reached end of vector slice
						if (sliceOff[sliceSortedIds[i]] >= sliceEnd[sliceSortedIds[i]])
						{
							//TODO
							// of << "Setting " << sliceSortedIds[i] << " to End\n";
							sliceIsDone[sliceSortedIds[i]] = true;
						}
					}
				}

				// endTime = std::chrono::high_resolution_clock::now();
				// runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
				// of << "\nAdvance runtime " << runtime << "\t\t";
				// if nothing has changed and at least one k-mer has windows left to process, test for match
				if (unchanged)
				{
					// startTime = std::chrono::high_resolution_clock::now();
					const uint32_t metaID = KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff[sliceSortedIds[0]]]);
					// TODO
					// of << "\nMatching Meta:\n\t" << metaID << "\t\t <--- \n";
					const bool metaStrand = ref.strandTable[sliceOff[sliceSortedIds[0]]];

					const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[metaID].start];
					const struct CpG& endCpg = ref.cpgTable[ref.metaCpGs[metaID].end];
					auto startIt = ref.fullSeq[startCpg.chrom].begin() + startCpg.pos;
					auto endIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2) + MyConst::MISCOUNT + MyConst::ADDMIS;

					// check if CpG was too near to the end
					if (endIt > ref.fullSeq[startCpg.chrom].end())
					{
						// if so move end iterator appropriately
						endIt = ref.fullSeq[startCpg.chrom].end();
					}

					// use shift and to find all matchings
					std::vector<uint64_t> matchings;
					std::vector<uint8_t> errors;
					// endTime = std::chrono::high_resolution_clock::now();
					// runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
					// of << "\nMatch preprocessing runtime " << runtime << "\t\t";
					// startTime = std::chrono::high_resolution_clock::now();
					if (metaStrand)
					{
						sa.querySeq(startIt, endIt, matchings, errors);
					} else {
						--startIt, --endIt;
						sa.queryRevSeq(endIt, startIt, matchings, errors);
					}
					// endTime = std::chrono::high_resolution_clock::now();
					// runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
					// of << "\nShiftAnd runtime " << runtime << "\t\t";
                    //
					// startTime = std::chrono::high_resolution_clock::now();

					size_t i = 0;
					// compare first found match with last found match of previous meta CpG
					if (matchings.size() > 0)
					{
						// compare chromosome and offset
						if (matchings[0] + ref.cpgTable[ref.metaCpGs[metaID].start].pos == prevOff && ref.cpgTable[ref.metaCpGs[metaID].start].chrom == prevChr)
						{
							++i;
						}
					}
					// go through matching and see if we had such a match (with that many errors) before - if so,
					// return to caller reporting no match
					for (; i < matchings.size(); ++i)
					{

						// check if we had a match with that many errors before
						if (multiMatch[errors[i]])
						{

							MATCH::match& match_2 = uniqueMatches[errors[i]];
							// const bool isStart = MATCH::isStart(match_2);
							const bool isFwd = MATCH::isFwd(match_2);
							// check if same k-mer (borders of meta CpGs)
							if (ref.cpgTable[ref.metaCpGs[MATCH::getMetaID(match_2)].start].pos + MATCH::getOffset(match_2) == startCpg.pos + matchings[i])
							{
								if ((isFwd && metaStrand) || (!isFwd && !metaStrand))
									continue;

							} else {

								// check if this is a match without errors
								if (!errors[i])
								{

									// if so, return without a match
									// of << "\tNonunique <---\n\n";
									return -1;

								}
								// set the number of matches with that many errors to 2
								// indicating that we do not have a unique match with that many errors
								multiMatch[errors[i]] = 2;
							}


						} else {

							// update qgram lemma
							uint16_t newQ = sa.size() - MyConst::KMERLEN - (MyConst::KMERLEN * errors[i]);
							// check for overflow and if we improved old q
							if (newQ < sa.size() && newQ > qThreshold)
								qThreshold = newQ;


							// we don't have such a match yet,
							// so save this match at the correct position
							if (metaStrand)
							{
								uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i], errors[i], 1, 0, metaID);
							} else {
								uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i], errors[i], 0, 0, metaID);
							}
							multiMatch[errors[i]] = 1;
						}
					}
					if (matchings.size() > 0)
					{

						prevChr = ref.cpgTable[ref.metaCpGs[metaID].start].chrom;
						prevOff = ref.cpgTable[ref.metaCpGs[metaID].start].pos + matchings[matchings.size() - 1];

					} else {

						prevChr = 0;
						prevOff = 0xffffffffffffffffULL;
					}
					// endTime = std::chrono::high_resolution_clock::now();
					// runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
					// of << "\nPostprocessing matches runtime " << runtime << "\t\t";
					// startTime = std::chrono::high_resolution_clock::now();
					// Advance all offsets of k-mers with matched window
					const bool qStrand = ref.strandTable[sliceOff[sliceSortedIds[qThreshold-1]]];
					for (size_t i = 0; i < qThreshold; ++i)
					{
						// of << sliceSortedIds[i] << "\t";
						// of << KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff[sliceSortedIds[i]]]) << "++";
						++sliceOff[sliceSortedIds[i]];
						// reached end of vector slice
						if (sliceOff[sliceSortedIds[i]] >= sliceEnd[sliceSortedIds[i]])
						{
							// of << "Setting " << sliceSortedIds[i] << " to End\n";
							sliceIsDone[sliceSortedIds[i]] = true;
						}
					}
					for (size_t i = qThreshold; i < kmerNum; ++i)
					{
						if (sliceIsDone[sliceSortedIds[i]])
							continue;
						if (KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff[sliceSortedIds[i]]]) == qWindow)
						{
							if (ref.strandTable[sliceOff[sliceSortedIds[i]]] == qStrand)
							{
								// of << sliceSortedIds[i] << "\t";
								// of << KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff[sliceSortedIds[i]]]) << "++";
								++sliceOff[sliceSortedIds[i]];
								// reached end of vector slice
								if (sliceOff[sliceSortedIds[i]] >= sliceEnd[sliceSortedIds[i]])
								{
									// of << "Setting " << sliceSortedIds[i] << " to End\n";
									sliceIsDone[sliceSortedIds[i]] = true;
								}
							}
						}
					}
					// endTime = std::chrono::high_resolution_clock::now();
					// runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
					// of << "Advance after match runtime " << runtime << "\t\t";
				}

				// startTime = std::chrono::high_resolution_clock::now();
				sort_by_n(sliceSortedIds.begin(), sliceSortedIds.begin() + qThreshold - 1, sliceSortedIds.end(), sliceOff, sliceIsDone);
				// endTime = std::chrono::high_resolution_clock::now();
				// runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
				// of << "Sort runtime " << runtime << "\n\n";
				// of << "SortedIds:\n\t";
				// for (size_t i = 0; i < kmerNum; ++i)
				// {
				// 	of << sliceSortedIds[i] << "\t";
				// }
				// of << "\n";

			}

			// of << "qGrams at end: \n";
			// for (size_t i = 0; i < kmerNum; ++i)
			// {
			// 	if (!sliceIsDone[sliceSortedIds[qThreshold-1]])
			// 		of << "\t" << KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff[sliceSortedIds[i]]]) << "\n";
			// 	else
			// 		of << "\tEnd\n";
			// }

            // go through found matches for each [0,maxErrorNumber] and see if it is unique
            for (size_t i = 0; i < multiMatch.size(); ++i)
            {
                // there is no match with that few errors, search the one with more errors
                if (multiMatch[i] == 0)
                {
                    continue;
                }
                mat = uniqueMatches[i];
                // if match is not unique, return unsuccessfull to caller
                if (multiMatch[i] > 1)
                {

					// of << "\tNonunique <---\n\n";
                    return -1;

                // exactly one with that many errors - return successfull
                } else {

					// of << "\tSuccessfull match <---\n\n";
                    return 1;
                }

            }
            // we have not a single match at all, return unsuccessfull to caller
			// of << "\tNo match <---\n\n";
            return 0;
		}
		// Karl's matching
		//
		// ARGUMENTS:
		//			seq1		read 1 sequence to match
		//			seq2		read 2 sequence to match
		//			qThreshold	minimum number of k-mers required in Meta CpG to test for matching
		//			sa1			ShiftAnd automaton for read 1
		//			sa2			ShiftAnd automaton for read 2
		//			mat			pair of matching positions for read1 and read2, if one is found
		//
		//	RETURN:
		//			-1			iff no successfull match (e.g. nonunique)
		//			0			no match at all
		//			1			successfull match
		//
		//	MODIFICATION:
		//			Adapts qThreshold if match is found.
		//
		// inline int matchPaired(const std::string& seq1, const std::string& seq2,uint16_t& qThreshold, ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa1, ShiftAnd<MyConst::MISCOUNT + MyConst::ADDMIS>& sa2, std::pair<MATCH::match, MATCH::match> mat)
		// {
        //
		// 	// std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
		// 	const size_t kmerNum1 = seq1.size() - MyConst::KMERLEN + 1;
		// 	const size_t kmerNum2 = seq2.size() - MyConst::KMERLEN + 1;
		// 	// const std::array<const size_t> kmerNum = {{kmerNum1, kmerNum2}};
        //
		// 	// slices of hash table currently looking at
		// 	std::vector<uint64_t> sliceOff1 (kmerNum1);
		// 	std::vector<uint64_t> sliceEnd1 (kmerNum1);
		// 	std::vector<uint64_t> sliceOff2 (kmerNum2);
		// 	std::vector<uint64_t> sliceEnd2 (kmerNum2);
		// 	// const std::array<std::vector<uint64_t>*, 2> sliceOff = {{&sliceOff1, &sliceOff2}};
		// 	// const std::array<std::vector<uint64_t>*, 2> sliceEnd = {{&sliceEnd1, &sliceEnd2}};
		// 	// referencing indices correspond to indices of sliceHashes
		// 	std::vector<unsigned int> sliceSortedIds1 (kmerNum1);
		// 	std::vector<unsigned int> sliceSortedIds2 (kmerNum2);
		// 	// const std::array<std::vector<uint64_t>*, 2> sliceSortedIds = {{&sliceSortedIds1, &sliceSortedIds2}};
		// 	// flags if whole slice is already processed
		// 	std::vector<bool> sliceIsDone1 (kmerNum1, false);
		// 	std::vector<bool> sliceIsDone2 (kmerNum2, false);
		// 	// const std::array<std::vector<uint64_t>*, 2> sliceIsDone = {{&sliceIsDone1, &sliceIsDone2}};
        //
        //
		// 	std::iota(sliceSortedIds1.begin(), sliceSortedIds1.end(), 0);
        //     // retrieve start and end of reference kmer list for initial read kmer of read 1
		// 	// uint64_t fhVal = ntHash::NTP64(seq.data());
		// 	// TODO: spaced
		// 	uint64_t fhVal;
		// 	uint64_t sfVal = ntHash::NTPS64(seq1.data(), MyConst::SEED, MyConst::KMERLEN, fhVal);
        //
		// 	sliceOff1[0] = ref.tabIndex[sfVal % MyConst::HTABSIZE];
		// 	sliceEnd1[0] = ref.tabIndex[(sfVal % MyConst::HTABSIZE) + 1];
        //
		// 	for (size_t i = 1; i < kmerNum1; ++i)
		// 	{
		// 		// ntHash::NTP64(fhVal, seq[i-1], seq[i-1+MyConst::KMERLEN]);
		// 		// TODO: spaced
		// 		sfVal = ntHash::NTPS64(seq1.data()+i, MyConst::SEED, seq1[i-1], seq1[i-1+MyConst::KMERLEN], MyConst::KMERLEN, fhVal);
		// 		sliceOff1[i] = ref.tabIndex[sfVal % MyConst::HTABSIZE];
		// 		sliceEnd1[i] = ref.tabIndex[(sfVal % MyConst::HTABSIZE) + 1];
		// 	}
		// 	for (size_t i = 0; i < kmerNum1; ++i)
		// 	{
		// 		if (sliceOff1[i] >= sliceEnd1[i])
		// 		{
		// 			sliceIsDone1[i] = true;
		// 			// of << "Happening! End for " << i << "\n";
		// 		}
        //
		// 	}
		// 	std::iota(sliceSortedIds2.begin(), sliceSortedIds2.end(), 0);
        //     // retrieve start and end of reference kmer list for initial read kmer of read 2
		// 	// uint64_t fhVal = ntHash::NTP64(seq.data());
		// 	// TODO: spaced
		// 	sfVal = ntHash::NTPS64(seq2.data(), MyConst::SEED, MyConst::KMERLEN, fhVal);
		// 	sliceOff2[0] = ref.tabIndex[sfVal % MyConst::HTABSIZE];
		// 	sliceEnd2[0] = ref.tabIndex[(sfVal % MyConst::HTABSIZE) + 1];
        //
		// 	for (size_t i = 1; i < kmerNum2; ++i)
		// 	{
		// 		// ntHash::NTP64(fhVal, seq[i-1], seq[i-1+MyConst::KMERLEN]);
		// 		// TODO: spaced
		// 		sfVal = ntHash::NTPS64(seq2.data()+i, MyConst::SEED, seq2[i-1], seq2[i-1+MyConst::KMERLEN], MyConst::KMERLEN, fhVal);
		// 		sliceOff2[i] = ref.tabIndex[sfVal % MyConst::HTABSIZE];
		// 		sliceEnd2[i] = ref.tabIndex[(sfVal % MyConst::HTABSIZE) + 1];
		// 	}
		// 	for (size_t i = 0; i < kmerNum2; ++i)
		// 	{
		// 		if (sliceOff2[i] >= sliceEnd2[i])
		// 		{
		// 			sliceIsDone2[i] = true;
		// 			// of << "Happening! End for " << i << "\n";
		// 		}
        //
		// 	}
            //
			// // std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
			// // auto runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
            // //
			// // of <<"Runtime init: " << runtime << "\n";
            //
			// // TODO
			// // for (size_t j = 0; j < kmerNum; ++j)
			// // {
			// // 	of << "Cell " << j << ": " << sliceEnd[j] - sliceOff[j] << "\n";
			// // 	// for (unsigned int l = sliceOff[j]; l < sliceEnd[j]; ++l)
			// // 	// {
			// // 	// 	of << KMER_S::getMetaCpG(ref.kmerTableSmall[l]) << " / " << ref.cpgTable[ref.metaCpGs[ref.kmerTableSmall[l]].start].pos << "\n";
			// // 	// }
			// // 	of << "\n---\n\n";
			// // }
			// // store the last match found in current MetaCpG
			// uint8_t prevChr = 0;
			// uint64_t prevOff = 0xffffffffffffffffULL;
            //
			// // order descending on window and reverse strand > fwd strand
			// // Order:
			// //		return false (i.e. id1 > id2) if
			// //			has no windows left to process in 1
			// //		return true (i.e. id1 < id2) if
			// //			window id 1 > window id 2 and
			// //				strand is rev for 1 and is fwd for 2
			// //			else return false
			// sort_by_n(sliceSortedIds1.begin(), sliceSortedIds1.begin() + qThreshold - 1, sliceSortedIds1.end(), sliceOff1, sliceIsDone1);
			// sort_by_n(sliceSortedIds2.begin(), sliceSortedIds2.begin() + qThreshold - 1, sliceSortedIds2.end(), sliceOff2, sliceIsDone2);
			// // of << "SortedIds:\n\t";
			// // for (size_t i = 0; i < kmerNum; ++i)
			// // {
			// // 	of << sliceSortedIds[i] << "\t";
			// // }
			// // of << "\n";
			// constexpr unsigned int contextWLen = (unsigned int)(((double)MyConst::MAXPDIST / MyConst::WINLEN) + 0.5);
            //
			// while (!sliceIsDone1[sliceSortedIds1[qThreshold-1]] && !sliceIsDone2[sliceSortedIds2[qThreshold-1]])
			// {
            //
			// 	const uint32_t qWindow1 = KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff1[sliceSortedIds1[qThreshold-1]]]);
			// 	const uint32_t qWindow2 = KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff2[sliceSortedIds2[qThreshold-1]]]);
            //
			// 	bool unchanged = true;
			// 	// if window of read 1 is smaller, i.e. is further to the right in the genome
			// 	if (qWindow1 <= qWindow2)
			// 	{
			// 		// advance all lists of read 1 before pivot 1 until at least pivot 1
			// 		for (size_t i = 0; i < qThreshold-1; ++i)
			// 		{
			// 			// of << sliceSortedIds[i] << "\t";
			// 			// advance pointers while window id is larger
			// 			while (KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff1[sliceSortedIds1[i]]]) > qWindow1)
			// 			{
			// 				//TODO
			// 				// of << KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff[sliceSortedIds[i]]]) << "++";
			// 				++sliceOff1[sliceSortedIds1[i]];
			// 				unchanged = false;
			// 				// reached end of vector slice
			// 				if (sliceOff1[sliceSortedIds1[i]] >= sliceEnd1[sliceSortedIds1[i]])
			// 				{
			// 					// TODO
			// 					// of << "Setting " << sliceSortedIds[i] << " to End\n";
			// 					sliceIsDone1[sliceSortedIds1[i]] = true;
			// 					break;
			// 				}
			// 			}
			// 			// of << "\n";
            //
			// 			// advance pointer if strand is fwd but qStrand is rev
			// 			if (KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff1[sliceSortedIds1[i]]]) == qWindow1 &&
			// 					!sliceIsDone1[sliceSortedIds1[i]] &&
			// 					ref.strandTable[sliceOff1[sliceSortedIds1[i]]] > ref.strandTable[sliceOff1[sliceSortedIds1[qThreshold-1]]])
			// 			{
			// 				//TODO
			// 				// of << "WrongStrand++" << KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff[sliceSortedIds[i]]]) << "\n";
					// 		++sliceOff1[sliceSortedIds1[i]];
					// 		unchanged = false;
					// 		// reached end of vector slice
					// 		if (sliceOff1[sliceSortedIds1[i]] >= sliceEnd1[sliceSortedIds1[i]])
					// 		{
					// 			//TODO
					// 			// of << "Setting " << sliceSortedIds[i] << " to End\n";
					// 			sliceIsDone1[sliceSortedIds1[i]] = true;
					// 		}
					// 	}
					// }
                    //
					// // advance all lists of read2 until max {pivot2, pivot1 + contextWLen}
					// const unsigned int minmaxpiv = std::max(qWindow2, qWindow1+contextWLen);
					// const bool minmaxstrand = !ref.strandTable[sliceOff1[sliceSortedIds1[qThreshold-1]]];
					// for (size_t i = 0; i < kmerNum2; ++i)
					// {
					// 	// of << sliceSortedIds[i] << "\t";
					// 	// advance pointers while window id is larger
					// 	while (KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff2[sliceSortedIds2[i]]]) > minmaxpiv)
					// 	{
					// 		//TODO
					// 		// of << KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff[sliceSortedIds[i]]]) << "++";
					// 		++sliceOff2[sliceSortedIds2[i]];
					// 		// reached end of vector slice
					// 		if (sliceOff2[sliceSortedIds2[i]] >= sliceEnd2[sliceSortedIds2[i]])
					// 		{
					// 			// TODO
					// 			// of << "Setting " << sliceSortedIds[i] << " to End\n";
					// 			sliceIsDone2[sliceSortedIds2[i]] = true;
					// 			break;
					// 		}
					// 	}
					// 	// of << "\n";
                    //
					// 	// advance pointer if strand is fwd but qStrand is rev
					// 	if (KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff2[sliceSortedIds2[i]]]) == minmaxpiv &&
					// 			!sliceIsDone2[sliceSortedIds2[i]] &&
					// 			ref.strandTable[sliceOff2[sliceSortedIds2[i]]] > minmaxstrand)
					// 	{
					// 		//TODO
					// 		// of << "WrongStrand++" << KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff[sliceSortedIds[i]]]) << "\n";
					// 		++sliceOff2[sliceSortedIds2[i]];
					// 		// reached end of vector slice
					// 		if (sliceOff2[sliceSortedIds2[i]] >= sliceEnd2[sliceSortedIds1[i]])
					// 		{
					// 			//TODO
					// 			// of << "Setting " << sliceSortedIds[i] << " to End\n";
					// 			sliceIsDone2[sliceSortedIds2[i]] = true;
					// 		}
					// 	}
					// }
					// // if nothing has changed, test for match
					// if (unchanged)
					// {
                    //
					// 	// test if any of the windows adjacent to p1 have enough k-mers
					// 	std::vector<uint16_t> candCounts(2*contextWLen + 1, 0);
                    //
					// 	for (unsigned int lIdx = 0; lIdx < kmerNum2; ++lIdx)
					// 	{
					// 		if (sliceIsDone2[lIdx])
					// 			continue;
					// 		for (unsigned int off = 0; off < 2*(2*contextWLen + 1) && sliceOff2[lIdx] + off < sliceEnd2[lIdx]; ++off)
					// 		{
					// 			if (ref.kmerTableSmall[sliceOff2[lIdx] + off] < qWindow1 - 2)
					// 				break;
                    //
					// 			// need to match on opposing strands
					// 			if (ref.strandTable[sliceOff2[lIdx] + off] == minmaxstrand)
					// 				++candCounts[qWindow1 - ref.kmerTableSmall[sliceOff2[lIdx] + off] + contextWLen];
                    //
					// 		}
					// 	}
                    //
                    //
                    //
                        //
                        //
						// // startTime = std::chrono::high_resolution_clock::now();
						// for (unsigned int candIdx = 0; candIdx < candCounts.size(); ++candIdx)
						// {
						// 	if (candCounts[candIdx] < qThreshold)
						// 		continue;
                        //
						// 	const uint32_t metaID1 = qWindow1;
						// 	const uint32_t metaID2 = qWindow1 - contextWLen + candIdx;
                        //
						// 	const bool metaStrand1 = !minmaxstrand;
						// 	const bool metaStrand2 = minmaxstrand;
						// 	// of << "\nMatching Meta:\n\t" << metaID1 << "/" << metaID2 << "\t\t <--- \n";
                        //
						// 	const struct CpG& startCpg1 = ref.cpgTable[ref.metaCpGs[metaID1].start];
						// 	const struct CpG& endCpg1 = ref.cpgTable[ref.metaCpGs[metaID1].end];
						// 	auto startIt1 = ref.fullSeq[startCpg1.chrom].begin() + startCpg1.pos;
						// 	auto endIt1 = ref.fullSeq[startCpg1.chrom].begin() + endCpg1.pos + (2*MyConst::READLEN - 2) + MyConst::MISCOUNT + MyConst::ADDMIS;
                        //
						// 	// check if CpG was too near to the end
						// 	if (endIt1 > ref.fullSeq[startCpg1.chrom].end())
						// 	{
						// 		// if so move end iterator appropriately
						// 		endIt1 = ref.fullSeq[startCpg1.chrom].end();
						// 	}
                        //
						// 	// use shift and to find all matchings
						// 	std::vector<uint64_t> matchings1;
						// 	std::vector<uint8_t> errors1;
						// 	// endTime = std::chrono::high_resolution_clock::now();
						// 	// runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
						// 	// of << "\nMatch preprocessing runtime " << runtime << "\t\t";
						// 	// startTime = std::chrono::high_resolution_clock::now();
						// 	if (metaStrand1)
						// 	{
						// 		sa1.querySeq(startIt1, endIt1, matchings1, errors1);
						// 	} else {
						// 		--startIt1, --endIt1;
						// 		sa1.queryRevSeq(endIt1, startIt1, matchings1, errors1);
						// 	}
						// 	// endTime = std::chrono::high_resolution_clock::now();
						// 	// runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
						// 	// of << "\nShiftAnd runtime " << runtime << "\t\t";
                        //
						// 	// startTime = std::chrono::high_resolution_clock::now();
                        //
						// 	// TODO: match of second read
                        //
                        //
                        //
                        //
                        //
                        //
						// 	// size_t i = 0;
						// 	// // compare first found match with last found match of previous meta CpG
						// 	// if (matchings1.size() > 0)
						// 	// {
						// 	// 	// compare chromosome and offset
						// 	// 	if (matchings1[0] + ref.cpgTable[ref.metaCpGs[metaID1].start].pos == prevOff && ref.cpgTable[ref.metaCpGs[metaID1].start].chrom == prevChr)
						// 	// 	{
						// 	// 		++i;
						// 	// 	}
						// 	// }
						// 	// // go through matching and see if we had such a match (with that many errors) before - if so,
						// 	// // return to caller reporting no match
						// 	// for (; i < matchings1.size(); ++i)
						// 	// {
                        //     //
						// 	// 	// check if we had a match with that many errors before
						// 	// 	if (multiMatch[errors1[i]])
						// 	// 	{
                        //     //
						// 	// 		MATCH::match& match_2 = uniqueMatches[errors[i]];
						// 	// 		// const bool isStart = MATCH::isStart(match_2);
						// 	// 		const bool isFwd = MATCH::isFwd(match_2);
						// 	// 		// check if same k-mer (borders of meta CpGs)
						// 	// 		if (ref.cpgTable[ref.metaCpGs[MATCH::getMetaID(match_2)].start].pos + MATCH::getOffset(match_2) == startCpg.pos + matchings[i])
						// 	// 		{
						// 	// 			if ((isFwd && metaStrand) || (!isFwd && !metaStrand))
						// 	// 				continue;
                        //     //
						// 	// 		} else {
                        //     //
						// 	// 			// check if this is a match without errors
						// 	// 			if (!errors[i])
						// 	// 			{
                        //     //
						// 	// 				// if so, return without a match
						// 	// 				// of << "\tNonunique <---\n\n";
						// 	// 				return -1;
                        //     //
						// 	// 			}
						// 	// 			// set the number of matches with that many errors to 2
						// 	// 			// indicating that we do not have a unique match with that many errors
						// 	// 			multiMatch[errors[i]] = 2;
						// 	// 		}
                        //     //
                        //     //
						// 	// 	} else {
                        //     //
						// 	// 		// update qgram lemma
						// 	// 		uint16_t newQ = sa.size() - MyConst::KMERLEN - (MyConst::KMERLEN * errors[i]);
						// 	// 		// check for overflow and if we improved old q
						// 	// 		if (newQ < sa.size() && newQ > qThreshold)
						// 	// 			qThreshold = newQ;
                        //     //
                        //     //
						// 	// 		// we don't have such a match yet,
						// 	// 		// so save this match at the correct position
						// 	// 		if (metaStrand)
						// 	// 		{
						// 	// 			uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i], errors[i], 1, 0, metaID);
						// 	// 		} else {
						// 	// 			uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i], errors[i], 0, 0, metaID);
						// 	// 		}
						// 	// 		multiMatch[errors[i]] = 1;
						// 	// 	}
						// 	// }
						// 	// if (matchings.size() > 0)
						// 	// {
                        //     //
						// 	// 	prevChr = ref.cpgTable[ref.metaCpGs[metaID].start].chrom;
						// 	// 	prevOff = ref.cpgTable[ref.metaCpGs[metaID].start].pos + matchings[matchings.size() - 1];
                        //     //
						// 	// } else {
                        //     //
						// 	// 	prevChr = 0;
						// 	// 	prevOff = 0xffffffffffffffffULL;
						// 	// }
						// 	// endTime = std::chrono::high_resolution_clock::now();
						// 	// runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
						// 	// of << "\nPostprocessing matches runtime " << runtime << "\t\t";
						// 	// startTime = std::chrono::high_resolution_clock::now();
                        //
						// 	// Advance all lists of read 1 to pivot1 - 1 and all lists of read 2 to pivot1+contextLen-1
						// 	for (size_t i = 0; i < kmerNum1; ++i)
						// 	{
						// 		if (sliceIsDone1[sliceSortedIds1[i]])
						// 			continue;
						// 		if (KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff1[sliceSortedIds1[i]]]) == qWindow1)
						// 		{
						// 			if (ref.strandTable[sliceOff1[sliceSortedIds1[i]]] == !minmaxstrand)
						// 			{
						// 				// of << sliceSortedIds[i] << "\t";
						// 				// of << KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff[sliceSortedIds[i]]]) << "++";
						// 				++sliceOff1[sliceSortedIds1[i]];
						// 				// reached end of vector slice
						// 				if (sliceOff1[sliceSortedIds1[i]] >= sliceEnd1[sliceSortedIds1[i]])
						// 				{
						// 					// of << "Setting " << sliceSortedIds[i] << " to End\n";
						// 					sliceIsDone1[sliceSortedIds1[i]] = true;
						// 				}
						// 			}
						// 		}
						// 	}
						// 	for (size_t i = 0; i < kmerNum2; ++i)
						// 	{
						// 		if (sliceIsDone2[sliceSortedIds2[i]])
						// 			continue;
			// 					if (KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff2[sliceSortedIds2[i]]]) > qWindow1 + contextWLen - 1)
			// 					{
			// 						// of << sliceSortedIds[i] << "\t";
			// 						// of << KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff[sliceSortedIds[i]]]) << "++";
			// 						++sliceOff2[sliceSortedIds2[i]];
			// 						// reached end of vector slice
			// 						if (sliceOff2[sliceSortedIds2[i]] >= sliceEnd2[sliceSortedIds2[i]])
			// 						{
			// 							// of << "Setting " << sliceSortedIds[i] << " to End\n";
			// 							sliceIsDone2[sliceSortedIds2[i]] = true;
			// 						}
			// 					}
			// 				}
			// 				// endTime = std::chrono::high_resolution_clock::now();
			// 				// runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
			// 				// of << "Advance after match runtime " << runtime << "\t\t";
			// 			}
			// 		}
			// 	// symmetric case qWindow2 < qWindow 1
			// 	} else {
			// 		// TODO
			// 	}
            //
			// 	// order descending on window and reverse strand > fwd strand
			// 	// Order:
			// 	//		return false (i.e. id1 > id2) if
			// 	//			has no windows left to process in 1
			// 	//		return true (i.e. id1 < id2) if
			// 	//			window id 1 > window id 2 and
			// 	//				strand is rev for 1 and is fwd for 2
			// 	//			else return false
			// 	// startTime = std::chrono::high_resolution_clock::now();
			// 	sort_by_n(sliceSortedIds1.begin(), sliceSortedIds1.begin() + qThreshold - 1, sliceSortedIds1.end(), sliceOff1, sliceIsDone1);
			// 	sort_by_n(sliceSortedIds2.begin(), sliceSortedIds2.begin() + qThreshold - 1, sliceSortedIds2.end(), sliceOff2, sliceIsDone2);
			// 	// endTime = std::chrono::high_resolution_clock::now();
			// 	// runtime = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
			// 	// of << "Sort runtime " << runtime << "\n\n";
			// 	// of << "SortedIds:\n\t";
			// 	// for (size_t i = 0; i < kmerNum; ++i)
			// 	// {
			// 	// 	of << sliceSortedIds[i] << "\t";
			// 	// }
			// 	// of << "\n";
        //     //
		// 	// }
        //     //
		// 	// of << "qGrams at end: \n";
		// 	// for (size_t i = 0; i < kmerNum; ++i)
		// 	// {
		// 	// 	if (!sliceIsDone[sliceSortedIds[qThreshold-1]])
		// 	// 		of << "\t" << KMER_S::getMetaCpG(ref.kmerTableSmall[sliceOff[sliceSortedIds[i]]]) << "\n";
		// 	// 	else
		// 	// 		of << "\tEnd\n";
		// 	// }
        //
        //     // go through found matches for each [0,maxErrorNumber] and see if it is unique
        //     // for (size_t i = 0; i < multiMatch.size(); ++i)
        //     // {
        //     //     // there is no match with that few errors, search the one with more errors
        //     //     if (multiMatch[i] == 0)
        //     //     {
        //     //         continue;
        //     //     }
        //     //     mat = uniqueMatches[i];
        //     //     // if match is not unique, return unsuccessfull to caller
        //     //     if (multiMatch[i] > 1)
        //     //     {
        //     //
		// 	// 		// of << "\tNonunique <---\n\n";
        //     //         return -1;
        //     //
        //     //     // exactly one with that many errors - return successfull
        //     //     } else {
        //     //
		// 	// 		// of << "\tSuccessfull match <---\n\n";
        //     //         return 1;
        //     //     }
        //     //
        //     // }
        //     // we have not a single match at all, return unsuccessfull to caller
		// 	// of << "\tNo match <---\n\n";
        //     return 0;
		// }



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

        // TODO: paired end DS for this????
        std::array<std::vector<uint16_t>, CORENUM> countsFwdStart;
        std::array<std::vector<uint16_t>, CORENUM> countsRevStart;
        std::array<google::dense_hash_map<uint32_t, uint16_t, MetaHash>, CORENUM> fwdMetaIDs;
        std::array<google::dense_hash_map<uint32_t, uint16_t, MetaHash>, CORENUM> revMetaIDs;
        // Holds counts for each thread for counting heuristic
        // KEY: Meta CpG ID
        // VALUE:
        //      1) K-mer count of first read
        //      2) K-mer count of second read
		//      3) Boolena flag that is ture iff second read has enough kmers in this or adjacent MetaCpGs
        //      4) Boolean flag that is true iff first read is matched to this or adjacent Meta CpGs
        //
        std::array<google::dense_hash_map<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>, MetaHash>, CORENUM> paired_fwdMetaIDs;
        std::array<google::dense_hash_map<uint32_t, std::tuple<uint8_t, uint8_t, bool, bool>, MetaHash>, CORENUM> paired_revMetaIDs;

        bool isPaired;

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

        // holds matching info
        std::array<uint64_t, CORENUM> matchStats;
        std::array<uint64_t, CORENUM> nonUniqueStats;
        std::array<uint64_t, CORENUM> noMatchStats;
        std::array<uint64_t, CORENUM> matchPairedStats;
        std::array<uint64_t, CORENUM> tooShortCounts;

		// buffer arrays for Karl's matching
		std::array<std::vector<uint64_t>, CORENUM> sliceOffThreads;
		std::array<std::vector<uint64_t>, CORENUM> sliceEndThreads;
		std::array<std::vector<unsigned int>, CORENUM> sliceSortedIdsThreads;
		std::array<std::vector<bool>, CORENUM> sliceIsDoneThreads;

		bool bothStrandsFlag;
		// counter for read 1 matches to fwd strand
		uint64_t r1FwdMatches;
		// counter for read 1 matches to rev strand
		uint64_t r1RevMatches;
		bool matchR1Fwd;

        // TODO
        std::ofstream of;
        inline void printMatch(std::ofstream& o, MATCH::match& mat)
        {
            o << "Match at offset " << MATCH::getOffset(mat) + ref.cpgTable[ref.metaCpGs[MATCH::getMetaID(mat)].start].pos << \
                " on chromosome " << static_cast<uint64_t>(ref.cpgTable[ref.metaCpGs[MATCH::getMetaID(mat)].start].chrom) << \
                " on strand " << (MATCH::isFwd(mat) ? "fwd" : "rev") << " with " << MATCH::getErrNum(mat) << " many errors. Meta CpG " << MATCH::getMetaID(mat);
        }


};

#endif /* READQUEUE_H */
