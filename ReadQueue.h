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
#include <unordered_map>
#include <unordered_set>
#include <array>
#include <algorithm> // reverse

#ifdef _OPENMP
#include <omp.h>
#endif

#include <sparsehash/dense_hash_map>

// #include "sparsepp/spp.h"

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
        ReadQueue(const char* filePath, RefGenome& ref, bool isGZ);

        // for paired end
        // ARGUMENTS:
        //          filePath    path to file containing read1 of paired reads in fastq format
        //          filePath2   path to file containing read2 of paired reads in fastq format
        //          ref         internal representation of reference genome
        //          isGZ        flag - true iff file is gzipped
        //
        // NOTE:
        //          provided files are ASSUMED to have equal number of reads in correct (paired) order!
        ReadQueue(const char* filePath, const char* filePath2, RefGenome& reference, bool isGZ);

        // -----------

        // Parses a chunk of the ifstream file/ igzstream file, respectively
        // Reads up to MyConst::CHUNKSIZE many reads and saves them
        // ARGUMENT:
        //          procReads   will contain number of reads that have been read into buffer
        // returns true if neither read error nor EOF occured, false otherwise
        bool parseChunk(unsigned int& procReads);
        bool parseChunkGZ(unsigned int& procReads);


        // Match all reads in readBuffer to reference genome
        // ARGUMENT:
        //          procReads   number of reads to match
        // First retrieve seeds using getSeeds(...)
        // Filter seeds using filterHeuSeeds(...) according to simple heuristic
        // Make full match using bit shift trick
        // Align match using banded levenshtein alignment to update methylation counts
        bool matchReads(const unsigned int& procReads, uint64_t& succMatch, uint64_t& nonUniqueMatch, uint64_t& unSuccMatch);
        bool matchPairedReads(const unsigned int& procReads, uint64_t& succMatch, uint64_t& nonUniqueMatch, uint64_t& unSuccMatch);

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

        // filters seeds according to simple counting criteria
        // #kmers of one metaCpG should be > READLEN - KMERLEN + 1 - (KMERLEN * MISCOUNT)
        inline void filterHeuSeeds(std::vector<std::vector<KMER_S::kmer> >& seedsK, std::vector<std::vector<bool> >& seedsS, const unsigned int readSize)
        {

            std::vector<uint16_t>& threadCountFwd = countsFwd[omp_get_thread_num()];
            std::vector<uint16_t>& threadCountRev = countsRev[omp_get_thread_num()];
            // fill with zeroes
            threadCountFwd.assign(ref.metaCpGs.size(), 0);
            threadCountRev.assign(ref.metaCpGs.size(), 0);

            // count occurences of meta CpGs
            for (unsigned int i = 0; i < seedsK.size(); ++i)
            {

                // last visited id in this table entry
                // avoid counting metaCpGs more then once per kmer
                // note that metaCpGs are hashed in reverse order
                uint64_t lastId = 0xffffffffffffffffULL;
                // strand of last visited id (true iff forward strand)
                bool lastStrand = false;

                for (size_t j = 0; j < seedsK[i].size(); ++j)
                {

                    const uint64_t metaId = KMER_S::getMetaCpG(seedsK[i][j]);
                    const bool metaStrand = seedsS[i][j];
                    // check if we visited meta CpG before
                    if (metaId == lastId && metaStrand == lastStrand)
                    {
                        continue;
                    }

                    lastId = metaId;
                    lastStrand = metaStrand;
                    if (metaStrand)
                    {
                        ++threadCountFwd[metaId];

                    } else {

                        ++threadCountRev[metaId];

                    }
                }
            }

            // More than cutoff many kmers are required per metaCpG
            // unsigned int countCut = readSize - MyConst::KMERLEN + 1 - (MyConst::KMERLEN * MyConst::MISCOUNT);
            // unsigned int countCut = readSize - MyConst::KMERLEN + 1 - (MyConst::KMERLEN * MyConst::MISCOUNT) - 10;
            // check if we overflowed
            // if (countCut > readSize)
            // {
            //     countCut = 0;
            // }
            uint16_t countCut = 20;


            // throw out rare metaCpGs
            for (size_t i = 0; i < seedsK.size(); ++i)
            {

                // iterator to the element that we process
                auto srcItK = seedsK[i].begin();
                auto srcItS = seedsS[i].begin();
                // iterator to the position one of the last inserted FILTERED element, always at most as far as srcIt
                auto filterItK = seedsK[i].begin();
                auto filterItS = seedsS[i].begin();

                for (size_t j = 0; j < seedsK[i].size(); ++j, ++srcItK, ++srcItS)
                {

                    // check strand
                    if (seedsS[i][j])
                    {
                        // test for strict heuristic criterias
                        if (threadCountFwd[KMER_S::getMetaCpG(seedsK[i][j])] >= countCut)
                        {

                            *filterItK = *srcItK;
                            *filterItS = *srcItS;
                            ++filterItK;
                            ++filterItS;

                        }

                    } else {

                        // test for strict heuristic criterias
                        if (threadCountRev[KMER_S::getMetaCpG(seedsK[i][j])] >= countCut)
                        {

                            *filterItK = *srcItK;
                            *filterItS = *srcItS;
                            ++filterItK;
                            ++filterItS;

                        }
                    }
                }
                seedsK[i].resize(filterItK - seedsK[i].begin());
                seedsS[i].resize(filterItK - seedsK[i].begin());
            }

        }
        // inline void filterHeuSeedsRef(std::vector<size_t>& seedsK, const unsigned int readSize)
        // {
        //
        //     std::vector<uint16_t>& threadCountFwd = countsFwd[omp_get_thread_num()];
        //     std::vector<uint16_t>& threadCountRev = countsRev[omp_get_thread_num()];
        //     std::vector<uint16_t>& threadCountFwdStart = countsFwdStart[omp_get_thread_num()];
        //     std::vector<uint16_t>& threadCountRevStart = countsRevStart[omp_get_thread_num()];
        //     // fill with zeroes
        //     threadCountFwd.assign(ref.metaCpGs.size(), 0);
        //     threadCountRev.assign(ref.metaCpGs.size(), 0);
        //     threadCountFwdStart.assign(ref.metaStartCpGs.size(), 0);
        //     threadCountRevStart.assign(ref.metaStartCpGs.size(), 0);
        //
        //     // count occurences of meta CpGs
        //     for (size_t i : seedsK)
        //     {
        //
        //         // last visited id in this table entry
        //         // avoid counting metaCpGs more then once per kmer
        //         // note that metaCpGs are hashed in reverse order
        //         uint64_t lastId = 0xffffffffffffffffULL;
        //         // strand of last visited id (true iff forward strand)
        //         bool wasFwd = false;
        //         bool wasStart = false;
        //
        //         for (size_t j = ref.tabIndex[i]; j < ref.tabIndex[i+1]; ++j)
        //         {
        //
        //             // retrieve seed information
        //             KMER_S::kmer& k = ref.kmerTableSmall[j];
        //             const bool isFwd = ref.strandTable[j];
        //             const uint64_t metaId = KMER_S::getMetaCpG(k);
        //             const bool isStart = KMER_S::isStartCpG(k);
        //             // check if we visited meta CpG before
        //             if (metaId == lastId && isFwd == wasFwd && isStart == wasStart)
        //             {
        //                 continue;
        //             }
        //
        //             // update vars for last checked metaCpG
        //             lastId = metaId;
        //             wasFwd = isFwd;
        //             wasStart = isStart;
        //             if (isStart)
        //             {
        //                 if (isFwd)
        //                 {
        //                     ++threadCountFwdStart[metaId];
        //
        //                 } else {
        //
        //                     ++threadCountRevStart[metaId];
        //
        //                 }
        //
        //             } else {
        //
        //                 if (isFwd)
        //                 {
        //                     ++threadCountFwd[metaId];
        //
        //                 } else {
        //
        //                     ++threadCountRev[metaId];
        //
        //                 }
        //             }
        //         }
        //     }
        // }

        // Do a bitmatching between the specified seeds of the reference and the read r or the reverse complement (Rev suffix)
        //
        // ARGUMENTS:
        //              r       read to match with
        //              seedsK  list of kmer positions in reference that should be checked
        //              seedsS  list of flags for each kmer in seedsK stating if it is from forward or reverse reference strand
        //
        // RETURN:
        //              void
        //
        // MODIFICATIONS:
        //              function will filter seedsK and seedsS to contain only reference kmers that match read kmer under bitmask
        //              comparison
        // inline void bitMatching(const Read& r, std::vector<std::vector<KMER_S::kmer> >& seedsK, std::vector<std::vector<bool> >& seedsS)
        // {
        //
        //     // masking for actual kmer bits
        //     constexpr uint64_t signiBits = 0xffffffffffffffffULL >> (64 - (2*MyConst::KMERLEN));
        //     // bit representation of current kmer of read
        //     uint64_t kmerBits = 0;
        //     // generate bit representation of first kmerlen - 1 letters of read
        //     for (unsigned int i = 0; i < (MyConst::KMERLEN - 1); ++i)
        //     {
        //
        //         kmerBits = kmerBits << 2;
        //         kmerBits |= BitFun::getBitRepr(r.seq[i]);
        //
        //     }
        //
        //     // Note that seedsK must have the same size as seedsS anyway, this way we may have a cache hit
        //     std::vector<std::vector<KMER_S::kmer> > newSeedsK(seedsK.size());
        //     std::vector<std::vector<bool> > newSeedsS(seedsK.size());
        //
        //     // go over each read kmer and compare with reference seeds
        //     for (unsigned int offset = 0; offset < (r.seq.size() - MyConst::KMERLEN + 1); ++offset)
        //     {
        //
        //         // retrieve seeds for current kmer
        //         std::vector<KMER_S::kmer>& localSeedsK = seedsK[offset];
        //         std::vector<bool>& localSeedsS = seedsS[offset];
        //         // reserve some space for result seedlist
        //         newSeedsK[offset].reserve(localSeedsK.size());
        //         newSeedsS[offset].reserve(localSeedsK.size());
        //
        //         // update current read kmer representation
        //         kmerBits = kmerBits << 2;
        //         kmerBits = (kmerBits | BitFun::getBitRepr(r.seq[offset + MyConst::KMERLEN - 1])) & signiBits;
        //
        //         // iterate over corresponding seeds for reference
        //         for (unsigned int i = 0; i < localSeedsK.size(); ++i)
        //         {
        //
        //             KMER_S::kmer& refKmer = localSeedsK[i];
        //
        //             // will hold the first CpG in metaCpG after retrieving the meta CpG info
        //             uint8_t chrom;
        //             // will hold the position of the meta CpG in the genome (chromosomal region specified by chrom of the genome)
        //             uint32_t pos = 0;
        //             // retrieve reference bit representation
        //             //
        //             // First retrieve meta CpG info
        //             if (KMER_S::isStartCpG(refKmer))
        //             {
        //
        //                 const uint32_t& cpgInd = ref.metaStartCpGs[KMER_S::getMetaCpG(refKmer)].start;
        //                 chrom = ref.cpgStartTable[cpgInd].chrom;
        //
        //             } else {
        //
        //                 const uint32_t& cpgInd = ref.metaCpGs[KMER_S::getMetaCpG(refKmer)].start;
        //                 chrom = ref.cpgTable[cpgInd].chrom;
        //                 pos = ref.cpgTable[cpgInd].pos;
        //
        //             }
        //             // will hold bit representation of seed
        //             uint64_t refKmerBit;
        //             // decide if forward or reverse strand of reference
        //             //
        //             // if forward
        //             if (localSeedsS[i])
        //             {
        //                 // retrieve sequence in forward strand
        //                 refKmerBit = ref.genomeBit[chrom].getSeqKmer(pos + KMER_S::getOffset(refKmer));
        //
        //             // is reverse
        //             } else {
        //
        //                 // retrieve sequence in forward strand
        //                 refKmerBit = ref.genomeBit[chrom].getSeqKmerRev(pos + KMER_S::getOffset(refKmer));
        //             }
        //
        //             // COMPARE read kmer and seed
        //             //  matching is 0 iff full match
        //             if ( !( refKmerBit ^ (kmerBits & BitFun::getMask(refKmerBit)) ) )
        //             {
        //
        //                 // if we have a match, keep this kmer and strand flag in list
        //                 newSeedsK[offset].emplace_back(refKmer);
        //                 newSeedsS[offset].emplace_back(localSeedsS[i]);
        //             }
        //         }
        //         newSeedsK[offset].shrink_to_fit();
        //         newSeedsS[offset].shrink_to_fit();
        //
        //     }
        //
        //     seedsK = std::move(newSeedsK);
        //     seedsS = std::move(newSeedsS);
        //
        // }
        // inline void bitMatchingRev(const Read& r, std::vector<std::vector<KMER_S::kmer> >& seedsK, std::vector<std::vector<bool> >& seedsS)
        // {
        //
            // const unsigned int readSize = r.seq.size();
            // // masking for actual kmer bits
            // constexpr uint64_t signiBits = 0xffffffffffffffffULL >> (64 - (2*MyConst::KMERLEN));
            // // bit representation of current kmer of read
            // uint64_t kmerBits = 0;
            // // generate bit representation of first kmerlen - 1 letters of reverse complement of read
            // // Not that we start reading from right
            // for (unsigned int i = readSize - 1; i > readSize - MyConst::KMERLEN; --i)
            // {
            //
            //     kmerBits = kmerBits << 2;
            //     kmerBits |= BitFun::getBitReprRev(r.seq[i]);
            //
            // }
            //
            // // Note that seedsK must have the same size as seedsS anyway, this way we may have a cache hit
            // std::vector<std::vector<KMER_S::kmer> > newSeedsK(seedsK.size());
            // std::vector<std::vector<bool> > newSeedsS(seedsK.size());
            //
            // // index for seed vector for current read kmer
            // unsigned int kmerInd = 0;
            // // go over each read kmer and compare with reference seeds
            // for (unsigned int offset = readSize - MyConst::KMERLEN; offset > 0; --offset, ++kmerInd)
            // {
            //
            //     std::vector<KMER_S::kmer>& localSeedsK = seedsK[kmerInd];
            //     std::vector<bool>& localSeedsS = seedsS[kmerInd];
            //     // reserve some space for seedlist
            //     newSeedsK[kmerInd].reserve(localSeedsK.size());
            //     newSeedsS[kmerInd].reserve(localSeedsK.size());
            //
            //     // update current read kmer representation
            //     kmerBits = kmerBits << 2;
            //     kmerBits = (kmerBits | BitFun::getBitReprRev(r.seq[offset - 1])) & signiBits;
            //
            //     // iterate over corresponding seeds for reference
            //     for (unsigned int i = 0; i < localSeedsK.size(); ++i)
            //     {
            //
            //         KMER_S::kmer& refKmer = localSeedsK[i];
            //
            //         // will hold the chromosome index in metaCpG after retrieving the meta CpG info
            //         uint8_t chrom;
            //         // will hold the position of the meta CpG in the genome (chromosomal region specified by chrom of the genome)
            //         uint32_t pos = 0;
            //         // retrieve reference bit representation
            //         //
            //         // First retrieve meta CpG info
            //         if (KMER_S::isStartCpG(refKmer))
            //         {
            //
            //             uint32_t& cpgInd = ref.metaStartCpGs[KMER_S::getMetaCpG(refKmer)].start;
            //             chrom = ref.cpgStartTable[cpgInd].chrom;
            //
            //         } else {
            //
            //             uint32_t& cpgInd = ref.metaCpGs[KMER_S::getMetaCpG(refKmer)].start;
            //             chrom = ref.cpgTable[cpgInd].chrom;
            //             pos = ref.cpgTable[cpgInd].pos;
            //
                    // }
                    // will hold bit representation of seed
                    // uint64_t refKmerBit;
                    // decide if forward or reverse strand of reference
                    //
                    // if forward
        //             if (localSeedsS[i])
        //             {
        //                 // retrieve sequence in forward strand
        //                 refKmerBit = ref.genomeBit[chrom].getSeqKmer(pos + KMER_S::getOffset(refKmer));
        //
        //             // is reverse
        //             } else {
        //
        //                 // retrieve sequence in forward strand
        //                 refKmerBit = ref.genomeBit[chrom].getSeqKmerRev(pos + KMER_S::getOffset(refKmer));
        //             }
        //
        //             // COMPARE read kmer and seed
        //             //  matching is 0 iff full match
        //             if ( !( refKmerBit ^ (kmerBits & BitFun::getMask(refKmerBit)) ) )
        //             {
        //
        //                 // if we have a match, keep this kmer and strand flag in list
        //                 newSeedsK[kmerInd].emplace_back(refKmer);
        //                 newSeedsS[kmerInd].emplace_back(localSeedsS[i]);
        //             }
        //         }
        //         newSeedsK[kmerInd].shrink_to_fit();
        //         newSeedsS[kmerInd].shrink_to_fit();
        //
        //     }
        //
        //     seedsK = std::move(newSeedsK);
        //     seedsS = std::move(newSeedsS);
        // }
        //

        // query seeds to a given shift-and automaton
        //
        // ARGUMENTS:
        //              sa          shift-and automaton
        //              seedsK      k-mer representation of seeds
        //              seedsS      strandedness of seeds
        //              mat         will contain a unique match if successfull (i.e. return is true)
        //
        // RETURN:
        //              0 if no match found
        //              1 if match is found
        //              -1 if multiple matches, no unique best
        //
        // MODIFICATIONS:
        //              updates the match mat
        //              shift-and automaton is used - needs a reset (implicitly done when querying a new sequence internally)
        //
        // inline int saQuerySeedSet(ShiftAnd<MyConst::MISCOUNT>& sa, std::vector<std::vector<KMER_S::kmer> >& seedsK, std::vector<std::vector<bool> >& seedsS, MATCH::match& mat)
        // {

            // use counters to flag what has been processed so far
            // 0 if not processed at all
            // 1 if start metaCpG has been processed before
            // 2 if metaCpG with this index has been proc. before
            // 3 if both
            // std::vector<uint16_t>& threadCountFwd = countsFwd[omp_get_thread_num()];
            // std::vector<uint16_t>& threadCountRev = countsRev[omp_get_thread_num()];
            //
            // // fill with zeroes
            // threadCountFwd.assign(ref.metaCpGs.size(), 0);
            // threadCountRev.assign(ref.metaCpGs.size(), 0);
            //
            // // counter for how often we had a match
            // std::array<uint8_t, MyConst::MISCOUNT + 1> multiMatch;
            // multiMatch.fill(0);
            //
            // // will contain matches iff match is found for number of errors specified by index
            // std::array<MATCH::match, MyConst::MISCOUNT + 1> uniqueMatches;
            //
            // for (size_t outerNdx = 0; outerNdx < seedsK.size(); ++outerNdx)
            // {
            //
            //     for (size_t innerNdx = 0; innerNdx < seedsK[outerNdx].size(); ++innerNdx)
            //     {
            //
            //         // retrieve kmer
            //         KMER_S::kmer& k = seedsK[outerNdx][innerNdx];
            //
            //         // retrieve meta CpG
            //         const uint64_t m = KMER_S::getMetaCpG(k);
            //
            //         // retrieve if start meta CpG
            //         const bool isStart = KMER_S::isStartCpG(k);
            //
            //         // retrieve strand
            //         const bool isFwd = seedsS[outerNdx][innerNdx];
            //
            //         if (isStart)
            //         {
            //
            //             if (isFwd)
            //             {
            //                 // check if we queried this meta CpG it before
            //                 if (threadCountFwd[m] == 1 || threadCountFwd[m] >= 3)
            //                 {
            //                     continue;
            //
            //                 } else {
            //
            //                     // state that we queried this
            //                     threadCountFwd[m] += 1;
            //                     // retrieve it
            //                     const struct CpG& startCpg = ref.cpgStartTable[ref.metaStartCpGs[m].start];
            //                     const struct CpG& endCpg = ref.cpgStartTable[ref.metaStartCpGs[m].end];
            //
            //                     auto startIt = ref.fullSeq[startCpg.chrom].begin();
            //                     auto endIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2);
            //
            //                     // check if CpG was too near to the end
            //                     if (endIt > ref.fullSeq[startCpg.chrom].end())
            //                     {
            //                         // if so move end iterator appropriately
            //                         endIt = ref.fullSeq[startCpg.chrom].end();
            //                     }
            //
                        //         // use shift and to find all matchings
                        //         std::vector<uint64_t> matchings;
                        //         std::vector<uint8_t> errors;
                        //         sa.querySeq(startIt, endIt, matchings, errors);
                        //
                        //         // go through matching and see if we had such a match (with that many errors) before - if so,
                        //         // return to caller reporting no match
                        //         for (size_t i = 0; i < matchings.size(); ++i)
                        //         {
                        //
                        //             // check if we had a match with that many errors before
                        //             if (multiMatch[errors[i]])
                        //             {
                        //
                        //                 MATCH::match& match_2 = uniqueMatches[errors[i]];
                        //                 // check if same k-mer (borders of meta CpGs)
                        //                 if ((MATCH::getChrom(match_2) == startCpg.chrom) && (MATCH::getOffset(match_2) == matchings[i]) && (MATCH::isFwd(match_2)))
                        //                 {
                        //                     continue;
                        //
                        //                 } else {
                        //
                        //                     // check if this is a match without errors
                        //                     if (!errors[i])
                        //                     {
                        //
                        //                         // if so, return without a match
                        //                         return -1;
                        //
                        //                     }
                        //                     // set the number of matches with that many errors to 2
                        //                     // indicating that we do not have a unique match with that many errors
                        //                     multiMatch[errors[i]] = 2;
                        //                 }
                        //
                        //
                        //             } else {
                        //
                        //                 // we don't have such a match yet,
                        //                 // so save this match at the correct position
                        //                 uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i], startCpg.chrom, errors[i], 1);
                        //                 multiMatch[errors[i]] = 1;
                        //             }
                        //
                        //
                        //         }
                        //
                        //     }
                        //
                        // // is not forward strand meta CpG
                        // } else {
                        //
                        //     // check if we queried it before
                        //     if (threadCountRev[m] == 1 || threadCountRev[m] >= 3)
                        //     {
                        //         continue;
                        //
                        //     } else {
                        //
                        //         threadCountRev[m] += 1;
                        //         // retrieve it
                        //         const struct CpG& startCpg = ref.cpgStartTable[ref.metaStartCpGs[m].start];
                        //         const struct CpG& endCpg = ref.cpgStartTable[ref.metaStartCpGs[m].end];
                        //
                        //         auto endIt = ref.fullSeq[startCpg.chrom].begin() - 1;
                        //         auto startIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2) - 1;
                        //
                        //         // check if CpG was too near to the end
                        //         if (startIt >= ref.fullSeq[startCpg.chrom].end())
                        //         {
                        //             // if so move end iterator appropriately
                        //             startIt = ref.fullSeq[startCpg.chrom].end() - 1;
                        //         }
                        //
                        //         // use shift and to find all matchings
                        //         std::vector<uint64_t> matchings;
                        //         std::vector<uint8_t> errors;
                        //         sa.queryRevSeq(startIt, endIt, matchings, errors);
                        //
                        //         // go through matching and see if we had such a match (with that many errors) before - if so,
                        //         // return to caller reporting no match
                        //         for (size_t i = 0; i < matchings.size(); ++i)
                        //         {
                        //
                        //             // check if we had a match with that many errors before
                        //             if (multiMatch[errors[i]])
                        //             {
                        //
                        //                 MATCH::match& match_2 = uniqueMatches[errors[i]];
                        //                 // check if same k-mer (borders of meta CpGs)
                        //                 if ((MATCH::getChrom(match_2) == startCpg.chrom) && (MATCH::getOffset(match_2) == matchings[i]) && !(MATCH::isFwd(match_2)))
                        //                 {
                        //                     continue;
                        //
                        //                 } else {
                        //
                        //                     // check if this is a match without errors
                        //                     if (!errors[i])
                        //                     {
                        //
                    //                             // if so, return without a match
                    //                             return -1;
                    //
                    //                         }
                    //                         // set the number of matches with that many errors to 2
                    //                         // indicating that we do not have a unique match with that many errors
                    //                         multiMatch[errors[i]] = 2;
                    //                     }
                    //
                    //
                    //                 } else {
                    //
                    //                     // we don't have such a match yet,
                    //                     // so save this match at the correct position
                    //                     uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i], startCpg.chrom, errors[i], 0);
                    //                     multiMatch[errors[i]] = 1;
                    //                 }
                    //
                    //             }
                    //
                    //         }
                    //
                    //     }
                    //
                    // // if is not start metaCpG
                    // } else {
                    //
                    //     if (isFwd)
                    //     {
                    //
                    //         // check if we queried this meta CpG it before
                    //         if (threadCountFwd[m] >= 2)
                    //         {
                    //             continue;
                    //
                    //         } else {
                    //
                    //
                    //             // state that we queried this
                    //             threadCountFwd[m] += 2;
                    //             // retrieve it
                    //             const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m].start];
                    //             const struct CpG& endCpg = ref.cpgTable[ref.metaCpGs[m].end];
                    //
                    //             auto startIt = ref.fullSeq[startCpg.chrom].begin() + startCpg.pos;
                    //             auto endIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2);
                    //
                    //             // check if CpG was too near to the end
                    //             if (endIt > ref.fullSeq[startCpg.chrom].end())
                    //             {
                    //                 // if so move end iterator appropriately
                    //                 endIt = ref.fullSeq[startCpg.chrom].end();
                    //             }
                    //
                    //             // use shift and to find all matchings
                    //             std::vector<uint64_t> matchings;
                    //             std::vector<uint8_t> errors;
                    //             sa.querySeq(startIt, endIt, matchings, errors);
                    //
                    //             // go through matching and see if we had such a match (with that many errors) before - if so,
                        //         // return to caller reporting no match
                        //         for (size_t i = 0; i < matchings.size(); ++i)
                        //         {
                        //
                        //             // check if we had a match with that many errors before
                        //             if (multiMatch[errors[i]])
                        //             {
                        //
                        //                 MATCH::match& match_2 = uniqueMatches[errors[i]];
                        //                 // check if same k-mer (borders of meta CpGs)
                        //                 if ((MATCH::getChrom(match_2) == startCpg.chrom) && (MATCH::getOffset(match_2) == matchings[i] + startCpg.pos) && (MATCH::isFwd(match_2)))
                        //                 {
                        //                     continue;
                        //
                        //                 } else {
                        //
                        //                     // check if this is a match without errors
                        //                     if (!errors[i])
                        //                     {
                        //
                        //                         // if so, return without a match
                        //                         // std::cout << "Nonunique no-error match\n";
                        //                         return -1;
                        //
                        //                     }
                        //                     // set the number of matches with that many errors to 2
                        //                     // indicating that we do not have a unique match with that many errors
                        //                     multiMatch[errors[i]] = 2;
                        //                 }
                        //
                        //
                        //             } else {
                        //
                        //                 // we don't have such a match yet,
                        //                 // so save this match at the correct position
                        //                 uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i] + startCpg.pos, startCpg.chrom, errors[i], 1);
                        //                 multiMatch[errors[i]] = 1;
                        //             }
                        //         }
                        //     }
                        //
                        // // kmer is on backward strand
                        // } else {
                        //
                        //     // check if we queried this meta CpG it before
                        //     if (threadCountRev[m] >= 2)
                        //     {
                        //         continue;
                        //
                        //     } else {
                        //
                        //         // state that we queried this
                        //         threadCountRev[m] += 2;
                        //         // retrieve it
                        //         const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m].start];
                        //         const struct CpG& endCpg = ref.cpgTable[ref.metaCpGs[m].end];
                        //
                        //         auto endIt = ref.fullSeq[startCpg.chrom].begin() + startCpg.pos - 1;
                        //         auto startIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2) - 1;
                        //
                        //         // check if CpG was too near to the end
                        //         if (startIt >= ref.fullSeq[startCpg.chrom].end())
                        //         {
                        //             // if so move end iterator appropriately
                        //             startIt = ref.fullSeq[startCpg.chrom].end() - 1;
                        //         }
                        //
                        //         // use shift and to find all matchings
                        //         std::vector<uint64_t> matchings;
                        //         std::vector<uint8_t> errors;
                        //         sa.queryRevSeq(startIt, endIt, matchings, errors);
                        //
                        //         // go through matching and see if we had such a match (with that many errors) before - if so,
                        //         // return to caller reporting no match
                        //         for (size_t i = 0; i < matchings.size(); ++i)
                        //         {
                        //
                        //             // check if we had a match with that many errors before
                        //             if (multiMatch[errors[i]])
                        //             {
                        //
                        //                 MATCH::match& match_2 = uniqueMatches[errors[i]];
                        //                 if ((MATCH::getChrom(match_2) == startCpg.chrom) && (MATCH::getOffset(match_2) == matchings[i] + startCpg.pos) && !(MATCH::isFwd(match_2)))
                        //                 {
                        //                     continue;
                        //
                        //                 } else {
                        //
                        //                     // check if this is a match without errors
                        //                     if (!errors[i])
                        //                     {
                        //
                        //                         // if so, return without a match
                        //                         return -1;
                        //
                        //                     }
                        //                     // set the number of matches with that many errors to 2
                        //                     // indicating that we do not have a unique match with that many errors
                        //                     multiMatch[errors[i]] = 2;
                        //                 }
        //
        //
        //                             } else {
        //
        //                                 // we don't have such a match yet,
        //                                 // so save this match at the correct position
        //                                 uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i] + startCpg.pos, startCpg.chrom, errors[i], 0);
        //                                 multiMatch[errors[i]] = 1;
        //                             }
        //                         }
        //
        //                     } // end else for visited meta CpG test
        //                 } // end else isFwd
        //             } // end else isStart
        //         } // end inner kmer loop
        //     } // end outer kmer loop
        //
        //
        //     // go through found matches for each [0,maxErrorNumber] and see if it is unique
        //     for (size_t i = 0; i < multiMatch.size(); ++i)
        //     {
        //         // there is no match with that few errors, search the one with more errors
        //         if (multiMatch[i] == 0)
        //         {
        //             continue;
        //         }
        //         // if match is not unique, return unsuccessfull to caller
        //         if (multiMatch[i] > 1)
        //         {
        //             // of << "Too bad, multimatch in internal\n";
        //             return -1;
        //
        //         } else {
        //
        //             mat = uniqueMatches[i];
        //             return 1;
        //         }
        //
        //     }
        //     // we have not a single match at all, return unsuccessfull to caller
        //     // of << "No match at all\t";
        //     return 0;
        // }
        inline int saQuerySeedSetRef(ShiftAnd<MyConst::MISCOUNT>& sa, MATCH::match& mat, uint16_t& qThreshold)
        {

            // use counters to flag what has been processed so far
            std::vector<uint16_t>& threadCountFwdStart = countsFwdStart[omp_get_thread_num()];
            std::vector<uint16_t>& threadCountRevStart = countsRevStart[omp_get_thread_num()];
            auto& fwdMetaIDs_t = fwdMetaIDs[omp_get_thread_num()];
            auto& revMetaIDs_t = revMetaIDs[omp_get_thread_num()];

            // counter for how often we had a match
            std::array<uint8_t, MyConst::MISCOUNT + 1> multiMatch;
            multiMatch.fill(0);

            // will contain matches iff match is found for number of errors specified by index
            std::array<MATCH::match, MyConst::MISCOUNT + 1> uniqueMatches;

            // check all fwd meta CpGs
            for (const auto& m : fwdMetaIDs_t)
            {
                // apply qgram lemma
                if (m.second < qThreshold)
                    continue;

                const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m.first].start];
                const struct CpG& endCpg = ref.cpgTable[ref.metaCpGs[m.first].end];
                auto startIt = ref.fullSeq[startCpg.chrom].begin() + startCpg.pos;
                auto endIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2) + MyConst::MISCOUNT;

                // std::cout << std::string(startIt + 220, endIt) << "\n";
                // check if CpG was too near to the end
                if (endIt > ref.fullSeq[startCpg.chrom].end())
                {
                    // if so move end iterator appropriately
                    endIt = ref.fullSeq[startCpg.chrom].end();
                }

                // use shift and to find all matchings
                std::vector<uint64_t> matchings;
                std::vector<uint8_t> errors;
                sa.querySeq(startIt, endIt, matchings, errors);

                // go through matching and see if we had such a match (with that many errors) before - if so,
                // return to caller reporting no match
                for (size_t i = 0; i < matchings.size(); ++i)
                {

                    // check if we had a match with that many errors before
                    if (multiMatch[errors[i]])
                    {

                        MATCH::match& match_2 = uniqueMatches[errors[i]];
                        const bool isStart = MATCH::isStart(match_2);
                        const bool isFwd = MATCH::isFwd(match_2);
                        // check if same k-mer (borders of meta CpGs)
                        if (isFwd && !isStart && ref.cpgTable[ref.metaCpGs[MATCH::getMetaID(match_2)].start].pos + MATCH::getOffset(match_2) == startCpg.pos + matchings[i])
                        {
                            continue;

                        } else {

                            // check if this is a match without errors
                            if (!errors[i])
                            {

                                // if so, return without a match
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
                        uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i], errors[i], 1, 0, m.first);
                        multiMatch[errors[i]] = 1;
                    }
                }
            }
            // go through reverse sequences
            for (const auto& m : revMetaIDs_t)
            {

                // apply qgram lemma
                if (m.second < qThreshold)
                    continue;

                // retrieve sequence
                const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m.first].start];
                const struct CpG& endCpg = ref.cpgTable[ref.metaCpGs[m.first].end];
                auto endIt = ref.fullSeq[startCpg.chrom].begin() + startCpg.pos - 1;
                auto startIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2) + MyConst::MISCOUNT - 1;

                // check if CpG was too near to the end
                if (startIt >= ref.fullSeq[startCpg.chrom].end())
                {
                    // if so move end iterator appropriately
                    startIt = ref.fullSeq[startCpg.chrom].end() - 1;
                }

                // use shift and to find all matchings
                std::vector<uint64_t> matchings;
                std::vector<uint8_t> errors;
                sa.queryRevSeq(startIt, endIt, matchings, errors);

                // go through matching and see if we had such a match (with that many errors) before - if so,
                // return to caller reporting no match
                for (size_t i = 0; i < matchings.size(); ++i)
                {

                    // check if we had a match with that many errors before
                    if (multiMatch[errors[i]])
                    {

                        MATCH::match& match_2 = uniqueMatches[errors[i]];
                        const bool isStart = MATCH::isStart(match_2);
                        const bool isFwd = MATCH::isFwd(match_2);
                        // check if same k-mer (borders of meta CpGs)
                        if (!isFwd && !isStart && ref.cpgTable[ref.metaCpGs[MATCH::getMetaID(match_2)].start].pos + MATCH::getOffset(match_2) == startCpg.pos + matchings[i])
                        {
                            continue;

                        } else {

                            // check if this is a match without errors
                            if (!errors[i])
                            {

                                // if so, return without a match
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
                        uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i], errors[i], 0, 0, m.first);
                        multiMatch[errors[i]] = 1;
                    }
                }
            }
            // check all fwd meta CpGs that are at start
            for (size_t i = 0; i < threadCountFwdStart.size(); ++i)
            {

                // check if we fulfill the qgram lemma
                // if not - continue with next meta CpG
                if (threadCountFwdStart[i] < qThreshold)
                {
                    continue;
                }
                // retrieve sequence
                const struct CpG& startCpg = ref.cpgStartTable[ref.metaStartCpGs[i].start];
                const struct CpG& endCpg = ref.cpgStartTable[ref.metaStartCpGs[i].end];
                auto startIt = ref.fullSeq[startCpg.chrom].begin();
                auto endIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2) + MyConst::MISCOUNT;

                // check if CpG was too near to the end
                if (endIt > ref.fullSeq[startCpg.chrom].end())
                {
                    // if so move end iterator appropriately
                    endIt = ref.fullSeq[startCpg.chrom].end();
                }

                // use shift and to find all matchings
                std::vector<uint64_t> matchings;
                std::vector<uint8_t> errors;
                sa.querySeq(startIt, endIt, matchings, errors);

                // go through matching and see if we had such a match (with that many errors) before - if so,
                // return to caller reporting no match
                for (size_t j = 0; j < matchings.size(); ++j)
                {

                    // check if we had a match with that many errors before
                    if (multiMatch[errors[j]])
                    {

                        MATCH::match& match_2 = uniqueMatches[errors[j]];
                        const bool isStart = MATCH::isStart(match_2);
                        // check if same k-mer (borders of meta CpGs)
                        if (isStart && MATCH::getOffset(match_2) == matchings[i])
                        {
                            continue;

                        } else {

                            // check if this is a match without errors
                            if (!errors[j])
                            {

                                // if so, return without a match
                                return -1;

                            }
                            // set the number of matches with that many errors to 2
                            // indicating that we do not have a unique match with that many errors
                            multiMatch[errors[j]] = 2;
                        }


                    } else {

                        // we don't have such a match yet,
                        // so save this match at the correct position
                        uniqueMatches[errors[j]] = MATCH::constructMatch(matchings[j], errors[j], 1, 1, i);
                        multiMatch[errors[j]] = 1;
                    }
                }
            }
            // go through reverse sequences of start meta CpGs
            for (size_t i = 0; i < threadCountRevStart.size(); ++i)
            {

                // check if we fulfill the qgram lemma
                // if not - continue with next meta CpG
                if (threadCountRevStart[i] < qThreshold)
                {
                    continue;
                }
                // retrieve sequence
                const struct CpG& startCpg = ref.cpgStartTable[ref.metaStartCpGs[i].start];
                const struct CpG& endCpg = ref.cpgStartTable[ref.metaStartCpGs[i].end];
                auto endIt = ref.fullSeq[startCpg.chrom].begin() - 1;
                auto startIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2) + MyConst::MISCOUNT - 1;

                // check if CpG was too near to the end
                if (startIt >= ref.fullSeq[startCpg.chrom].end())
                {
                    // if so move end iterator appropriately
                    startIt = ref.fullSeq[startCpg.chrom].end() - 1;
                }

                // use shift and to find all matchings
                std::vector<uint64_t> matchings;
                std::vector<uint8_t> errors;
                sa.queryRevSeq(startIt, endIt, matchings, errors);

                // go through matching and see if we had such a match (with that many errors) before - if so,
                // return to caller reporting no match
                for (size_t j = 0; j < matchings.size(); ++j)
                {

                    // check if we had a match with that many errors before
                    if (multiMatch[errors[j]])
                    {

                        MATCH::match& match_2 = uniqueMatches[errors[j]];
                        const bool isStart = MATCH::isStart(match_2);
                        // check if same k-mer (borders of meta CpGs)
                        if (isStart && MATCH::getOffset(match_2) == matchings[i])
                        {
                            continue;

                        } else {

                            // check if this is a match without errors
                            if (!errors[j])
                            {

                                // if so, return without a match
                                return -1;

                            }
                            // set the number of matches with that many errors to 2
                            // indicating that we do not have a unique match with that many errors
                            multiMatch[errors[j]] = 2;
                        }


                    } else {

                        // we don't have such a match yet,
                        // so save this match at the correct position
                        uniqueMatches[errors[j]] = MATCH::constructMatch(matchings[j], errors[j], 0, 1, i);
                        multiMatch[errors[j]] = 1;
                    }
                }
            }



            // go through found matches for each [0,maxErrorNumber] and see if it is unique
            for (size_t i = 0; i < multiMatch.size(); ++i)
            {
                // there is no match with that few errors, search the one with more errors
                if (multiMatch[i] == 0)
                {
                    continue;
                }
                // if match is not unique, return unsuccessfull to caller
                if (multiMatch[i] > 1)
                {
                    // of << "Too bad, multimatch in internal\n";
                    return -1;

                } else {

                    mat = uniqueMatches[i];
                    return 1;
                }

            }
            // we have not a single match at all, return unsuccessfull to caller
            return 0;
        }
        inline void saQuerySeedSetRefPaired(ShiftAnd<MyConst::MISCOUNT>& sa, std::vector<MATCH::match>& mats, uint16_t& qThreshold)
        {

            // use counters to flag what has been processed so far
            std::vector<uint16_t>& threadCountFwdStart = countsFwdStart[omp_get_thread_num()];
            std::vector<uint16_t>& threadCountRevStart = countsRevStart[omp_get_thread_num()];
            auto& fwdMetaIDs_t = fwdMetaIDs[omp_get_thread_num()];
            auto& revMetaIDs_t = revMetaIDs[omp_get_thread_num()];

            // counter for how often we had a match
            // std::array<uint8_t, MyConst::MISCOUNT + 1> multiMatch;
            // multiMatch.fill(0);

            // will contain matches iff match is found for number of errors specified by index
            // std::array<MATCH::match, MyConst::MISCOUNT + 1> uniqueMatches;

            // check all fwd meta CpGs
            for (const auto& m : fwdMetaIDs_t)
            {
                // apply qgram lemma
                if (m.second < qThreshold)
                    continue;

                const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m.first].start];
                const struct CpG& endCpg = ref.cpgTable[ref.metaCpGs[m.first].end];
                auto startIt = ref.fullSeq[startCpg.chrom].begin() + startCpg.pos;
                auto endIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2) + MyConst::MISCOUNT;

                // check if CpG was too near to the end
                if (endIt > ref.fullSeq[startCpg.chrom].end())
                {
                    // if so move end iterator appropriately
                    endIt = ref.fullSeq[startCpg.chrom].end();
                }

                // use shift and to find all matchings
                std::vector<uint64_t> matchings;
                std::vector<uint8_t> errors;
                sa.querySeq(startIt, endIt, matchings, errors);

                // translate found matchings
                for (size_t i = 0; i < matchings.size(); ++i)
                {
                    mats.push_back(std::move(MATCH::constructMatch(matchings[i], errors[i], 1, 0, m.first)));
                }
            }
            // go through reverse sequences
            for (const auto& m : revMetaIDs_t)
            {

                // apply qgram lemma
                if (m.second < qThreshold)
                    continue;

                // retrieve sequence
                const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m.first].start];
                const struct CpG& endCpg = ref.cpgTable[ref.metaCpGs[m.first].end];
                auto endIt = ref.fullSeq[startCpg.chrom].begin() + startCpg.pos - 1;
                auto startIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2) + MyConst::MISCOUNT - 1;

                // check if CpG was too near to the end
                if (startIt >= ref.fullSeq[startCpg.chrom].end())
                {
                    // if so move end iterator appropriately
                    startIt = ref.fullSeq[startCpg.chrom].end() - 1;
                }

                // use shift and to find all matchings
                std::vector<uint64_t> matchings;
                std::vector<uint8_t> errors;
                sa.queryRevSeq(startIt, endIt, matchings, errors);

                for (size_t i = 0; i < matchings.size(); ++i)
                {
                    mats.push_back(std::move(MATCH::constructMatch(matchings[i], errors[i], 0, 0, m.first)));
                }
            }
            // check all fwd meta CpGs that are at start
            for (size_t i = 0; i < threadCountFwdStart.size(); ++i)
            {

                // check if we fulfill the qgram lemma
                // if not - continue with next meta CpG
                if (threadCountFwdStart[i] < qThreshold)
                {
                    continue;
                }
                // retrieve sequence
                const struct CpG& startCpg = ref.cpgStartTable[ref.metaStartCpGs[i].start];
                const struct CpG& endCpg = ref.cpgStartTable[ref.metaStartCpGs[i].end];
                auto startIt = ref.fullSeq[startCpg.chrom].begin();
                auto endIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2) + MyConst::MISCOUNT;

                // check if CpG was too near to the end
                if (endIt > ref.fullSeq[startCpg.chrom].end())
                {
                    // if so move end iterator appropriately
                    endIt = ref.fullSeq[startCpg.chrom].end();
                }

                // use shift and to find all matchings
                std::vector<uint64_t> matchings;
                std::vector<uint8_t> errors;
                sa.querySeq(startIt, endIt, matchings, errors);

                // go through matching and see if we had such a match (with that many errors) before - if so,
                // return to caller reporting no match
                for (size_t j = 0; j < matchings.size(); ++j)
                {
                    mats.push_back(std::move(MATCH::constructMatch(matchings[j], errors[j], 1, 1, i)));
                }
            }
            // go through reverse sequences of start meta CpGs
            for (size_t i = 0; i < threadCountRevStart.size(); ++i)
            {

                // check if we fulfill the qgram lemma
                // if not - continue with next meta CpG
                if (threadCountRevStart[i] < qThreshold)
                {
                    continue;
                }
                // retrieve sequence
                const struct CpG& startCpg = ref.cpgStartTable[ref.metaStartCpGs[i].start];
                const struct CpG& endCpg = ref.cpgStartTable[ref.metaStartCpGs[i].end];
                auto endIt = ref.fullSeq[startCpg.chrom].begin() - 1;
                auto startIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2) + MyConst::MISCOUNT - 1;

                // check if CpG was too near to the end
                if (startIt >= ref.fullSeq[startCpg.chrom].end())
                {
                    // if so move end iterator appropriately
                    startIt = ref.fullSeq[startCpg.chrom].end() - 1;
                }

                // use shift and to find all matchings
                std::vector<uint64_t> matchings;
                std::vector<uint8_t> errors;
                sa.queryRevSeq(startIt, endIt, matchings, errors);

                // go through matching and see if we had such a match (with that many errors) before - if so,
                // return to caller reporting no match
                for (size_t j = 0; j < matchings.size(); ++j)
                {
                    mats.push_back(std::move(MATCH::constructMatch(matchings[j], errors[j], 0, 1, i)));
                }
            }
        }

        // count all metaCpG occurences of k-mers appearing in seq
        //
        // ARGUMENTS:
        //          seq             sequence of the read to query to hash table
        //
        // MODIFICATION:
        //          The threadCount* fields are modified such that they have the count of metaCpGs after
        //          a call to this function
        inline void getSeedRefs(const std::string& seq, const size_t& readSize, const uint16_t qThreshold)
        {

            // std::vector<uint16_t>& threadCountFwd = countsFwd[omp_get_thread_num()];
            // std::vector<uint16_t>& threadCountRev = countsRev[omp_get_thread_num()];
            std::vector<uint16_t>& threadCountFwdStart = countsFwdStart[omp_get_thread_num()];
            std::vector<uint16_t>& threadCountRevStart = countsRevStart[omp_get_thread_num()];
            // fill with zeroes
            // threadCountFwd.assign(ref.metaCpGs.size(), 0);
            // threadCountRev.assign(ref.metaCpGs.size(), 0);
            threadCountFwdStart.assign(ref.metaStartCpGs.size(), 0);
            threadCountRevStart.assign(ref.metaStartCpGs.size(), 0);

            auto& fwdMetaIDs_t = fwdMetaIDs[omp_get_thread_num()];
            auto& revMetaIDs_t = revMetaIDs[omp_get_thread_num()];
            fwdMetaIDs_t.clear();
            revMetaIDs_t.clear();
            // fwdMetaIDs_t.reserve(50000);
            // revMetaIDs_t.reserve(50000);

            // retrieve kmers for first hash
            uint64_t fhVal = ntHash::NTP64(seq.data());

            uint64_t key = fhVal % MyConst::HTABSIZE;

            uint64_t lastId = 0xffffffffffffffffULL;
            bool wasFwd = false;
            bool wasStart = false;

            // maximum position until we can insert completely new meta cpgs
            uint32_t maxQPos = seq.size() - MyConst::KMERLEN - qThreshold;

            for (uint64_t i = ref.tabIndex[key]; i < ref.tabIndex[key+1]; ++i)
            {

                const uint32_t metaId = KMER_S::getMetaCpG(ref.kmerTableSmall[i]);
                const bool isFwd = ref.strandTable[i];
                const bool isStart = KMER_S::isStartCpG(ref.kmerTableSmall[i]);
                // check if we visited meta CpG before
                if (metaId == lastId && isFwd == wasFwd && isStart == wasStart)
                {
                    continue;
                }

                // update vars for last checked metaCpG
                lastId = metaId;
                wasFwd = isFwd;
                wasStart = isStart;
                if (isStart)
                {
                    if (isFwd)
                    {
                        ++threadCountFwdStart[metaId];

                    } else {

                        ++threadCountRevStart[metaId];

                    }

                } else {

                    if (isFwd)
                    {
                        ++fwdMetaIDs_t[metaId];

                    } else {

                        ++revMetaIDs_t[metaId];

                    }
                }
            }

            for (unsigned int cIdx = 0; cIdx < (seq.size() - MyConst::KMERLEN); ++cIdx)
            {

                // use rolling hash
                ntHash::NTP64(fhVal, seq[cIdx], seq[cIdx + MyConst::KMERLEN]);

                key = fhVal % MyConst::HTABSIZE;

                lastId = 0xffffffffffffffffULL;
                wasFwd = false;
                wasStart = false;

                for (uint64_t i = ref.tabIndex[key]; i < ref.tabIndex[key+1]; ++i)
                {

                    const uint32_t metaId = KMER_S::getMetaCpG(ref.kmerTableSmall[i]);
                    const bool isFwd = ref.strandTable[i];
                    const bool isStart = KMER_S::isStartCpG(ref.kmerTableSmall[i]);
                    // check if we visited meta CpG before
                    if (metaId == lastId && isFwd == wasFwd && isStart == wasStart)
                    {
                        continue;
                    }

                    // update vars for last checked metaCpG
                    lastId = metaId;
                    wasFwd = isFwd;
                    wasStart = isStart;
                    if (isStart)
                    {

                        if (isFwd)
                        {
                            ++threadCountFwdStart[metaId];

                        } else {

                            ++threadCountRevStart[metaId];

                        }

                    } else {

                        if (isFwd)
                        {
                            // check if it is at all possible to have newly inserted element passing q
                            if (cIdx < maxQPos)
                            {
                                ++fwdMetaIDs_t[metaId];

                            } else {

                                auto it = fwdMetaIDs_t.find(metaId);
                                if (it != fwdMetaIDs_t.end())
                                {
                                    ++it->second;
                                }
                            }

                        } else {

                            if (cIdx < maxQPos)
                            {
                                ++revMetaIDs_t[metaId];

                            } else {

                                auto it = revMetaIDs_t.find(metaId);
                                if (it != revMetaIDs_t.end())
                                {
                                    ++it->second;
                                }
                            }

                        }
                    }
                }
            }
        }



        // Extract single match for given lists of matches of fwd and reverse complement of a single read
        // Internally updates methylation counts
        //
        // ARGUMENTS:
        //          fwdMatches  array of matches retrieved for original sequence
        //          revMatches  array of matches retrieved for reverse complement of sequence
        //          r           read representation
        //          revSeq      reverse complement sequence of read
        //          *MatchT     Counting objects for matching type
        inline void extractSingleMatch(std::vector<MATCH::match>& fwdMatches, std::vector<MATCH::match>& revMatches, Read& r, std::string& revSeq, uint64_t& succMatchT, uint64_t& unSuccMatchT, uint64_t& nonUniqueMatchT)
        {
            // Construct artificial best match
            MATCH::match bestMat = MATCH::constructMatch(0, MyConst::MISCOUNT + 1,0,0,0);
            bool isUnique = true;
            // Extract best match from list of found matches of fwd reads
            for (MATCH::match mat : fwdMatches)
            {
                if (MATCH::getErrNum(mat) == MATCH::getErrNum(bestMat))
                {
                    // Check if same match
                    uint32_t matPos = ref.cpgTable[ref.metaCpGs[MATCH::getMetaID(mat)].start].pos + MATCH::getOffset(mat);
                    uint32_t bestMatPos = ref.cpgTable[ref.metaCpGs[MATCH::getMetaID(bestMat)].start].pos + MATCH::getOffset(bestMat);

                    // Dealing with large offsets, we need unsigned. Hence check both directions.
                    // check within offset of one for security reasons (insertions deletions etc)
                    if (matPos - bestMatPos > 1 && bestMatPos - matPos > 1)
                    {
                        isUnique = false;
                    }
                } else if (MATCH::getErrNum(mat) < MATCH::getErrNum(bestMat))
                {
                    // update bestMatch
                    bestMat = mat;
                    isUnique = true;
                }

            }
            for (MATCH::match mat : revMatches)
            {
                if (MATCH::getErrNum(mat) == MATCH::getErrNum(bestMat))
                {
                    // Check if same match
                    uint32_t matPos = ref.cpgTable[ref.metaCpGs[MATCH::getMetaID(mat)].start].pos + MATCH::getOffset(mat);
                    uint32_t bestMatPos = ref.cpgTable[ref.metaCpGs[MATCH::getMetaID(bestMat)].start].pos + MATCH::getOffset(bestMat);

                    // Dealing with large offsets, we need unsigned. Hence check both directions.
                    // check within offset of one for security reasons (insertions deletions etc)
                    if (matPos - bestMatPos > 1 && bestMatPos - matPos > 1)
                    {
                        isUnique = false;
                    }
                } else if (MATCH::getErrNum(mat) < MATCH::getErrNum(bestMat))
                {
                    // update bestMatch
                    bestMat = mat;
                    isUnique = true;
                }

            }

            // check if found match is unique and not the artficial initialization
            if (isUnique && MATCH::getErrNum(bestMat) < MyConst::MISCOUNT + 1)
            {
                r.mat = bestMat;
                if (MATCH::isFwd(bestMat))
                {
                    computeMethLvl(bestMat, r.seq);

                } else {

                    computeMethLvl(bestMat, revSeq);
                }
                ++succMatchT;

            } else {

                r.isInvalid = true;
                if (!isUnique)
                {
                    ++nonUniqueMatchT;

                } else if (MATCH::getErrNum(bestMat) < MyConst::MISCOUNT + 1)
                {
                    ++unSuccMatchT;
                }
            }
        }



        // compute the qgram threshold for a given sequence
        // inline uint16_t computeQgramThresh(std::string& seq)
        // {
        //
        //     uint16_t qThresh = seq.size() - MyConst::KMERLEN - (MyConst::KMERLEN * MyConst::MISCOUNT);
        //     of << qThresh << "\t";
        //
        //     uint64_t kSeq = 0;
        //
        //     // read the first k-1 letters
        //     for (unsigned int i = 0; i < MyConst::KMERLEN - 1; ++i)
        //     {
        //         kSeq = kSeq << 2;
        //
        //         switch (seq[i])
        //         {
        //
        //             case 'C':
        //             case 'T':
        //
        //                 kSeq += 3;
        //                 break;
        //
        //             case 'G':
        //
        //                 kSeq += 2;
        //                 break;
        //
        //         }
        //     }
        //
        //     for (unsigned int i = MyConst::KMERLEN - 1; i < seq.size(); ++i)
        //     {
        //
        //         kSeq = kSeq << 2;
        //
        //         switch (seq[i])
        //         {
        //
        //             case 'C':
        //             case 'T':
        //
        //                 kSeq += 3;
        //                 break;
        //
        //             case 'G':
        //
        //                 kSeq += 2;
        //                 break;
        //
        //         }
        //         // test if kmer was filtered -> if yes, reduce qgram lemma threshold
        //         if (ref.filteredKmers.count(kSeq & MyConst::KMERMASK))
        //         {
        //             --qThresh;
        //             of << "okay now it happened\t";
        //         }
        //
        //     }
        //     // if underflow, reset
        //     if (qThresh > seq.size())
        //         qThresh = 0;
        //     of << qThresh << "\n\n";
        //     return qThresh;
        // }


        // compute the methylation levels for the given read by traversing the CpGs of the matched meta CpG
        //
        // ARGUMENTS:
        //          mat     match to process
        //          seq     sequence that was matched (i.e. r.seq or revSeq in main query routine)
        //
        // MODIFICATIONS:
        //              will modify internal methLevel counters
        inline void computeMethLvl(MATCH::match& mat, std::string& seq)
        {

            // retrieve matched metaId
            bool isFwd = MATCH::isFwd(mat);
            bool isStart = MATCH::isStart(mat);
            uint32_t metaID = MATCH::getMetaID(mat);
            uint16_t offset = MATCH::getOffset(mat);
            uint8_t errNum = MATCH::getErrNum(mat);

            // if no errors -> simple lookup of sequences
            if (errNum == 0)
            {

                // retrieve chromosome and position of match
                if (isStart)
                {
                    struct metaCpG& m = ref.metaStartCpGs[metaID];
                    uint8_t chrom = ref.cpgStartTable[m.start].chrom;
                    for (uint32_t cpgId = m.start; cpgId <= m.end; ++cpgId)
                    {
                        // check if CpG is too far downstream of read match
                        // i.e. no overlap
                        if (ref.cpgStartTable[cpgId].pos < offset - seq.size())
                            continue;
                        // check if too far upstream
                        if (ref.cpgStartTable[cpgId].pos + 1 > offset)
                            break;


                        // position of CpG in read
                        // last term represents position of start of read in reference sequence
                        uint32_t readCpGPos = ref.cpgStartTable[cpgId].pos - (offset - seq.size() + 1);
                        if (isFwd)
                        {
                            if (seq[readCpGPos] == 'T')
#ifdef _OPENMP
#pragma omp atomic
#endif
                                ++methLevelsStart[cpgId].unmethFwd;
                            else
#ifdef _OPENMP
#pragma omp atomic
#endif
                                ++methLevelsStart[cpgId].methFwd;
                        } else {

                            if (seq[seq.size() - readCpGPos - 1] == 'T')
#ifdef _OPENMP
#pragma omp atomic
#endif
                                ++methLevelsStart[cpgId].unmethRev;
                            else
#ifdef _OPENMP
#pragma omp atomic
#endif
                                ++methLevelsStart[cpgId].methRev;
                        }
                    }


                // normal match
                } else {

                    struct metaCpG& m = ref.metaCpGs[metaID];
                    uint8_t chrom = ref.cpgTable[m.start].chrom;
                    uint32_t metaPos = ref.cpgTable[m.start].pos;
                    const uint32_t minPos = metaPos + offset - (seq.size() - 1);
                    const uint32_t maxPos = metaPos + offset;
                    for (uint32_t cpgId = m.start; cpgId <= m.end; ++cpgId)
                    {
                        // check if CpG is too far downstream of read match
                        // i.e. no overlap
                        if (ref.cpgTable[cpgId].pos + MyConst::READLEN - 2 < minPos)
                            continue;
                        // check if too far upstream
                        if (isFwd)
                        {
                            if (ref.cpgTable[cpgId].pos + MyConst::READLEN - 2 > maxPos)
                                break;
                        } else {
                            if (ref.cpgTable[cpgId].pos + MyConst::READLEN - 1 > maxPos)
                                break;
                        }



                        // position of CpG in read
                        uint32_t readCpGPos = ref.cpgTable[cpgId].pos + MyConst::READLEN - 2 - (metaPos + offset - (seq.size() - 1));
                        if (isFwd)
                        {
                            if (seq[readCpGPos] == 'T')
                            {
#ifdef _OPENMP
#pragma omp atomic
#endif
                                ++methLevels[cpgId].unmethFwd;
                            }
                            else if (seq[readCpGPos] == 'C')
                            {
#ifdef _OPENMP
#pragma omp atomic
#endif
                                ++methLevels[cpgId].methFwd;
                            }
                            // else std::cout << "This should not happen 1!\n";

                        } else {

                            if (seq[seq.size() - readCpGPos - 2] == 'T')
                            {
#ifdef _OPENMP
#pragma omp atomic
#endif
                                ++methLevels[cpgId].unmethRev;
                            }
                            else  if (seq[seq.size() - readCpGPos - 2] == 'C')
                            {
#ifdef _OPENMP
#pragma omp atomic
#endif
                                ++methLevels[cpgId].methRev;
                            }
                            // else
                            // {
                            //     std::cout << "This should not happen 2!\n";
                            //     std::cout << "Reference sequence (top) and read (bot):\n";
                            //     std::cout << std::string(ref.fullSeq[chrom].data() + metaPos + offset - seq.size() + 1, seq.size()) << "\n";
                            //     std::reverse(seq.begin(), seq.end());
                            //     std::cout << seq << "\n\n";
                            // }
                        }
                    }
                }

            // if match was errornous
            } else {

                if (isStart)
                {
                    struct metaCpG& m = ref.metaStartCpGs[metaID];
                    uint8_t chrom = ref.cpgStartTable[m.start].chrom;

                    const char* refSeq = ref.fullSeq[chrom].data() + offset;

                    // init levenshtein DP algo
                    LevenshtDP<uint16_t, MyConst::MISCOUNT> lev(seq, refSeq);
                    std::vector<ERROR_T> alignment;

                    // minimum position for overlap
                    uint32_t minPos = offset - (seq.size() - 1);
                    uint32_t maxPos = offset - 1;
                    // positions of first/ last overlapping CpG
                    uint32_t minIndex = m.start;
                    uint32_t maxIndex = m.end;
                    for (uint32_t cpgID = m.start; cpgID <= m.end; ++cpgID)
                    {
                        if (ref.cpgStartTable[cpgID].pos < minPos)
                        {
                            ++minIndex;
                        } else if (ref.cpgStartTable[cpgID].pos > maxPos)
                        {
                            maxIndex = cpgID - 1;
                            break;
                        }
                    }

                    if (isFwd)
                    {

                        // compute alignment
                        lev.runDPFill<CompiFwd>(cmpFwd);
                        lev.backtrackDP<CompiFwd>(cmpFwd, alignment);
                        uint32_t refSeqPos = offset;
                        uint32_t readSeqPos = seq.size() - 1;
                        uint32_t alignPos = alignment.size() - 1;

                        // go through all overlapping CpGs from back,
                        // move through read and reference according to alignment
                        // if position of CpG is hit, compare and count
                        for (uint32_t cpgID = maxIndex; cpgID >= minIndex; --cpgID)
                        {
                            // align until this CpG
                            while (ref.cpgStartTable[cpgID].pos < refSeqPos)
                            {
                                switch (alignment[alignPos--])
                                {
                                    case (MATCHING):
                                    case (MISMATCH):
                                        --readSeqPos;
                                        --refSeqPos;
                                        break;
                                    case (DELETION):
                                        --refSeqPos;
                                        break;
                                    case(INSERTION):
                                        --readSeqPos;
                                        break;
                                }
                                // check if we have a CpG aligned to the reference CpG
                                if (seq[readSeqPos + 1] == 'G')
                                {
                                    // check for methylated C
                                    if (seq[readSeqPos] == 'C')
                                    {
#ifdef _OPENMP
#pragma omp atomic
#endif
                                        ++methLevelsStart[cpgID].methFwd;
                                    }
                                    else if (seq[readSeqPos] == 'T')
                                    {
#ifdef _OPENMP
#pragma omp atomic
#endif
                                        ++methLevelsStart[cpgID].unmethFwd;
                                    }

                                }
                            }
                        }

                    } else {

                        lev.runDPFillRev<CompiRev>(cmpRev);
                        lev.backtrackDPRev<CompiRev>(cmpRev, alignment);
                        uint32_t refSeqPos = offset;
                        uint32_t readSeqPos = seq.size() - 1;
                        uint32_t alignPos = alignment.size() - 1;

                        // go through all overlapping CpGs from back,
                        // move through read and reference according to alignment
                        // if position of CpG is hit, compare and count
                        for (uint32_t cpgID = maxIndex; cpgID >= minIndex; --cpgID)
                        {
                            // align until this CpG
                            while (ref.cpgStartTable[cpgID].pos < refSeqPos)
                            {
                                switch (alignment[alignPos--])
                                {
                                    case (MATCHING):
                                    case (MISMATCH):
                                        --readSeqPos;
                                        --refSeqPos;
                                        break;
                                    case (DELETION):
                                        --refSeqPos;
                                        break;
                                    case(INSERTION):
                                        --readSeqPos;
                                        break;
                                }
                                // check if we have a CpG aligned to the reference CpG
                                if (seq[readSeqPos] == 'G')
                                {
                                    // check for methylated C (on reverse!)
                                    if (seq[readSeqPos + 1] == 'C')
                                    {
#ifdef _OPENMP
#pragma omp atomic
#endif
                                        ++methLevelsStart[cpgID].methRev;
                                    }
                                    else if (seq[readSeqPos] == 'T')
                                    {
#ifdef _OPENMP
#pragma omp atomic
#endif
                                        ++methLevelsStart[cpgID].unmethRev;
                                    }

                                }
                            }
                        }
                    }


                // not start CpG (normal)
                } else {

                    struct metaCpG& m = ref.metaCpGs[metaID];
                    uint8_t chrom = ref.cpgTable[m.start].chrom;
                    uint32_t metaPos = ref.cpgTable[m.start].pos;

                    const char* refSeq = ref.fullSeq[chrom].data() + metaPos + offset;

                    // init levenshtein DP algo
                    LevenshtDP<uint16_t, MyConst::MISCOUNT> lev(seq, refSeq);
                    std::vector<ERROR_T> alignment;

                    // minimum position for overlap
                    const uint32_t minPos = metaPos + offset - (seq.size() - 1);
                    const uint32_t maxPos = metaPos + offset;
                    // positions of first/ last overlapping CpG
                    int32_t minIndex = m.start;
                    int32_t maxIndex = m.end;
                    for (int32_t cpgID = m.start; cpgID <= m.end; ++cpgID)
                    {
                        if (ref.cpgTable[cpgID].pos + MyConst::READLEN - 2 < minPos)
                        {
                            ++minIndex;

                        } else if (isFwd)
                        {
                            if (ref.cpgTable[cpgID].pos + MyConst::READLEN - 2 > maxPos)
                            {
                                maxIndex = cpgID - 1;
                                break;
                            }
                        } else {
                            if (ref.cpgTable[cpgID].pos + MyConst::READLEN - 1 > maxPos)
                            {
                                maxIndex = cpgID - 1;
                                break;
                            }
                        }
                    }

                    if (isFwd)
                    {

                        // compute alignment
                        lev.runDPFill<CompiFwd>(cmpFwd);
                        lev.backtrackDP<CompiFwd>(cmpFwd, alignment);
                        // current position in read and reference
                        uint32_t refSeqPos = metaPos + offset;
                        int32_t readSeqPos = seq.size() - 1;
                        // current position in alignment (note that we align from right to left)
                        int32_t alignPos = alignment.size() - 1;
                        // sanity check
                        // if (lev.getEditDist() != errNum)
                        // {
                        //     // Report error and print out found alignment
                        //     std::cout << "Editdist: " << lev.getEditDist() << " shiftand: " << errNum << "\n";
                            // std::string readAl (alignment.size(),'+');
                            // std::string refAl (alignment.size(), '+');
                            // for (auto rIt = alignment.rbegin(); rIt != alignment.rend(); ++rIt, --alignPos)
                            // {
                            //     switch (*rIt)
                            //     {
                            //         case (MATCHING):
                            //         case (MISMATCH):
                            //             readAl[alignPos] = seq[readSeqPos];
                            //             refAl[alignPos] = ref.fullSeq[chrom][refSeqPos];
                            //             --readSeqPos;
                            //             --refSeqPos;
                            //             break;
                            //         case (DELETION):
                            //             readAl[alignPos] = '-';
                            //             refAl[alignPos] = ref.fullSeq[chrom][refSeqPos];
                            //             --refSeqPos;
                            //             break;
                            //         case(INSERTION):
                            //             readAl[alignPos] = seq[readSeqPos];
                            //             refAl[alignPos] = '-';
                            //             --readSeqPos;
                            //             break;
                            //     }
                            // }
                            // std::cout << "Alignment seems to be wrong! (Read top, reference bottom)\n" << readAl << "\n" << refAl << "\n\n";
                            // exit(1);
                        // }

                        // go through all overlapping CpGs from back,
                        // move through read and reference according to alignment
                        // if position of CpG is hit, compare and count
                        for (int32_t cpgID = maxIndex; cpgID >= minIndex; --cpgID)
                        {
                            // align until this CpG
                            while (ref.cpgTable[cpgID].pos + MyConst::READLEN - 2 < refSeqPos && alignPos >= 0)
                            {
                                switch (alignment[alignPos])
                                {
                                    case (MATCHING):
                                    case (MISMATCH):
                                        --readSeqPos;
                                        --refSeqPos;
                                        break;
                                    case (DELETION):
                                        --refSeqPos;
                                        break;
                                    case(INSERTION):
                                        --readSeqPos;
                                        break;
                                }
                                if (readSeqPos < 0)
                                    break;
                                if (readSeqPos == seq.size() - 1)
                                    continue;
                                --alignPos;
                            }
                            if (readSeqPos < 0)
                            {
                                break;
                            }
                            // check if we have a CpG aligned to the reference CpG
                            // if (seq[readSeqPos + 1] == 'G')
                            // {
                                // check for unmethylated C
                                if (seq[readSeqPos] == 'C')
                                {
#ifdef _OPENMP
#pragma omp atomic
#endif
                                    ++methLevels[cpgID].methFwd;
                                }
                                else if (seq[readSeqPos] == 'T')
                                {
#ifdef _OPENMP
#pragma omp atomic
#endif
                                    ++methLevels[cpgID].unmethFwd;
                                }

                            // }
                        }

                    } else {

                        lev.runDPFillRev<CompiRev>(cmpRev);
                        lev.backtrackDPRev<CompiRev>(cmpRev, alignment);
                        uint32_t refSeqPos = metaPos + offset;
                        int32_t readSeqPos = seq.size() - 1;
                        int32_t alignPos = alignment.size() - 1;
                        std::reverse(seq.begin(),seq.end());
                        // sanity check
                        // if (lev.getEditDist() != errNum)
                        // {
                        //     // Report error and print out found alignment
                        //     std::cout << "Editdist: " << lev.getEditDist() << " shiftand: " << static_cast<uint16_t>(errNum) << "\n";
                            // std::string readAl (alignment.size(),'+');
                            // std::string refAl (alignment.size(), '+');
                            // for (auto rIt = alignment.rbegin(); rIt != alignment.rend(); ++rIt, --alignPos)
                            // {
                            //     switch (*rIt)
                            //     {
                            //         case (MATCHING):
                            //         case (MISMATCH):
                            //             readAl[alignPos] = seq[readSeqPos];
                            //             refAl[alignPos] = ref.fullSeq[chrom][refSeqPos];
                            //             --readSeqPos;
                            //             --refSeqPos;
                            //             break;
                            //         case (DELETION):
                            //             readAl[alignPos] = '-';
                            //             refAl[alignPos] = ref.fullSeq[chrom][refSeqPos];
                            //             --refSeqPos;
                            //             break;
                            //         case(INSERTION):
                            //             readAl[alignPos] = seq[readSeqPos];
                            //             refAl[alignPos] = '-';
                            //             --readSeqPos;
                            //             break;
                            //     }
                            // }
                            // std::cout << "Alignment seems to be wrong! (Read top, reference bottom)\n" << readAl << "\n" << refAl << "\n";
                            // std::cout << "Full reference: " << std::string(ref.fullSeq[chrom].begin() + metaPos + offset - 102, ref.fullSeq[chrom].begin() + metaPos + offset) << "\n\n";
                        //     exit(1);
                        // }

                        // go through all overlapping CpGs from back,
                        // move through read and reference according to alignment
                        // if position of CpG is hit, compare and count
                        for (int32_t cpgID = maxIndex; cpgID >= minIndex; --cpgID)
                        {
                            // align until this CpG
                            while (ref.cpgTable[cpgID].pos + MyConst::READLEN - 2 < refSeqPos && alignPos >= 0)
                            {
                                switch (alignment[alignPos])
                                {
                                    case (MATCHING):
                                    case (MISMATCH):
                                        --readSeqPos;
                                        --refSeqPos;
                                        break;
                                    case (DELETION):
                                        --refSeqPos;
                                        break;
                                    case(INSERTION):
                                        --readSeqPos;
                                        break;
                                }
                                if (readSeqPos < 0)
                                    break;
                                if (readSeqPos == seq.size() - 1)
                                    continue;
                                --alignPos;
                            }
                            // check if we have a CpG aligned to the reference CpG
                            // TODO
                            // if (seq[readSeqPos] == 'G')
                            // {
                                // check for unmethylated C
                                if (seq[readSeqPos + 1] == 'C')
                                {
#ifdef _OPENMP
#pragma omp atomic
#endif
                                    ++methLevels[cpgID].methRev;
                                }
                                else if (seq[readSeqPos + 1] == 'T')
                                {
#ifdef _OPENMP
#pragma omp atomic
#endif
                                    ++methLevels[cpgID].unmethRev;
                                }

                            // }
                        }
                    }
                }
            }
        }

        // print statistics over seed set to statFile and countFile
        //
        // statFile contains statistics over how many times (at most n) the same meta CpG appears in the seed list of one kmer
        // blocks of 4 lines show the change before countfilter, after countfilter, after bitmatch, after second countfilter
        //
        // OUTPUT FORMAT (tsv):
        // #occurences of same meta cpg   1   2   3   ...    n
        // read 1 firstLayer
        // read 1 secondLayer
        // read 1 thirdlLayer
        // read 1 fourthLayer
        // read 2 firstLayer
        // read 2 ...
        //
        //
        // countFile contains counts on how many seeds were found for the given read, 4 columns forming the layers
        // before countfilter, after countfilter, after bitmatch, after second countfilter
        //
        // OUTPUT FORMAT (tsv):
        // Layer    1   2   3   4
        // read1
        // read2
        // read3
        // ...
        //
        // function called for each layer
        //
        // ARGUMENTS:
        //              SeedsK      (current) seed set
        //
        // RETURN:
        //              void
        //
        // MODIFICATIONS:
        //              none
        //
        // TODO
        // static constexpr unsigned int n = 400;
        // std::ofstream statFile;
        // std::ofstream countFile;
        // void printStatistics(const std::vector<std::vector<KMER_S::kmer> > SeedsK);

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
        std::array<uint8_t, 256> lmap;

        // holds counts for each thread for counting heuristic
        // for forward and reverse strand metaCpGs, respectively
        std::array<std::vector<uint16_t>, CORENUM> countsFwd;
        std::array<std::vector<uint16_t>, CORENUM> countsRev;
        std::array<std::vector<uint16_t>, CORENUM> countsFwdStart;
        std::array<std::vector<uint16_t>, CORENUM> countsRevStart;
        // std::array<std::unordered_map<uint32_t, uint16_t, MetaHash>, CORENUM> fwdMetaIDs;
        // std::array<std::unordered_map<uint32_t, uint16_t, MetaHash>, CORENUM> revMetaIDs;
        // std::array<spp::sparse_hash_map<uint32_t, uint16_t, MetaHash>, CORENUM> fwdMetaIDs;
        // std::array<spp::sparse_hash_map<uint32_t, uint16_t, MetaHash>, CORENUM> revMetaIDs;
        std::array<google::dense_hash_map<uint32_t, uint16_t, MetaHash>, CORENUM> fwdMetaIDs;
        std::array<google::dense_hash_map<uint32_t, uint16_t, MetaHash>, CORENUM> revMetaIDs;

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


        // TODO
        std::ofstream of;

};

#endif /* READQUEUE_H */
