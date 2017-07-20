#ifndef READQUEUE_H
#define READQUEUE_H

#include <string>
#include <fstream>
#include <unordered_map>
#include <array>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "gzstream/gzstream.h"
#include "CONST.h"
#include "RefGenome.h"
#include "Read.h"
#include "ShiftAnd.h"

class ReadQueue
{

    public:

        // Ctor -----


        ReadQueue() = delete;

        ReadQueue(const char* filePath, RefGenome& ref, bool isGZ);

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
        // Extend remaining seeds with BitMatch(...)
        bool matchReads(const unsigned int& procReads);


    private:

        // filters seeds according to simple counting criteria
        // #kmers of one metaCpG should be > READLEN - KMERLEN + 1 - (KMERLEN * MISCOUNT)
        inline void filterHeuSeeds(std::vector<std::vector<KMER::kmer> >& seedsK, std::vector<std::vector<bool> >& seedsS, const unsigned int readSize)
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

                    const uint64_t metaId = KMER::getMetaCpG(seedsK[i][j]);
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
                        if (threadCountFwd[KMER::getMetaCpG(seedsK[i][j])] >= countCut)
                        {

                            *filterItK = *srcItK;
                            *filterItS = *srcItS;
                            ++filterItK;
                            ++filterItS;

                        }

                    } else {

                        // test for strict heuristic criterias
                        if (threadCountRev[KMER::getMetaCpG(seedsK[i][j])] >= countCut)
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
        //             KMER::kmer& k = ref.kmerTable[j];
        //             const bool isFwd = ref.strandTable[j];
        //             const uint64_t metaId = KMER::getMetaCpG(k);
        //             const bool isStart = KMER::isStartCpG(k);
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
        inline void bitMatching(const Read& r, std::vector<std::vector<KMER::kmer> >& seedsK, std::vector<std::vector<bool> >& seedsS)
        {

            // masking for actual kmer bits
            constexpr uint64_t signiBits = 0xffffffffffffffffULL >> (64 - (2*MyConst::KMERLEN));
            // bit representation of current kmer of read
            uint64_t kmerBits = 0;
            // generate bit representation of first kmerlen - 1 letters of read
            for (unsigned int i = 0; i < (MyConst::KMERLEN - 1); ++i)
            {

                kmerBits = kmerBits << 2;
                kmerBits |= BitFun::getBitRepr(r.seq[i]);

            }

            // Note that seedsK must have the same size as seedsS anyway, this way we may have a cache hit
            std::vector<std::vector<KMER::kmer> > newSeedsK(seedsK.size());
            std::vector<std::vector<bool> > newSeedsS(seedsK.size());

            // go over each read kmer and compare with reference seeds
            for (unsigned int offset = 0; offset < (r.seq.size() - MyConst::KMERLEN + 1); ++offset)
            {

                // retrieve seeds for current kmer
                std::vector<KMER::kmer>& localSeedsK = seedsK[offset];
                std::vector<bool>& localSeedsS = seedsS[offset];
                // reserve some space for result seedlist
                newSeedsK[offset].reserve(localSeedsK.size());
                newSeedsS[offset].reserve(localSeedsK.size());

                // update current read kmer representation
                kmerBits = kmerBits << 2;
                kmerBits = (kmerBits | BitFun::getBitRepr(r.seq[offset + MyConst::KMERLEN - 1])) & signiBits;

                // iterate over corresponding seeds for reference
                for (unsigned int i = 0; i < localSeedsK.size(); ++i)
                {

                    KMER::kmer& refKmer = localSeedsK[i];

                    // will hold the first CpG in metaCpG after retrieving the meta CpG info
                    uint8_t chrom;
                    // will hold the position of the meta CpG in the genome (chromosomal region specified by chrom of the genome)
                    uint32_t pos = 0;
                    // retrieve reference bit representation
                    //
                    // First retrieve meta CpG info
                    if (KMER::isStartCpG(refKmer))
                    {

                        const uint32_t& cpgInd = ref.metaStartCpGs[KMER::getMetaCpG(refKmer)].start;
                        chrom = ref.cpgStartTable[cpgInd].chrom;

                    } else {

                        const uint32_t& cpgInd = ref.metaCpGs[KMER::getMetaCpG(refKmer)].start;
                        chrom = ref.cpgTable[cpgInd].chrom;
                        pos = ref.cpgTable[cpgInd].pos;

                    }
                    // will hold bit representation of seed
                    uint64_t refKmerBit;
                    // decide if forward or reverse strand of reference
                    //
                    // if forward
                    if (localSeedsS[i])
                    {
                        // retrieve sequence in forward strand
                        refKmerBit = ref.genomeBit[chrom].getSeqKmer(pos + KMER::getOffset(refKmer));

                    // is reverse
                    } else {

                        // retrieve sequence in forward strand
                        refKmerBit = ref.genomeBit[chrom].getSeqKmerRev(pos + KMER::getOffset(refKmer));
                    }

                    // COMPARE read kmer and seed
                    //  matching is 0 iff full match
                    if ( !( refKmerBit ^ (kmerBits & BitFun::getMask(refKmerBit)) ) )
                    {

                        // if we have a match, keep this kmer and strand flag in list
                        newSeedsK[offset].emplace_back(refKmer);
                        newSeedsS[offset].emplace_back(localSeedsS[i]);
                    }
                }
                newSeedsK[offset].shrink_to_fit();
                newSeedsS[offset].shrink_to_fit();

            }

            seedsK = std::move(newSeedsK);
            seedsS = std::move(newSeedsS);

        }
        inline void bitMatchingRev(const Read& r, std::vector<std::vector<KMER::kmer> >& seedsK, std::vector<std::vector<bool> >& seedsS)
        {

            const unsigned int readSize = r.seq.size();
            // masking for actual kmer bits
            constexpr uint64_t signiBits = 0xffffffffffffffffULL >> (64 - (2*MyConst::KMERLEN));
            // bit representation of current kmer of read
            uint64_t kmerBits = 0;
            // generate bit representation of first kmerlen - 1 letters of reverse complement of read
            // Not that we start reading from right
            for (unsigned int i = readSize - 1; i > readSize - MyConst::KMERLEN; --i)
            {

                kmerBits = kmerBits << 2;
                kmerBits |= BitFun::getBitReprRev(r.seq[i]);

            }

            // Note that seedsK must have the same size as seedsS anyway, this way we may have a cache hit
            std::vector<std::vector<KMER::kmer> > newSeedsK(seedsK.size());
            std::vector<std::vector<bool> > newSeedsS(seedsK.size());

            // index for seed vector for current read kmer
            unsigned int kmerInd = 0;
            // go over each read kmer and compare with reference seeds
            for (unsigned int offset = readSize - MyConst::KMERLEN; offset > 0; --offset, ++kmerInd)
            {

                std::vector<KMER::kmer>& localSeedsK = seedsK[kmerInd];
                std::vector<bool>& localSeedsS = seedsS[kmerInd];
                // reserve some space for seedlist
                newSeedsK[kmerInd].reserve(localSeedsK.size());
                newSeedsS[kmerInd].reserve(localSeedsK.size());

                // update current read kmer representation
                kmerBits = kmerBits << 2;
                kmerBits = (kmerBits | BitFun::getBitReprRev(r.seq[offset - 1])) & signiBits;

                // iterate over corresponding seeds for reference
                for (unsigned int i = 0; i < localSeedsK.size(); ++i)
                {

                    KMER::kmer& refKmer = localSeedsK[i];

                    // will hold the chromosome index in metaCpG after retrieving the meta CpG info
                    uint8_t chrom;
                    // will hold the position of the meta CpG in the genome (chromosomal region specified by chrom of the genome)
                    uint32_t pos = 0;
                    // retrieve reference bit representation
                    //
                    // First retrieve meta CpG info
                    if (KMER::isStartCpG(refKmer))
                    {

                        uint32_t& cpgInd = ref.metaStartCpGs[KMER::getMetaCpG(refKmer)].start;
                        chrom = ref.cpgStartTable[cpgInd].chrom;

                    } else {

                        uint32_t& cpgInd = ref.metaCpGs[KMER::getMetaCpG(refKmer)].start;
                        chrom = ref.cpgTable[cpgInd].chrom;
                        pos = ref.cpgTable[cpgInd].pos;

                    }
                    // will hold bit representation of seed
                    uint64_t refKmerBit;
                    // decide if forward or reverse strand of reference
                    //
                    // if forward
                    if (localSeedsS[i])
                    {
                        // retrieve sequence in forward strand
                        refKmerBit = ref.genomeBit[chrom].getSeqKmer(pos + KMER::getOffset(refKmer));

                    // is reverse
                    } else {

                        // retrieve sequence in forward strand
                        refKmerBit = ref.genomeBit[chrom].getSeqKmerRev(pos + KMER::getOffset(refKmer));
                    }

                    // COMPARE read kmer and seed
                    //  matching is 0 iff full match
                    if ( !( refKmerBit ^ (kmerBits & BitFun::getMask(refKmerBit)) ) )
                    {

                        // if we have a match, keep this kmer and strand flag in list
                        newSeedsK[kmerInd].emplace_back(refKmer);
                        newSeedsS[kmerInd].emplace_back(localSeedsS[i]);
                    }
                }
                newSeedsK[kmerInd].shrink_to_fit();
                newSeedsS[kmerInd].shrink_to_fit();

            }

            seedsK = std::move(newSeedsK);
            seedsS = std::move(newSeedsS);
        }


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
        inline int saQuerySeedSet(ShiftAnd<MyConst::MISCOUNT>& sa, std::vector<std::vector<KMER::kmer> >& seedsK, std::vector<std::vector<bool> >& seedsS, MATCH::match& mat)
        {

            // use counters to flag what has been processed so far
            // 0 if not processed at all
            // 1 if start metaCpG has been processed before
            // 2 if metaCpG with this index has been proc. before
            // 3 if both
            std::vector<uint16_t>& threadCountFwd = countsFwd[omp_get_thread_num()];
            std::vector<uint16_t>& threadCountRev = countsRev[omp_get_thread_num()];

            // fill with zeroes
            threadCountFwd.assign(ref.metaCpGs.size(), 0);
            threadCountRev.assign(ref.metaCpGs.size(), 0);

            // counter for how often we had a match
            std::array<uint8_t, MyConst::MISCOUNT + 1> multiMatch;
            multiMatch.fill(0);

            // will contain matches iff match is found for number of errors specified by index
            std::array<MATCH::match, MyConst::MISCOUNT + 1> uniqueMatches;

            for (size_t outerNdx = 0; outerNdx < seedsK.size(); ++outerNdx)
            {

                for (size_t innerNdx = 0; innerNdx < seedsK[outerNdx].size(); ++innerNdx)
                {

                    // retrieve kmer
                    KMER::kmer& k = seedsK[outerNdx][innerNdx];

                    // retrieve meta CpG
                    const uint64_t m = KMER::getMetaCpG(k);

                    // retrieve if start meta CpG
                    const bool isStart = KMER::isStartCpG(k);

                    // retrieve strand
                    const bool isFwd = seedsS[outerNdx][innerNdx];

                    if (isStart)
                    {

                        if (isFwd)
                        {
                            // check if we queried this meta CpG it before
                            if (threadCountFwd[m] == 1 || threadCountFwd[m] >= 3)
                            {
                                continue;

                            } else {

                                // state that we queried this
                                threadCountFwd[m] += 1;
                                // retrieve it
                                const struct CpG& startCpg = ref.cpgStartTable[ref.metaStartCpGs[m].start];
                                const struct CpG& endCpg = ref.cpgStartTable[ref.metaStartCpGs[m].end];

                                auto startIt = ref.fullSeq[startCpg.chrom].begin();
                                auto endIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2);

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

                                        MATCH::match& m_2 = uniqueMatches[errors[i]];
                                        // check if same k-mer (borders of meta CpGs)
                                        if ((MATCH::getChrom(m_2) == startCpg.chrom) && (MATCH::getOffset(m_2) == matchings[i]) && (MATCH::isFwd(m_2)))
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

                                        // we don't have such a match yet,
                                        // so save this match at the correct position
                                        uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i], startCpg.chrom, errors[i], 1);
                                        multiMatch[errors[i]] = 1;
                                    }


                                }

                            }

                        // is not forward strand meta CpG
                        } else {

                            // check if we queried it before
                            if (threadCountRev[m] == 1 || threadCountRev[m] >= 3)
                            {
                                continue;

                            } else {

                                threadCountRev[m] += 1;
                                // retrieve it
                                const struct CpG& startCpg = ref.cpgStartTable[ref.metaStartCpGs[m].start];
                                const struct CpG& endCpg = ref.cpgStartTable[ref.metaStartCpGs[m].end];

                                auto endIt = ref.fullSeq[startCpg.chrom].begin() - 1;
                                auto startIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2) - 1;

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

                                        MATCH::match& m_2 = uniqueMatches[errors[i]];
                                        // check if same k-mer (borders of meta CpGs)
                                        if ((MATCH::getChrom(m_2) == startCpg.chrom) && (MATCH::getOffset(m_2) == matchings[i]) && !(MATCH::isFwd(m_2)))
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

                                        // we don't have such a match yet,
                                        // so save this match at the correct position
                                        uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i], startCpg.chrom, errors[i], 0);
                                        multiMatch[errors[i]] = 1;
                                    }

                                }

                            }

                        }

                    // if is not start metaCpG
                    } else {

                        if (isFwd)
                        {

                            // check if we queried this meta CpG it before
                            if (threadCountFwd[m] >= 2)
                            {
                                continue;

                            } else {


                                // state that we queried this
                                threadCountFwd[m] += 2;
                                // retrieve it
                                const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m].start];
                                const struct CpG& endCpg = ref.cpgTable[ref.metaCpGs[m].end];

                                auto startIt = ref.fullSeq[startCpg.chrom].begin() + startCpg.pos;
                                auto endIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2);

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

                                        MATCH::match& m_2 = uniqueMatches[errors[i]];
                                        // check if same k-mer (borders of meta CpGs)
                                        if ((MATCH::getChrom(m_2) == startCpg.chrom) && (MATCH::getOffset(m_2) == matchings[i] + startCpg.pos) && (MATCH::isFwd(m_2)))
                                        {
                                            continue;

                                        } else {

                                            // check if this is a match without errors
                                            if (!errors[i])
                                            {

                                                // if so, return without a match
                                                // std::cout << "Nonunique no-error match\n";
                                                return -1;

                                            }
                                            // set the number of matches with that many errors to 2
                                            // indicating that we do not have a unique match with that many errors
                                            multiMatch[errors[i]] = 2;
                                        }


                                    } else {

                                        // we don't have such a match yet,
                                        // so save this match at the correct position
                                        uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i] + startCpg.pos, startCpg.chrom, errors[i], 1);
                                        multiMatch[errors[i]] = 1;
                                    }
                                }
                            }

                        // kmer is on backward strand
                        } else {

                            // check if we queried this meta CpG it before
                            if (threadCountRev[m] >= 2)
                            {
                                continue;

                            } else {

                                // state that we queried this
                                threadCountRev[m] += 2;
                                // retrieve it
                                const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[m].start];
                                const struct CpG& endCpg = ref.cpgTable[ref.metaCpGs[m].end];

                                auto endIt = ref.fullSeq[startCpg.chrom].begin() + startCpg.pos - 1;
                                auto startIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2) - 1;

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

                                        MATCH::match& m_2 = uniqueMatches[errors[i]];
                                        if ((MATCH::getChrom(m_2) == startCpg.chrom) && (MATCH::getOffset(m_2) == matchings[i] + startCpg.pos) && !(MATCH::isFwd(m_2)))
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

                                        // we don't have such a match yet,
                                        // so save this match at the correct position
                                        uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i] + startCpg.pos, startCpg.chrom, errors[i], 0);
                                        multiMatch[errors[i]] = 1;
                                    }
                                }

                            } // end else for visited meta CpG test
                        } // end else isFwd
                    } // end else isStart
                } // end inner kmer loop
            } // end outer kmer loop


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
            // of << "No match at all\t";
            return 0;
        }
        inline int saQuerySeedSetRef(ShiftAnd<MyConst::MISCOUNT>& sa, MATCH::match& mat)
        {

            // use counters to flag what has been processed so far
            std::vector<uint16_t>& threadCountFwd = countsFwd[omp_get_thread_num()];
            std::vector<uint16_t>& threadCountRev = countsRev[omp_get_thread_num()];
            std::vector<uint16_t>& threadCountFwdStart = countsFwdStart[omp_get_thread_num()];
            std::vector<uint16_t>& threadCountRevStart = countsRevStart[omp_get_thread_num()];

            // counter for how often we had a match
            std::array<uint8_t, MyConst::MISCOUNT + 1> multiMatch;
            multiMatch.fill(0);

            // will contain matches iff match is found for number of errors specified by index
            std::array<MATCH::match, MyConst::MISCOUNT + 1> uniqueMatches;

            // uint16_t qThreshold = sa.size() - MyConst::KMERLEN - (MyConst::KMERLEN * MyConst::MISCOUNT);
            uint16_t qThreshold = 1;
            // check for overflow (i.e. read is to small for lemma)
            if (qThreshold > sa.size())
            {
                qThreshold = 0;
            }

            // check all fwd meta CpGs
            for (size_t i = 0; i < threadCountFwd.size(); ++i)
            {

                // check if we fulfill the qgram lemma
                // if not - continue with next meta CpG
                if (threadCountFwd[i] < qThreshold)
                {
                    continue;
                }
                // retrieve sequence
                const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[i].start];
                const struct CpG& endCpg = ref.cpgTable[ref.metaCpGs[i].end];
                auto startIt = ref.fullSeq[startCpg.chrom].begin() + startCpg.pos;
                auto endIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2);

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

                        MATCH::match& m_2 = uniqueMatches[errors[i]];
                        // check if same k-mer (borders of meta CpGs)
                        if ((MATCH::getChrom(m_2) == startCpg.chrom) && (MATCH::getOffset(m_2) == matchings[i]) && (MATCH::isFwd(m_2)))
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

                        // we don't have such a match yet,
                        // so save this match at the correct position
                        uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i] + startCpg.pos, startCpg.chrom, errors[i], 1);
                        multiMatch[errors[i]] = 1;
                    }
                }
            }
            // go through reverse sequences
            for (size_t i = 0; i < threadCountRev.size(); ++i)
            {

                // check if we fulfill the qgram lemma
                // if not - continue with next meta CpG
                if (threadCountRev[i] < qThreshold)
                {
                    continue;
                }
                // retrieve sequence
                const struct CpG& startCpg = ref.cpgTable[ref.metaCpGs[i].start];
                const struct CpG& endCpg = ref.cpgTable[ref.metaCpGs[i].end];
                auto endIt = ref.fullSeq[startCpg.chrom].begin() + startCpg.pos - 1;
                auto startIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2) - 1;

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

                        MATCH::match& m_2 = uniqueMatches[errors[i]];
                        // check if same k-mer (borders of meta CpGs)
                        if ((MATCH::getChrom(m_2) == startCpg.chrom) && (MATCH::getOffset(m_2) == matchings[i]) && !(MATCH::isFwd(m_2)))
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

                        // we don't have such a match yet,
                        // so save this match at the correct position
                        uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i] + startCpg.pos, startCpg.chrom, errors[i], 0);
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
                auto endIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2);

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

                        MATCH::match& m_2 = uniqueMatches[errors[i]];
                        // check if same k-mer (borders of meta CpGs)
                        if ((MATCH::getChrom(m_2) == startCpg.chrom) && (MATCH::getOffset(m_2) == matchings[i]) && (MATCH::isFwd(m_2)))
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

                        // we don't have such a match yet,
                        // so save this match at the correct position
                        uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i], startCpg.chrom, errors[i], 1);
                        multiMatch[errors[i]] = 1;
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
                auto startIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + (2*MyConst::READLEN - 2) - 1;

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

                        MATCH::match& m_2 = uniqueMatches[errors[i]];
                        // check if same k-mer (borders of meta CpGs)
                        if ((MATCH::getChrom(m_2) == startCpg.chrom) && (MATCH::getOffset(m_2) == matchings[i]) && !(MATCH::isFwd(m_2)))
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

                        // we don't have such a match yet,
                        // so save this match at the correct position
                        uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i], startCpg.chrom, errors[i], 0);
                        multiMatch[errors[i]] = 1;
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
            // of << "No match at all\t";
            return 0;
        }

        // count all metaCpG occurences of k-mers appearing in seq
        //
        // ARGUMENTS:
        //          seq             sequence of the read to query to hash table
        //
        // RETURN:
        //          false iff too many matchable positions
        //
        // MODIFICATION:
        //          The threadCount* fields are modified such that they have the count of metaCpGs after
        //          a call to this function
        inline bool getSeedRefs(std::vector<char>& seq, const size_t& readSize)
        {

            std::vector<uint16_t>& threadCountFwd = countsFwd[omp_get_thread_num()];
            std::vector<uint16_t>& threadCountRev = countsRev[omp_get_thread_num()];
            std::vector<uint16_t>& threadCountFwdStart = countsFwdStart[omp_get_thread_num()];
            std::vector<uint16_t>& threadCountRevStart = countsRevStart[omp_get_thread_num()];
            // fill with zeroes
            threadCountFwd.assign(ref.metaCpGs.size(), 0);
            threadCountRev.assign(ref.metaCpGs.size(), 0);
            threadCountFwdStart.assign(ref.metaStartCpGs.size(), 0);
            threadCountRevStart.assign(ref.metaStartCpGs.size(), 0);

            // retrieve kmers for first hash
            uint64_t fhVal = ntHash::NTP64(seq.data());

            uint64_t key = fhVal % MyConst::HTABSIZE;

            uint64_t lastId = 0xffffffffffffffffULL;
            bool wasFwd = false;
            bool wasStart = false;

            for (uint64_t i = ref.tabIndex[key]; i < ref.tabIndex[key+1]; ++i)
            {

                const uint64_t metaId = KMER::getMetaCpG(ref.kmerTable[i]);
                const bool isFwd = ref.strandTable[i];
                const bool isStart = KMER::isStartCpG(ref.kmerTable[i]);
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
                        ++threadCountFwd[metaId];

                    } else {

                        ++threadCountRev[metaId];

                    }
                }
            }

            for (unsigned int i = 0; i < (seq.size() - MyConst::KMERLEN); ++i)
            {

                // use rolling hash
                ntHash::NTP64(fhVal, seq[i], seq[i + MyConst::KMERLEN]);

                lastId = 0xffffffffffffffffULL;
                wasFwd = false;
                wasStart = false;

                for (uint64_t i = ref.tabIndex[key]; i < ref.tabIndex[key+1]; ++i)
                {

                    const uint64_t metaId = KMER::getMetaCpG(ref.kmerTable[i]);
                    const bool isFwd = ref.strandTable[i];
                    const bool isStart = KMER::isStartCpG(ref.kmerTable[i]);
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
                            ++threadCountFwd[metaId];

                        } else {

                            ++threadCountRev[metaId];

                        }
                    }
                }
            }


            // COUNTING HEURISTIC
            // if we retrieve too many seed positions passing the q-gram lemma filter, return flag to skip read


            // // set q-gram lemma threshold
            // uint16_t qThreshold = readSize - MyConst::KMERLEN - (MyConst::KMERLEN * MyConst::MISCOUNT);
            // // if too small we get an overflow -> set to zero
            // if (qThreshold > readSize)
            //     qThreshold = 0;
            //
            // unsigned int passCount = 0;
            // for (uint16_t& c : threadCountFwd)
            // {
            //
            //     if (c >= qThreshold)
            //     {
            //         ++passCount;
            //
            //     } else {
            //
            //         c = 0;
            //     }
            //
            // }
            // for (uint16_t& c : threadCountRev)
            // {
            //
            //     if (c >= qThreshold)
            //     {
            //         ++passCount;
            //
            //     } else {
            //
            //         c = 0;
            //     }
            //
            // }
            // // test if we have too many k-mers passing filter
            // if (passCount > MyConst::QUERYTHRESHOLD)
            // {
            //     return false;
            // } else {
            //     return true;
            // }
            return true;
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
        static constexpr unsigned int n = 400;
        std::ofstream statFile;
        std::ofstream countFile;
        void printStatistics(const std::vector<std::vector<KMER::kmer> > SeedsK);

        // input stream of file given as path to Ctor
        std::ifstream file;
        igzstream igz;

        // representation of the reference genome
        RefGenome& ref;

        // buffer holding MyConst::CHUNKSIZE many reads
        std::vector<Read> readBuffer;

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
        // TODO
        std::ofstream of;
};

#endif /* READQUEUE_H */
