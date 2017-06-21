#ifndef READQUEUE_H
#define READQUEUE_H

#include <string>
#include <fstream>
#include <unordered_map>
#include <array>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "CONST.h"
#include "RefGenome.h"
#include "Read.h"
#include "ShiftAnd.h"

class ReadQueue
{

    public:

        ReadQueue(const char* filePath, RefGenome& ref);


        // Parses a chunk of the ifstream file
        // Reads up to MyConst::CHUNKSIZE many reads and saves them
        // ARGUMENT:
        //          procReads   will contain number of reads that have been read into buffer
        // returns true if neither read error nor EOF occured, false otherwise
        bool parseChunk(unsigned int& procReads);

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
                    if (lastStrand)
                    {
                        ++threadCountFwd[metaId];

                    } else {

                        ++threadCountRev[metaId];

                    }
                }
            }

            // More than cutoff many kmers are required per metaCpG
            const unsigned int countCut = readSize - MyConst::KMERLEN + 1 - (MyConst::KMERLEN * MyConst::MISCOUNT);


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
        //              true iff unique best match found
        //              false otherwise
        //
        // MODIFICATIONS:
        //              updates the matches field of the read set
        //
        inline bool saQuerySeedSet(ShiftAnd<MyConst::MISCOUNT>& sa, std::vector<std::vector<KMER::kmer> >& seedsK, std::vector<std::vector<bool> >& seedsS, MATCH::match mat)
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

            // flags stating if we had a match with x errors before
            // if we had a match with i errors then the MyConst::MISCOUNT - i 'th lower most bit and all lower bits are set
            // that is, we have the number of allowed mismatches + 1 lower bits blocked and the most significant of those
            // indicates 0 errors. All lower bits indicate 1, 2... and so on many errors.
            // uint64_t multiMatch = 0;

            // counter for how often we had a match
            std::vector<uint8_t> multiMatch(MyConst::MISCOUNT + 1, 0);

            // will contain matches iff match is found for number of errors specified by index
            std::vector<MATCH::match> uniqueMatches(MyConst::MISCOUNT + 1);

            for (size_t outerNdx = 0; outerNdx < seedsK.size(); ++outerNdx)
            {

                for (size_t innerNdx = 0; innerNdx < seedsK[0].size(); ++innerNdx)
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
                            if (threadCountFwd[m] == 1 || threadCountFwd[m] == 3)
                            {
                                continue;

                            } else {

                                // state that we queried this
                                threadCountFwd[m] += 1;
                                // retrieve it
                                const struct CpG& startCpg = ref.cpgStartTable[ref.metaStartCpGs[m].start];
                                const struct CpG& endCpg = ref.cpgStartTable[ref.metaStartCpGs[m].end];

                                auto startIt = ref.fullSeq[startCpg.chrom].begin();
                                auto endIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + MyConst::READLEN;

                                // check if CpG was too near to the end
                                if (endIt > ref.fullSeq[startCpg.chrom].end())
                                {
                                    // if so move end iterator appropriately
                                    endIt = ref.fullSeq[startCpg.chrom].end();
                                }

                                // use shift and to find all matchings
                                std::vector<uint64_t> matchings;
                                std::vector<uint16_t> errors;
                                sa.querySeq(startIt, endIt, matchings, errors);

                                // go through matching and see if we had such a match (with that many errors) before - if so,
                                // return to caller reporting no match
                                for (size_t i = 0; i < matchings.size(); ++i)
                                {

                                    // check if we had a match with that many errors before
                                    if (multiMatch[errors[i]])
                                    {

                                        // check if this is a match without errors
                                        if (!errors[i])
                                        {

                                            // if so, return without a match
                                            return false;

                                        }
                                        // set the number of matches with that many errors to 2
                                        // indicating that we do not have a unique match with that many errors
                                        multiMatch[errors[i]] = 2;


                                    } else {

                                        // we don't have such a match yet,
                                        // so save this match at the correct position
                                        uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i], startCpg.chrom, errors[i], 1);
                                        ++multiMatch[errors[i]];
                                    }


                                }

                            }

                        // is not forward strand meta CpG
                        } else {

                            // check if we queried it before
                            if (threadCountRev[m] == 1 || threadCountRev[m] == 3)
                            {
                                continue;

                            } else {

                                threadCountRev[m] += 1;
                                // retrieve it
                                const struct CpG& startCpg = ref.cpgStartTable[ref.metaStartCpGs[m].start];
                                const struct CpG& endCpg = ref.cpgStartTable[ref.metaStartCpGs[m].end];

                                auto endIt = ref.fullSeq[startCpg.chrom].begin() - 1;
                                auto startIt = ref.fullSeq[startCpg.chrom].begin() + endCpg.pos + MyConst::READLEN - 1;

                                // check if CpG was too near to the end
                                if (startIt >= ref.fullSeq[startCpg.chrom].end())
                                {
                                    // if so move end iterator appropriately
                                    startIt = ref.fullSeq[startCpg.chrom].end() - 1;
                                }

                                // use shift and to find all matchings
                                std::vector<uint64_t> matchings;
                                std::vector<uint16_t> errors;
                                sa.queryRevSeq(startIt, endIt, matchings, errors);

                                // go through matching and see if we had such a match (with that many errors) before - if so,
                                // return to caller reporting no match
                                for (size_t i = 0; i < matchings.size(); ++i)
                                {

                                    // check if we had a match with that many errors before
                                    if (multiMatch[errors[i]])
                                    {

                                        // check if this is a match without errors
                                        if (!errors[i])
                                        {

                                            // if so, return without a match
                                            return false;

                                        }
                                        // set the number of matches with that many errors to 2
                                        // indicating that we do not have a unique match with that many errors
                                        multiMatch[errors[i]] = 2;


                                    } else {

                                        // we don't have such a match yet,
                                        // so save this match at the correct position
                                        uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i], startCpg.chrom, errors[i], 1);
                                        ++multiMatch[errors[i]];
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
                                const struct CpG& startCpg = ref.cpgTable[ref.metaStartCpGs[m].start];
                                const struct CpG& endCpg = ref.cpgTable[ref.metaStartCpGs[m].end];

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
                                std::vector<uint16_t> errors;
                                sa.querySeq(startIt, endIt, matchings, errors);

                                // go through matching and see if we had such a match (with that many errors) before - if so,
                                // return to caller reporting no match
                                for (size_t i = 0; i < matchings.size(); ++i)
                                {

                                    // check if we had a match with that many errors before
                                    if (multiMatch[errors[i]])
                                    {

                                        // check if this is a match without errors
                                        if (!errors[i])
                                        {

                                            // if so, return without a match
                                            return false;

                                        }
                                        // set the number of matches with that many errors to 2
                                        // indicating that we do not have a unique match with that many errors
                                        multiMatch[errors[i]] = 2;


                                    } else {

                                        // we don't have such a match yet,
                                        // so save this match at the correct position
                                        uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i] + startCpg.pos, startCpg.chrom, errors[i], 1);
                                        ++multiMatch[errors[i]];
                                    }
                                }
                            }

                        // kmer is on backward strand
                        } else {

                            // check if we queried this meta CpG it before
                            if (threadCountFwd[m] >= 2)
                            {
                                continue;

                            } else {

                                // state that we queried this
                                threadCountFwd[m] += 2;
                                // retrieve it
                                const struct CpG& startCpg = ref.cpgTable[ref.metaStartCpGs[m].start];
                                const struct CpG& endCpg = ref.cpgTable[ref.metaStartCpGs[m].end];

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
                                std::vector<uint16_t> errors;
                                sa.queryRevSeq(startIt, endIt, matchings, errors);

                                // go through matching and see if we had such a match (with that many errors) before - if so,
                                // return to caller reporting no match
                                for (size_t i = 0; i < matchings.size(); ++i)
                                {

                                    // check if we had a match with that many errors before
                                    if (multiMatch[errors[i]])
                                    {

                                        // check if this is a match without errors
                                        if (!errors[i])
                                        {

                                            // if so, return without a match
                                            return false;

                                        }
                                        // set the number of matches with that many errors to 2
                                        // indicating that we do not have a unique match with that many errors
                                        multiMatch[errors[i]] = 2;


                                    } else {

                                        // we don't have such a match yet,
                                        // so save this match at the correct position
                                        uniqueMatches[errors[i]] = MATCH::constructMatch(matchings[i] + startCpg.pos, startCpg.chrom, errors[i], 1);
                                        ++multiMatch[errors[i]];
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
                // there is no much with that few errors, search the one with more errors
                if (multiMatch[i] == 0)
                {
                    continue;
                }
                // if match is not unique, return unsuccessfull to caller
                if (multiMatch[i] > 1)
                {
                    return false;

                } else {

                    mat = uniqueMatches[i];
                    return true;
                }

            }
            // we have not a single match at all, return unsuccessfull to caller
            return false;
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
        static constexpr unsigned int n = 400;
        std::ofstream statFile;
        std::ofstream countFile;
        void printStatistics(const std::vector<std::vector<KMER::kmer> > SeedsK);

        // input stream of file given as path to Ctor
        std::ifstream file;

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
};

#endif /* READQUEUE_H */
