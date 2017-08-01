
#include <chrono>


#include "RefReader_istr.h"
#include "RefGenome.h"
#include "ReadQueue.h"

void queryRoutine(ReadQueue& rQue);
void printHelp();

// --------------- MAIN -----------------
//
//
int main(int argc, char** argv)
{

    std::string indexFile = "";
    std::string genomeFile = "";
    char* readFile = NULL;

    // true iff index should be loaded from file
    bool loadIndexFlag = false;
    // true iff index should be stored
    bool storeIndexFlag = false;
    // ture iff reads are in .gz format
    bool readsGZ = false;

    if (argc == 1)
    {
        printHelp();
        return 0;
    }

    for (int i = 1; i < argc; ++i)
    {

        if (std::string(argv[i]) == "-h" || std::string(argv[i]) == "--help")
        {
            printHelp();
            return 0;
        }

        if (std::string(argv[i]) == "-r")
        {
            if (i + 1 < argc)
            {
                readFile = argv[++i];

            } else {

                std::cerr << "No filepath for option \"" << argv[i] << "\" provided! Terminating...\n\n";
                exit(1);
            }
            continue;
        }

        if (std::string(argv[i]) == "--genome")
        {
            if (i + 1 < argc)
            {
                genomeFile = std::string(argv[++i]);

            } else {

                std::cerr << "No filepath for option \"" << argv[i] << "\" provided! Terminating...\n\n";
                exit(1);
            }
            continue;
        }

        if (std::string(argv[i]) == "--load_index")
        {
            if (i + 1 < argc)
            {
                if (storeIndexFlag)
                {
                    std::cerr << "Cannot store AND load index! Terminating...\n\n";
                }
                indexFile = std::string(argv[++i]);
                loadIndexFlag = true;

            } else {

                std::cerr << "No filepath for option \"" << argv[i] << "\" provided! Terminating...\n\n";
                exit(1);
            }
            continue;
        }

        if (std::string(argv[i]) == "--store_index")
        {
            if (i + 1 < argc)
            {
                if (loadIndexFlag)
                {
                    std::cerr << "Cannot store AND load index! Terminating...\n\n";
                }
                indexFile = std::string(argv[++i]);
                storeIndexFlag = true;

            } else {

                std::cerr << "No filepath for option \"" << argv[i] << "\" provided! Terminating...\n\n";
                exit(1);
            }
            continue;
        }

        if (std::string(argv[i]) == "--gzip_reads")
        {
            readsGZ = true;
            continue;
        }

    }

    if (loadIndexFlag)
    {
        RefGenome ref(indexFile);

        if (readFile == NULL)
        {

            std::cerr << "No read file provided! Use \"-r\" option to specify read file path. Terminating...\n\n";
            exit(1);

        }
        ReadQueue rQue(readFile, ref, readsGZ);

        queryRoutine(rQue);

    } else {

        if (genomeFile == "")
        {
            std::cerr << "No reference genome file provided! Use \"--genome\" option to specify reference genome file path. Terminating...\n\n";
            exit(1);
        }

        std::vector<struct CpG> cpgTab;
        std::vector<struct CpG> cpgStartTab;
        std::vector<std::vector<char>> genSeq;

        readReference(genomeFile, cpgTab, cpgStartTab, genSeq);
        RefGenome ref(std::move(cpgTab), std::move(cpgStartTab), genSeq);

        if (storeIndexFlag)
        {
            ref.save(indexFile);

        } else if (readFile == NULL)
        {

            std::cerr << "No read file provided and no index storing option! Use \"-r\" option to specify read file path. Terminating...\n\n";
            exit(1);

        }

        if (readFile != NULL)
        {
            ReadQueue(readFile, ref, readsGZ);
        }
    }


    // std::vector<struct CpG> cpgTab;
    // std::vector<struct CpG> cpgStartTab;
    // std::vector<std::vector<char>> genSeq;
    //
    // readReference(argv[1], cpgTab, cpgStartTab, genSeq);
    // RefGenome ref(std::move(cpgTab), std::move(cpgStartTab), genSeq);
    // // ref.save(binFile);
    // // RefGenome refCopy(binFile);
    // // ReadQueue rQue(argv[2], refCopy, false);
    // // ReadQueue rQue(argv[2], refCopy, true);
    // ReadQueue rQue(argv[2], ref, false);
    //
    // unsigned int readCounter = 0;
    // unsigned int i = 0;
    //
    // rQue.parseChunk(readCounter);
    // rQue.matchReads(readCounter);
    // // rQue.parseChunkGZ(readCounter);
    // // rQue.matchReads(readCounter);
    // // std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
    // // while(rQue.parseChunk(readCounter))
    // // {
    // //     ++i;
    // //     if (i >= 10)
    // //         break;
    // //     rQue.matchReads(readCounter);
    // //     std::cout << "Processed " << MyConst::CHUNKSIZE * i << " many reads\n";
    // // }
    // // // match remaining reads
    // // rQue.matchReads(readCounter);
    // //
    // // std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
    // // auto runtime = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();
    // //
    // // std::cout << "Done processing in " << runtime << "s\n";
    //
    return 0;
}

//
//
//---------------------------------------




void queryRoutine(ReadQueue& rQue)
{

    unsigned int readCounter = 0;
    unsigned int i = 0;
    std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
    while(rQue.parseChunk(readCounter))
    {
        ++i;
        // TODO
        if (i >= 11)
            break;
        rQue.matchReads(readCounter);
        std::cout << "Processed " << MyConst::CHUNKSIZE * i << " many reads\n";
    }
    // match remaining reads
    rQue.matchReads(readCounter);

    std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
    auto runtime = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();

    std::cout << "Done processing in " << runtime << "s\n";

}

void printHelp()
{

    std::cout << "\n\nSUMMARY\n\n";

    std::cout << "\tThis program is designed for the computation of methylation levels\n";
    std::cout << "\tin (large) mammalian genomes. Please specify the desired values for\n";
    std::cout << "\tthe parameters in the file CONST.h, and rebuild the project using\n";
    std::cout << "\tthe provided Makefile (first \"make clean\" then \"make\").\n";
    std::cout << "\tAll files are provided as command line arguments by the user.\n\n\n";


    std::cout << "OPTIONS\n\nOptions followed by [.] require an additional argument\n\n";

    std::cout << "\t--help\n";
    std::cout << "\t-h               \t\tHelp\n\n";

    std::cout << "\t--genome      [.]\t\tSpecification of a filepath to a reference genome\n";
    std::cout << "\t                 \t\tin fasta format.\n\n";

    std::cout << "\t-r            [.]\t\tSpecification of a filepath to a set of reads in\n";
    std::cout << "\t                 \t\tfastq format. If not specified, index is built and\n";
    std::cout << "\t                 \t\tsaved in file provided via --store_index.\n\n";

    std::cout << "\t--gzip_reads     \t\tRead file specified by -r is treated as gzipped\n";
    std::cout << "\t                 \t\tfile (.gz file ending).\n\n";

    std::cout << "\t--store_index [.]\t\tStore index in provided file in binary format.\n\n";

    std::cout << "\t--load_index  [.]\t\tLoad index from provided file. Note that all\n";
    std::cout << "\t                 \t\tparameters used to build the index must be the same\n";
    std::cout << "\t                 \t\tas used in the current CONST.h. This will be checked\n";
    std::cout << "\t                 \t\twhile loading.\n\n";


    std::cout << "\nEXAMPLES\n\n";

    std::cout << "Setting: Read a reference genome and map a set of reads, save index for\n";
    std::cout << "consecutive usages.\n\n";
    std::cout << "\t /path/to/Metal --genome /path/to/reference.fasta --store_index index.bin \\\n";
    std::cout << "\t\t-r /path/to/reads.fastq\n\n\n";

    std::cout << "Setting: Load index from previously stored index, maap reads stored in .gz\n";
    std::cout << "format.\n\n";
    std::cout << "\t /path/to/Metal --load_index index.bin -r /path/to/reads.fastq.gz\n\n\n";

    std::cout << "\n\n";
}
