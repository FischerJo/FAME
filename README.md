# FAME

**F**ast and **A**ccurate **ME**thylation Aligner for large mammalian genomes.
Carries out alignment and methylation calling for CpGs of Whole Genome Bisulfite Sequencing (WGBS) reads in one go without the need of intermediate alignment or buffer files. 
The algorithm is working on the full alphabet (A,C,G,T), resolving the asymmetric mapping problem\* correctly.
Reference genomes are expected to be .fasta files and reads in .fastq (or .fastq.gz) format.

The code is written in C++ and parallelized using OpenMP and is licensed under GPL3.

\*In WGBS experiments, unmethylated Cytosines are converted to Thymines in Reads, thus it is necessary to allow
Read Cytosines to map to Reference Thymines, but not vice versa. This is termed asymmetric mapping problem and is commonly tackled by working only on the reduced alphabet (A,T,G), which results in false matchings of reference Ts to read Cs.

## 0) News

This tool is still under development. We will update the news section once we published the algorithm and a stable version of the code is available.


## 1) Getting Started

The following steps explain how to install and setup FAME .


### A) Dependencies

FAME is dependent on the following libraries, which are all shipped with this repository.
So no additional work is necessary here.

* [ntHash](https://github.com/bcgsc/ntHash) - a fast library for genomic rolling hash functions.
* [sparsehash](https://github.com/sparsehash/sparsehash) - an efficient hash map implementation from Google.
* [gzstream](https://www.cs.unc.edu/Research/compgeom/gzstream/) - a C++ stream interface for working with gzip files

To compile the program, a recent version of the GNU Compiler Collection (gcc) is needed.
The gzstream library is dependent on [zlib](https://zlib.net/), which is usually installed on Linux systems.


### B) Installation


To retrieve the code base, just clone this repository:
```
git clone https://github.com/FischerJo/Metal
```

To install FAME on your machine, just use the shipped Makefile by typing
```
make
```
in the top level directory of the cloned repository.


### C) Simple example

FAME, as most other tools, uses a sohpisticated index structure of the reference genome to accelerate alignment computation.
The tool has two separate functions, the construction of the index and the alignment of the reads.
After construction of the index, it is stored in an efficient binary format and can thus be reused for different experiments.

To construct an index with the parameters specified in `CONST.h` (see 2.A):
```
./FAME --genome /Path/To/Genome/genome.fasta --store_index /Path/To/produced_index
```
where `genome.fasta` is a genome stored in fasta format and produced_index is the name of the output file holding the index after execution. For a full set of options see 2.A and 2.B.

To align reads to an index call:
```
./FAME -r /Path/To/reads.fastq --load_index /Path/To/produced_index -o /Path/To/output
```
with `/Path/To/output` being the output file path for the CpG report.
For a full set of options we refer to section 2.B.



## 2) Manual


### A) Parameter settings

All parameters are set in the file `CONST.h`.
If you change any of those parameters, make sure to rebuild the program:
```
make clean
make
```
An index is dependet on the parameters, that is, if you call the program
with an index that was built with different parameters, it will throw an error.

Here is a list of the external parameters:

| Parameter     | Definition       | Recommended value  | Location (line number) |
| ------------- |-------------| -----:| :----: |
| READLEN      | (Maximum) Length of reads queried to the index | 100 | 34 |
| CORENUM      | Number of threads spawned by the program. Should be number of free cores on the system. | 16 | 37 |
| MINPDIST | Minimum distance between a read pair in paired end mode. Measured from end to first read to beginning of second read.| 50 | 40 |
| MAXPDIST | Maximum distance between a read pair in paired end mode. Measured from end to first read to beginning of second read.| 400 | 41 |
| CHROMNUM | Number of chromosomes of reference organism. | 24 | 44 |

Here is a list of important internal parameters, we strongly recommend not to change them:

| Parameter     | Definition       | Recommended value  | Location (line number) |
| ------------- |-------------| -----:| :----: |
| CHUNKSIZE      | Number of reads (or read pairs) read to buffer. | 300000 | 64 |
| KMERLEN     | k, the length of a k-mer for the index. This is a very sensitive parameter. | 25 | 68 |
| MISCOUNT | Number of errors considered for k-mer filters. | 2 | 88 |
| ADDMIS | Number of errors additionally (to MISCOUNT) allowed in alignment | 3 | 90 |
| KMERCUTOFF | Controls hash collisions in index. Lowe value means more lossy but faster filter. Not considered if `--no_loss` flag is set. | 1500 | 93 |


### B) FAME command line arguments

The following is a description of all arguments accepted by FAME.
Examples on how to use FAME in the command line are given in 2.D.

| Flag    | Argument       | Description  |
| ------------- |-------------| :-----:|
| -h      | None | Lists all available options with a description. |
| --help | None | see -h |
| --genome | Filepath | Forces the tool to build an index for the specified .fasta reference file |
| --store_index | Filepath | Writes output of index construction to filepath (~32GB for human genome). NOTE: Directory must exist. |
| --no_loss | None | Builds the index with lossless filter (NOT RECOMMENDED). |
| --load_index | Filepath | Loads the index constructed before. NOTE: Parameters in `CONST.h` must be the same. |
| -r | Filepath | Forces the tool to query the specified single end read .fastq file to a loaded index. |
| -r1 | Filepath | Path to file with first reads of a paired read set. Read format must be .fastq. |
| -r2 | Filepath | Path to file with second reads of a paired read set. Read format must be .fastq. |
|--both_strands | None | Always try both, the read and its reverse complement for matching against the genome. This is required if reads are could potentially be C->T and G->A converted. |
| --gzip_reads | None | Treats the read files passed to -r or -r1 and -r2 as gzipped files. |
| -o | Filepath | Base name for output file. Contains CpG methylation levels after processing. |
| --out_basename | Filepath | see -o |


### C) Output format

The alignment of a read set to an index produces an output file (see `-o` command line argument) that contains
methylation levels of all CpGs in the reference genome.
The file is a tab-separated value file (`.tsv`) with the following structure:
```
Chromosome  Position   #Meth_Cs_fwd   #Unmeth_Cs_fwd   #Meth_Cs_rev   #Unmeth_Cs_rev
```
where #Meth_Cs_fwd means number of methylated Cytosines on forward strand CpG and so on.
The position is a (zero based) count of the bases of a chromosome, indicating the position of the C of a CpG.
The forward strand is the strand provided in the reference file, the reverse complement strand the strand not provided in the reference file.


### D) Extended examples

TODO



## 3) Contact

We appreciate any feedback to our tool.
Feel free to contact us by mailing to fischer 'at' mpi-inf.mpg.de.

If you found a bug, please contact us with a description of how you called the tool and which data you used.


## 4) Authors

* **Jonas Fischer** - *Main Contributor* - Max-Planck-Institute for Informatics, Saarbruecken, Germany
* **Marcel H. Schulz** Max-Planck-Institute for Informatics, Saarbruecken, Germany


## 5) License

This project is licensed under the GPL3 - see the [LICENSE_GPL_3_0](LICENSE_GPL_3_0) file for details.


## 6) References

* Our paper: COMING SOON
* ntHash: [ntHash: recursive nucleotide hashing. Mohamadi H, Chu J, Vandervalk BP, Birol I.  Bioinformatics 2016](https://www.ncbi.nlm.nih.gov/pubmed/27423894)

