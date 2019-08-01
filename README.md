# FAME

**F**ast and **A**ccurate **ME**thylation Aligner for large mammalian genomes.
Carries out alignment and methylation calling for CpGs of Whole Genome Bisulfite Sequencing (WGBS) reads in one go without the need of intermediate alignment or buffer files. 
The algorithm is working on the full alphabet (A,C,G,T), resolving the asymmetric mapping problem\* correctly.
Reference genomes are expected to be .fasta files and reads in .fastq (or .fastq.gz) format.

The code is written in C++ and parallelized using OpenMP and is licensed under GPL3.

\*In WGBS experiments, unmethylated Cytosines are converted to Thymines in Reads, thus it is necessary to allow
Read Cytosines to map to Reference Thymines, but not vice versa. This is termed asymmetric mapping problem and is commonly tackled by working only on the reduced alphabet (A,T,G), which results in false matchings of reference Ts to read Cs.

## 0) News

This tool is still under development.
There is now an alpha release for single cell support.

## 1) Getting Started

The following steps explain how to install and setup FAME .


### A) Dependencies

FAME is dependent on the following libraries, which are all shipped with this repository.
So no additional work is necessary here.

* [ntHash](https://github.com/bcgsc/ntHash) - a fast library for genomic rolling hash functions.
* [sparsehash](https://github.com/sparsehash/sparsehash) - an efficient hash map implementation from Google.
* [hopscotch-map](https://github.com/Tessil/hopscotch-map) - an efficient hash map implementation by Tessil.
* [gzstream](https://www.cs.unc.edu/Research/compgeom/gzstream/) - a C++ stream interface for working with gzip files.

To compile the program, a recent version of the GNU Compiler Collection (gcc) or LLVM Clang (clang) is required.
The gzstream library is dependent on [zlib](https://zlib.net/), which is usually installed on Linux systems.


### B) Installation


To retrieve the code base, just clone this repository:
```
git clone https://github.com/FischerJo/FAME
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
An index is dependent on the parameters, that is, if you call the program
with an index that was built with different parameters, it will throw an error.

Here is a list of the external parameters:

| Parameter     | Definition       | Recommended value  | Location (line number) |
| ------------- |-------------| -----:| :----: |
| READLEN      | (Maximum) Length of reads queried to the index | 100 | 35 |
| CORENUM      | Number of threads spawned by the program. Should be number of free cores on the system. | 16 | 38 |
| MINPDIST | Minimum distance between a read pair in paired end mode. Measured from end to first read to beginning of second read.| 20 | 41 |
| MAXPDIST | Maximum distance between a read pair in paired end mode. Measured from end to first read to beginning of second read.| 400 | 42 |
| CHROMNUM | Number of chromosomes of reference organism. | 24 | 45 |

Here is a list of some important internal parameters, we strongly recommend *NOT* to change them:

| Parameter     | Definition       | Recommended value  | Location (line number) |
| ------------- |-------------| -----:| :----: |
| SEED | The gapped q-gram (aka seed) to use as array of bits | \[see below\] | 57 |
| SEEDBITS | The gapped q-gram as bitstring | 0b11011110111111011111111111111101 | 58 |
| CHUNKSIZE      | Number of reads (or read pairs) read to buffer. | 300000 | 69 |
| KMERLEN     | k, the length of a k-mer for the index. This is a very sensitive parameter. Must be the length of the seed. | 32 | 73 |
| QTHRESH | Minimum number of k-mer matches required for match verification | 5 | 80 |
| WINLEN | Window length for the index data structure. | 2048 | 90 |
| MISCOUNT | Number of errors considered for k-mer filters. | 2 | 94 |
| ADDMIS | Number of errors additionally (to MISCOUNT) allowed in alignment | 4 | 97 |
| KMERCUTOFF | Controls hash collisions in index. Low value means more lossy but faster filter. Not considered if `--no_loss` flag is set during index construction. | 1500 | 100 |
| KMERDIST | Controls pruning after matching a read to a Window. Minimum distance of count of window to prune and count of matched window. | 10 | 104 |
| SKIPMOD | Hash only every SKIPMODth k-mer of reference. | 2 | 107 |


### B) FAME command line arguments

The following is a description of all arguments accepted by FAME.
Examples on how to use FAME in the command line are given in 2.D.

| Flag    | Argument       | Description  |
| ------------- |-------------| :-----:|
| --genome | Filepath | Forces the tool to build an index for the specified .fasta reference file |
| --gzip_reads | None | Treats the read files passed to -r or -r1 and -r2 as gzipped files. |
| -h      | None | Lists all available options with a description. |
| --help | None | see -h |
| --human_opt | None | Uses optimizations for human reference genomes to prune away unlocalized contigs etc |
| --load_index | Filepath | Loads the index constructed before. NOTE: Parameters in `CONST.h` must be the same. |
| --no_loss | None | Builds the index with lossless filter (NOT RECOMMENDED). |
| -o | Filepath | Base name for output file. Contains CpG methylation levels after processing, bulked values for single cell mode. |
| --out_basename | Filepath | see -o |
| --paired | None | Flag for single cell mode that indicates that cells are paired-end sequenced. |
| -r | Filepath | Forces the tool to query the specified single end read .fastq file to a loaded index. |
| -r1 | Filepath | Path to file with first reads of a paired read set. Read format must be .fastq. |
| -r2 | Filepath | Path to file with second reads of a paired read set. Read format must be .fastq. |
| --sc_input | Filepath | Path to file containing meta information on single cell data (see Section 2E). |
| --sc_output | Filepath | Name for output file of single cell mode. |
| --store_index | Filepath | Writes output of index construction to filepath (~38GB for human genome). NOTE: Directory must exist. |
| --unord_reads | None | Disable optimization to find stranding of reads. |



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

Suppose you want to match paired reads to a human genome and save it to the file with prefix pairedResults. The reads are in gzip format.

First build the index:

```
./FAME --genome /Path/To/Genome/genome.fasta --store_index /Path/To/produced_index --human_opt
```

Align the reads:
```
./FAME --load_index /Path/To/produced_index -r1 /Path/To/r1.fastq.gz -r2 /Path/To/r2.fastq.gz --gzip_reads -o pairedResults
```
A file called pairedResults_cpg.tsv is generated containing the CpGs with corresponding methylation counts.


Align single cell files:
```
./FAME --load_index /Path/To/produced_index --sc_input /Path/To/Metafile --sc_output stratified_results -o bulk_results
```
Produces counts for each single cell specified in the Metafile and outputs them as separate rows into 'stratified_results'.
A summary of summed up methylation counts for each CpG across all cells (i.e. bulk values) are printed to bulk_results.

### E) Single Cell Meta File

To process single cell data, FAME requires a simple tsv file with meta information with 2 (3) columns for single-end (paired-end) single cell experiments.
The first column contains a unique identifier for the cells, which is also used for the output.
The second column contains the filepath to the .fastq read file for the cell.
In case of paired-end experiments, the second read file is expected in the third column.
Such a metafile could thus look like this:
```
Cell_ID1	/MyPath/Cell_ID1.read1.fastq	/MyPath/Cell_ID1.read2.fastq
Cell_ID2	/MyPath/Cell_ID2.read1.fastq	/MyPath/Cell_ID2.read2.fastq
```


The single cell output is a tsv file consisting of 4 rows per cell. The columns are CpGs, which locations are indicated by a two line header.
On the first line of the header, the chromosome of each CpG is given. The second line specifies the position within the chromosome, starting at 0.
The first two columns contain the cell id and the count type, respectively. There are 4 count types for each CpG, similar to bulk processing, the count of
mapped methylated Cs to the forward strand, unmethylated Cs to forward strand, methylated Cs to the reverse strand and unmethylated Cs to the reverse strand.
A dummy output file could look like this:

```
SC_ID	Count_Type	chr1	chr1	chr1	chr2
SC_ID	Count_Type	4405	4410	4443	1123
Cell_ID1	methFwd	1	4	5	5	3
Cell_ID1	unmethFwd	0	20	22	4	0
Cell_ID1	methRev	1	4	6	2	2
Cell_ID1	unmethRev	1	18	15	3	1
Cell_ID2	methFwd	0	4	7	5	3
Cell_ID2	unmethFwd	0	2	1	4	1
Cell_ID2	methRev	1	4	22	13	2
Cell_ID2	unmethRev	1	2	3	3	0
```

## 3) Contact

We appreciate any feedback to our tool.
Feel free to contact us by mailing to fischer 'at' mpi-inf.mpg.de.

If you found a bug, please create an issue in this github repo with a description of how you called the tool and which data you used.


## 4) Authors

* **Jonas Fischer** - *Main Contributor* - Max-Planck-Institute for Informatics, Saarbruecken, Germany
* **Marcel H. Schulz** Max-Planck-Institute for Informatics, Saarbruecken, Germany


## 5) License

This project is licensed under the GPL3 - see the [LICENSE_GPL_3_0](LICENSE_GPL_3_0) file for details.


## 6) References

* Our paper: COMING SOON
* ntHash: [ntHash: recursive nucleotide hashing. Mohamadi H, Chu J, Vandervalk BP, Birol I.  Bioinformatics 2016](https://www.ncbi.nlm.nih.gov/pubmed/27423894)

