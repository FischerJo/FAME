# FAME

Fast and Accurate Methylation Aligner for large mammalian genomes.
Carries out alignment and methylation calling of Whole Genome Bisulfite Sequencing (WGBS) reads in one go without the need of intermediate alignment or buffer files. 
The algorithm is working on the full alphabet (A,C,G,T), resolving the asymmetric mapping problem\* correctly.
Reference genomes are expected to be .fasta files and reads in .fastq (or .fastq.gz) format.

The code is written in C++ and parallelized using OpenMP and is licensed under GPL3.

\*In WGBS experiments, unmethylated Cytosines are converted to Thymines in Reads, thus it is necessary to allow
Read Cytosines to Reference Thymines, but not vice versa. This is termed asymmetric mapping problem and is commonly tackled by working only on the reduced alphabet (A,T,G), which results in false matchings of reference Ts to read Cs.

## Getting Started

The following steps explain how to install, setup, and run FAME.

### Dependencies

FAME is dependent on the following libraries, which are all shipped with this repository.
So no additional work is necessary here.

*[ntHash](https://github.com/bcgsc/ntHash) - a fast library for genomic rolling hash functions.
*[sparsehash](https://github.com/sparsehash/sparsehash) - an efficient hash map implementation from Google.
*[gzstream](https://www.cs.unc.edu/Research/compgeom/gzstream/) - a C++ stream interface for working with gzip files

To compile the program, a recent version of the GNU Compiler Collection (gcc) is needed.
The gzstream library is dependent on [zlib](https://zlib.net/), which is usually installed on a Linux system.


### Installing


To retrieve the code base, just clone this repository:
```
git clone https://github.com/FischerJo/Metal
```

To install FAME on your machine, just use the shipped Makefile by typing
```
make
```
in the top level directory of the cloned repository.

### Simple example

TODO
## Manual
TODO



## Authors

* **Jonas Fischer** Max-Planck-Institute for Informatics, Saarbruecken, Germany
* **Marcel H. Schulz** Max-Planck-Institute for Informatics, Saarbruecken, Germany


## License

This project is licensed under the GPL3 - see the [LICENSE_GPL_3_0](LICENSE_GPL_3_0) file for details.

## References

* Our paper: COMING SOON
* ntHash: [ntHash: recursive nucleotide hashing. Mohamadi H, Chu J, Vandervalk BP, Birol I.  Bioinformatics 2016](https://www.ncbi.nlm.nih.gov/pubmed/27423894)

