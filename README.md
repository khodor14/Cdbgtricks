[INTRODUCTION](#ccdbgupdater)   
[INSTALLATION](#installation)   
[USAGE](#binary-usage)  
[EXAMPLES](#examples)   
[ALGORITHM](#algorithm)   
[CITATION](#citation)   
[CONTACT](#contact)

# Cdbgtricks
Cdbgtricks is a modular tool for indexing, updating and querying reads on a compacted de Bruijn graph. 
### Augmenting a compacted de Bruijn Graph by a new genome or a set of sequences
We support at the moment a fasta output format. 
## Requirements

ccdbgUpdater is developed in C++11
* C++11 compiler:
    * [GCC](https://gcc.gnu.org/) >= 4.8.5

* It relies on kmtricks to find the set of absent k-mers which are in the genome to be added but not in the graph

* [Boost](https://www.boost.org/)
* [ZSTD](https://anaconda.org/conda-forge/zstd)
## Installation
Installing kmtricks: First Install kmtricks [version v1.2.1]
* From [Bioconda](https://bioconda.github.io):

  ```
  conda install -c conda-forge -c tlemane kmtricks (latest version)
  ```
Installing ccdbupdater: Currently we have only the source version
* From source

  ```
  git clone --recursive https://github.com/khodor14/Cdbgtricks.git
  cd Cdbgtricks && mkdir build && cd build
  cmake -S ../ -B .
  make
  ```

## Binary usage:
* First change the number of files that can be oppened

  ```
  ulimit -n 2000
  ```

```
./cdbgtricks
```
displays the command line interface:
```
Usage: ./cdbgtricks [COMMAND] [PARAMETERS]

[COMMAND]:
	update 		 update a compacted de Bruijn graph by a new genome
	index 		 index a compacted de Bruijn graph by (k-1)-mers
	convert 		 convert a graph from GFA/FASTA/binary to binary/Fasta

	query 		 Query reads in Fasta/Fastq format (output presence/absence of reads)

	map		 Map reads in Fasta/Fastq format (output uni-MEMs)

[PARAMETERS]: update

	-h[--help]	 prints this help message
	--input_graph[-ig]	 the path to the pangenome graph in gfa or fasta format
	--input_genome	 the path to the genome used to augment the input graph
	--k_mer_size[-k] the size of the k-mer.
			 It must be the same value used when constructing the input graph
	--minimizer_size[-m] the size of the minimizer (m<k)
	--k_mer_file	 the file of absent k-mers from the graph if already computed
	--smallest_merge[-s]	 the threshold for merging buckets (all buckets with size smaller this threshold are merged)
	--log_super_bucket[-l]	 log2 of the number of files to be used
	--multiplier_super_bucket[-msb]	 size of super-bucket in terms of small bucket
	--output_file_name[-o] the name of the output file
	--update_index[-u] index the constructed funitigs
	--load_index[-li] the path to the saved index
	--output_index[-oi] write the index to a binary file
	--output_graph_binary[-ogb] write the graph in binary format
	--load_graph_binary[-lgb] load the graph from binary file
	-v verbosity
[PARAMETERS]: index

	-h[--help]	 prints this help message
	--input_graph[-ig]	 the path to the pangenome graph in gfa/fasta/binary format
	--k_mer_size[-k] the size of the k-mer
	--minimizer_size[-m] the size of the minimizer (m<k)
	--smallest_merge[-s]	 the threshold for merging buckets (all buckets with size smaller this threshold are merged)
	--log_super_bucket[-l]	 log2 of the number of files to be used
	--multiplier_super_bucket[-msb]	 size of super-bucket in terms of small bucket
	--log_super_bucket[-l]	 log2 of the number of files to be used
	 -it	 the input is in text format (fasta/GFA)
	 -ib	 the input is in binary format
	--output_file_name[-o] the name of the output file
	-v verbosity
[PARAMETERS]: convert

	-h[--help]	 prints this help message
	--input_graph[-ig]	 the path to the pangenome graph in gfa/fasta/binary format
	 -it	 the input is in text format (fasta/GFA)
	 -ib	 the input is in binary format
	 -of	 output in fasta format
	 -ob	 output in binary format
	--output_file_name[-o] the name of the output file
	-v verbosity
[PARAMETERS]: map

	-h[--help]	 prints this help message
	--input_graph[-ig]	 the path to the pangenome graph in gfa/fasta/binary format
	--query_reads[-qr]	 the path to the to the read query set in fasta format
	--k_mer_size[-k] the size of the k-mer
	--minimizer_size[-m] the size of the minimizer (m<k)
	--smallest_merge[-s]	 the threshold for merging buckets (all buckets with size smaller this threshold are merged)
	--log_super_bucket[-l]	 log2 of the number of files to be used
	--multiplier_super_bucket[-msb]	 size of super-bucket in terms of small bucket
	--log_super_bucket[-l]	 log2 of the number of files to be used
	 -it	 the input is in text format (fasta/GFA)
	 -ib	 the input is in binary format
	--load_index[-li] the path to the saved index
	--ratio[-r]	 the ratio of read k-mers that should be found in the graph
	--output_file_name[-o]the name of the output file
	-v verbosity
[PARAMETERS]: query

	-h[--help]	 prints this help message
	--input_graph[-ig]	 the path to the pangenome graph in gfa/fasta/binary format
	--query_reads[-qr]	 the path to the to the read query set in fasta format
	--k_mer_size[-k] the size of the k-mer
	--minimizer_size[-m] the size of the minimizer (m<k)
	--smallest_merge[-s]	 the threshold for merging buckets (all buckets with size smaller this threshold are merged)
	--log_super_bucket[-l]	 log2 of the number of files to be used
	--multiplier_super_bucket[-msb]	 size of super-bucket in terms of small bucket
	--log_super_bucket[-l]	 log2 of the number of files to be used
	 -it	 the input is in text format (fasta/GFA)
	 -ib	 the input is in binary format
	--load_index[-li] the path to the saved index
	--ratio[-r]	 the ratio of read k-mers that should be found in the graph
	--output_file_name[-o]the name of the output file
	-v verbosity
```
### Examples
  1. **Indexing a compacted de Bruijn Graph**
     ```
     ./cdbgtricks index -ig ../test/graph.gfa -k 31 -s 5000 -msb 1 -m 11 -o graph_index -it
     ```
     The 31-mers (`-k 31`) of the compacted de Bruijn graph *graph.gfa* is indexed  using a minimizer size 11 (`-m 11`). We set here rho to be 5000 (`-s 5000`).
     It means that we should have at least 5000 31-mers in a bucket sharing the same 11-mer minimizer to create an MPHF f. We set gamma to be 1 (`-msb 1`). It means that the size of the super-bucket is 5000. (`-it`) means that the graph is in text format. The output is written to graph_index_index.bin

  2. **Update a compacted de Bruijn graph from a reference genome file with loading the graph index**
     ```
     ./cdbgtricks update -ig ../test/graph.gfa --input_genome ../test/genome.fa -k 31 -li graph_index_index.bin -o graph_update -u -oi -t 32
     ```
     The compacted de Bruijn graph is updated with the absent 31-mers (`-k 31`) of file *genome.fa* (`--input_genome genome.fa`). Here [kmtricks](https://github.com/tlemane/kmtricks) is used to find the absent k-mers. The index is loaded (`-li`) and used to update the graph. The graph is written to file *graph_update.fa* (`-o graph_update`). The index is updated (`-u`) and written to graph_update_index.bin.
  3. **Update a compacted de Bruijn graph without loading its index**
     ```
     ./cdbgtricks update -ig ../test/graph.gfa --input_genome ../test/genome.fa -k 31 -o graph_update -u -oi -t 32
     ```
     The compacted de Bruijn graph is updated with the absent 31-mers (`-k 31`) of file *genome.fa* (`--input_genome genome.fa`). Here [kmtricks](https://github.com/tlemane/kmtricks) is used to find the absent k-mers. The graph is written to file *graph_update.fa* (`-o graph_update`). The index is updated (`-u`) and written to graph_update_index.bin.
  4. **Query a set of reads (output presence/absence)**
    ```
    ./cdbgtricks query -ig ../test/graph.gfa -qr queries.fa -k 31 -r 1 -li graph_index_index.bin -o results_query
    ```
    This tries to find if the reads are present or absent. A read is present if it shares a ratio of k-mers equal to r `-r 1`.
    The output is written to results_query.tsv
  5. **Map a set of reads (output uni-MEMs)**
  ```
  ./cdbgtricks map -ig ../test/graph.gfa -qr unitigs_fasta.fa -k 31 -r 1 -li graph_index_index.bin -o results_mapping
  ```
  Instead of outputing if a read is present or absent, it outputs the uni-MEMs of the read. A uni-MEM is a 5-tuple (unitig id, unitig start, unitig end,read start, read end). The output is written to results_mapping.txt
### Algorithm
  1. **Using kmtricks**
    Three commands of kmtricks are used to find the absent kmers
     ```
     kmtricks pipeline --file fof.txt --run-dir graph_dir  --kmer-size 31 --hard-min 1 --mode kmer:pa:bin --nb-partitions 100 --cpr --threads 32
     ```
     ```
     kmtricks filter --in-matrix graph_dir --key query.txt --output filter --hard-min 1 --out-types k --cpr-in --cpr-out -t 32
     ```
     ```
     kmtricks aggregate --run-dir filter --count A1:kmer --format text --cpr-in --output kmers.txt
     ```
*For more details on the documentation of kmtricks visit [wiki](https://github.com/tlemane/kmtricks/wiki)*

  2. **Indexing the input Graph**
    The k-mers of the graph are distributed into buckets/super-buckets. Each bucket/super-bucket is associated to a minimal perfect hash function (MPHF) f.
    The MPHF f provides an identifier i for a k-mer x (i=f(x)) which points to the position of the k-mer in the graph (unitig id,unitig offset,orientation) 
  3. **Constructing the future unitigs (funitigs) from absent k-mers with respect to the graph**
     We go over the absent k-mers, we try to extend each by its in neighbor k-mer from left and right. This process is implemented by following
     the node centric definition of a de bruijn graph and by ensuring that the extended unitig does not branch to graph. 
     If it branches or it breaks the definition, we break the extension. We use a flag to tell if a k-mer was used before or not
    
  4. **Updating the input graph**
    The construction of the funitigs provides the target unitigs that need to be split or merged with some funitigs. We go over theses unitigs and split or join them with the target funitigs.
  5. **Updating the index of the graph**
    If new k-mers are added to some buckets/super-buckets, their MPHFs are re-computed.


## Contact

For any question, feedback or problem, please feel free to write an issue on this GitHub repository and we will get back to you as soon as possible.
