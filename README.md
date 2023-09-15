[INTRODUCTION](#ccdbgupdater)   
[INSTALLATION](#installation)   
[USAGE](#binary-usage)  
[EXAMPLES](#examples)   
[ALGORITHM](#algorithm)   
[CITATION](#citation)   
[CONTACT](#contact)

# ccdbgUpdater
ccdbgUpdater is a modular tool for updating a colored and compacted de Bruijn Graph. 
This tool is on progress. 
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
  git clone --recursive https://github.com/khodor14/ccdbgUpdater.git
  cd ccdbgUpdater && mkdir build && cd build
  cmake -S ../ -B .
  make
  ```

## Binary usage:

```
./ccdbgupdater
```

displays the command line interface:
```
Usage: ./ccdbgupdater [COMMAND] [PARAMETERS]

./ccdbgupdater --input_graph <value> --input_genome <value> --k_mer_size <value>
OR
./ccdbgupdater --input_graph <value> --k_mer_file <value> --k_mer_size <value>



[COMMAND]:
	build 		 build a compacted and colored de Bruijn Graph
	update 		 update a compacted and colored de Bruijn graph by a new genome
	index 		 index a compacted de Bruijn graph by (k-1)-mers
	convert 	 convert a graph from GFA/FASTA/binary to binary/Fasta

[PARAMETERS]: build

	-h[--help]	 prints this help message
	--input_genome [-igr]	 the path to the genome used to construct the graph (can be repeated)
	--k_mer_size[-k] the size of the k-mer.
			 It must be the same value used when constructing the input graph
	--output_file_name[-o] the name of the output file
	--output_index[-oi] write the index to a binary file
	--output_graph_binary[-ogb] write the graph in binary format
	-v verbosity
[PARAMETERS]: update

	-h[--help]	 prints this help message
	--input_graph[-ig]	 the path to the pangenome graph in gfa or fasta format
	--input_genome [-igr]	 the path to the genome used to augment the input graph
	--k_mer_size[-k] the size of the k-mer.
			 It must be the same value used when constructing the input graph
	--k_mer_file	 the file of absent k-mers from the graph if already computed
	--output_file_name[-o] the name of the output file
	--update_index[-u] index the constructed funitigs
	--load_index[-li] the path to the saved index
	--output_index[-oi] write the index to a binary file
	--output_graph_binary[-ogb] write the graph in binary format
	--load_graph_binary[-lgb] load the graph from binary file
	--load_Colors[-lc] read ccdbgupdater colors
	-v verbosity
[PARAMETERS]: index

	-h[--help]	 prints this help message
	--input_graph[-ig]	 the path to the pangenome graph in gfa/fasta/binary format
	--k_mer_size[-k] the size of the k-mer
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
  
```
### Examples
  1. **Update a compacted de bruijn graph from a kmer file**
     ```
     ./ccdbupdater --input_graph graph.gfa --k_mer_file kmers.txt -k 31 -o output_prefix
     ```
     The compacted de Bruijn graph *graph.gfa* is updated  from the 31-mers (`-k 31`) of files *kmers.txt*. The updated graph is stored in *output_prefix.fa*

  2. **Update a compacted de Bruijn graph from a reference genome file**
     ```
     ./ccdbgupdater --input_graph graph.gfa --input_genome B.fa -k 31 -o output_prefix
     ```
     The compacted de Bruijn graph is updated with the absent 31-mers (`-k 31`) of file *B.fa* (`--input_genome B.fa`). Here [kmtricks](https://github.com/tlemane/kmtricks) is used to find the absent k-mers. The graph is written to file *output_prefix.fa* (`-o output_prefix`).



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
     We store (k-1)-mers of the unitig of the graph as (k-1)-mer->Vector[(unitig id,position,orientation)]

  3. **Constructing the unitigs from absent k-mers with respect to the graph**
     We go over the absent k-mers, we try to extend each by its in neighbor k-mer from left and right. This process is implemented by following
     the node centric definition of a de bruijn graph and by ensuring that the extended unitig does not branch to graph. 
     If it branches or it breaks the definition, we break the extension. We use a flag to tell if a k-mer was used before or not
    
  4. **Indexing the contructed unitigs**
    Once the unitigs are constructed, we index them by their (k-1)-mer prefixes and suffixes
  5. **Updating the input graph**
    Finally we go over the index built in 4, we determine the decision (split,merge) by evaluating the position of the prefix and its multiplicity 
    both in the unitigs of the graph and the constructed unitigs

   6. **Updating the index of the graph**
    We update the index of the graph on the fly.

## Citation
This work is on progress, but it relies on kmtricks. If you use it cite kmtricks:
```
@article {10.1093/bioadv/vbac029,
    author = {Lemane, Téo and Medvedev, Paul and Chikhi, Rayan and Peterlongo, Pierre},
    title = "{kmtricks: efficient and flexible construction of Bloom filters for large sequencing data collections}",
    journal = {Bioinformatics Advances},
    volume = {2},
    number = {1},
    year = {2022},
    month = {04},
}
```

## Contact

For any question, feedback or problem, please feel free to file an issue on this GitHub repository and we will get back to you as soon as possible.
