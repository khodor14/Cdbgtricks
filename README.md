[INTRODUCTION](#ccdbgupdater)   
[INSTALLATION](#installation)   
[USAGE](#binary-usage)  
[EXAMPLES](#examples)   
[ALGORITHM](#algorithm)   
[CITATION](#citation)   
[CONTACT](#contact)

# Cdbgtricks
Cdbgtricks is a C++11 modular tool, mainly designed to update a compacted de Bruijn graph, when adding nex sequences. 
In addition it indexes the graph, updates the index when adding sequences and so it enables querying sequences on the graph. 


## Requirements

* C++11 compiler: [GCC](https://gcc.gnu.org/) >= 4.8.5

* [kmtricks](https://github.com/tlemane/kmtricks) (to find the set of absent k-mers which are in the genome to be added but not in the graph)

## Installation
* Installing kmtricks: First Install kmtricks [version v1.2.1]
	* From [Bioconda](https://bioconda.github.io):

  ```bash
  conda install -c conda-forge -c tlemane kmtricks #(latest version)
  ```

* Installing ccdbupdater: Currently we have only the source version
	* From sources

  ```bash
  git clone --recursive https://github.com/khodor14/Cdbgtricks.git
  cd Cdbgtricks && mkdir build && cd build
  cmake -S ../ -B .
  make
  ```

## Binary usage:
### Setup the ulimit
First change the number of files that can be oppened

  ```
  ulimit -n 2000
  ```

### Running Cdbgtricks
```
./cdbgtricks
```
displays the command line interface:
```
Usage: ./cdbgtricks [COMMAND] [PARAMETERS]
XXXTODO copy the help once stabilizedXXX
```
### Examples
  1. **Indexing a compacted de Bruijn Graph**
     ```
     ./cdbgtricks index -ig ../test/graph.gfa -k 31 -s 5000 -msb 1 -m 11 -o graph_index -it
     ```
     The 31-mers (`-k 31`) of the compacted de Bruijn graph *graph.gfa* is indexed  using a minimizer size 11 (`-m 11`). We set here rho to be 5000 (`-s 5000`).
     It means that we should have at least 5000 31-mers in a bucket sharing the same 11-mer minimizer to create an MPHF f. We set gamma to be 1 (`-msb 1`). It means that the size of the super-bucket is 5000. (`-it`) means that the graph is in text format. The output is written to *graph_index_index.bin*

  2. **Update a compacted de Bruijn graph from a reference genome file with loading the graph index**
     ```
     ./cdbgtricks update -ig ../test/graph.gfa --input_genome ../test/genome.fa -k 31 -li graph_index_index.bin -o graph_update -u -oi -t 32
     ```
     The compacted de Bruijn graph is updated with the absent 31-mers (`-k 31`) of file *genome.fa* (`--input_genome genome.fa`). Here [kmtricks](https://github.com/tlemane/kmtricks) is used to find the absent k-mers. The index is loaded (`-li`) and used to update the graph. The graph is written to file *graph_update.fa* (`-o graph_update`). The index is updated (`-u`) and written to *graph_update_index.bin*.
  3. **Update a compacted de Bruijn graph without loading its index**
     ```
     ./cdbgtricks update -ig ../test/graph.gfa --input_genome ../test/genome.fa -k 31 -o graph_update -u -oi -t 32
     ```
     The compacted de Bruijn graph is updated with the absent 31-mers (`-k 31`) of file *genome.fa* (`--input_genome genome.fa`). Here [kmtricks](https://github.com/tlemane/kmtricks) is used to find the absent k-mers. The graph is written to file *graph_update.fa* (`-o graph_update`). The index is created (`-u`) and written to *graph_update_index.bin*.
  4. **Query a set of reads (output presence/absence)**
    ```
    ./cdbgtricks query -ig ../test/graph.gfa -qr queries.fa -k 31 -r 1 -li graph_update_index.bin -o results_query
    ```
    This tries to find if the reads are present or absent. A read is present if it shares a ratio of k-mers equal to r `-r 1`.
    The output is written to *results_query.tsv*
  5. **Map a set of reads (output uni-MEMs)**
  ```
  ./cdbgtricks map -ig ../test/graph.gfa -qr unitigs_fasta.fa -k 31 -r 1 -li graph_update_index.bin -o results_mapping
  ```
  Instead of outputing if a read is present or absent, it outputs the uni-MEMs of the read. A uni-MEM is a 5-tuple (unitig id, unitig start, unitig end,read start, read end). The output is written to results_mapping.txt



## Contact

For any question, feedback or problem, please feel free to write an issue on this GitHub repository and we will get back to you as soon as possible.
