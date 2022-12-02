#!/bin/bash
# to be used for calling kmtricks
# it supports for now one genome
kmer=$1
graph_link=$2
genome_link=$3

#kmtricks run for values of k in the range [8,255]
if [ $kmer -lt 8 ] || [ $kmer -gt 255];
then
  echo "The value of k should be in the range [8,255]"
  exit 0
fi

#creating file of file and query
cmdf="touch fof.txt"
cmdq="touch query.txt"
$cmdf
$cmdq

#writing path of the genome and the graph
echo "G : $graph_link">"fof.txt"
echo "A1 : ${genome_link}" > "query.txt"

#process the pipeline on the graph
cmd1="kmtricks pipeline --file fof.txt --run-dir graph_dir  --kmer-size ${kmer} --hard-min 1 --mode kmer:pa:bin --nb-partitions 100 --cpr --threads 32"
#process the aggregate on the genome
cmd2="kmtricks filter --in-matrix graph_dir --key query.txt --output filter --hard-min 1 --out-types k --cpr-in --cpr-out -t 32"
# find the difference
cmd3="kmtricks aggregate --run-dir filter --count A1:kmer --format text --cpr-in --output kmers.txt"

#check if kmtricks is installed
if command -v kmtricks &> /dev/null
then 
    #execute kmtricks command
    $cmd1
    $cmd2
    $cmd3
    #clean directories and txt files
    rm -r graph_dir
    rm -r filter
    rm fof.txt
    rm query.txt
else
    echo "kmtricks is not installed"
    echo "please visit https://github.com/tlemane/kmtricks for a detailed tutorial on how to install it"
fi