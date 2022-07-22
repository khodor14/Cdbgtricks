#!/bin/bash
# to be used for calling kmtricks
# it supports for now one genome
kmer=$1
graph_link=$2
genome_link=$3

cat "G : ${graph_link}" > "fof.txt"
cat "A1 : ${genome_link}" > "query.txt"

#process the pipeline on the graph
cmd1="kmtricks pipeline --file fof.list --run-dir graph_dir  --kmer-size ${kmer} --hard-min 1 --mode kmer:pa:bin --nb-partitions 100 --cpr --threads 32"
#process the aggregate on the genome
cmd2="kmtricks filter --in-matrix graph_dir --key query.txt --output filter --hard-min 1 --out-types k --cpr-in --cpr-out -t 32"
# find the difference
cmd3="kmtricks aggregate --run-dir filter --count <A1>:kmer --format text --cpr-in --output kmers.txt"
if command -v kmtricks &> /dev/null
then 
    $cmd1
    $cmd2
    $cmd3
else
    echo "kmtricks is not installed"
fi

#clean directories

rm -r graph_dir
rm -r filter