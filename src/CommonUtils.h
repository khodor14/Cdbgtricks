#ifndef COMMON_UTILS_H
#define COMMON_UTILS_H
#include <iostream>
#include <cstdint>
#include <string>
#include <unordered_map>
#include "unitig.h"
#include "bit_vector.hpp"
#include "ParseGFA.h"
#include "index_kmer.h"
size_t baseToInt(char base);
uint64_t hash(std::string kmer);
uint64_t hash(std::string_view kmer);
char complement(char c);
void createHashTable(std::ifstream& readStructFile,std::vector<uint64_t> & hashes);
std::string reverseComplement(const std::string& s);
std::string getCanonical(const std::string& s);
bool isCanonical(const std::string& seq);
std::unordered_map<uint64_t,bool> createHashTable(std::string file_name);
void write_unitigs_to_fasta(std::unordered_map<int,Unitig> unitigs,std::string filename);
std::string to_string(uint64_t kmer_bits,int k);
std::tuple<uint64_t,bool> reverseComplementCanonical(uint64_t kmer_bits, int k);
uint8_t bit_encoding(std::string_view seq);
std::string bits_to_seq_4(uint8_t encoding,int length);
uint64_t reverse_complement(const uint64_t kmer,const int k);
uint64_t canonical_bits(const uint64_t kmer,const int k);
bool is_canonical(const uint64_t kmer,const int k);
void insert_kmers_to_map(std::string filename,std::unordered_map<uint64_t,std::pair<bool,BitVector>> & kmers,size_t num_bits,int color,int k);
std::unordered_map<uint64_t,std::pair<bool,BitVector>> prepare_kmers_for_construction(std::vector<std::string> filenames,int k);
#endif // COMMON_UTILS_H
