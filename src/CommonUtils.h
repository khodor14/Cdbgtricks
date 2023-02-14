#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <sparsehash/sparse_hash_map>
size_t baseToInt(char base);
uint64_t hash(uint64_t key);
uint64_t hash(std::string kmer);
char complement(char c);
void createHashTable(std::ifstream& readStructFile,std::vector<uint64_t> & hashes);
std::string reverseComplement(const std::string& s);
std::string getCanonical(const std::string& s);
bool isCanonical(const std::string& seq);
google::sparse_hash_map<std::string,bool> createHashTable(std::string file_name);
void write_unitigs_to_fasta(google::sparse_hash_map<int,std::string> unitigs,std::string filename);
std::string to_string(uint64_t kmer_bits,int k);
uint64_t reverseComplement(uint64_t kmer_bits, int k);