#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
size_t baseToInt(char base);
uint64_t hash(uint64_t key);
uint64_t hash(std::string kmer);
char complement(char c);
void createHashTable(std::ifstream& readStructFile,std::vector<uint64_t> & hashes);
std::string reverseComplement(const std::string& s);
std::string getCanonical(const std::string& s);
bool isCanonical(const std::string& seq);
std::unordered_map<std::string,bool> createHashTable(std::string file_name);
void write_unitigs_constructed_to_fasta(std::vector<std::string> unitigs);