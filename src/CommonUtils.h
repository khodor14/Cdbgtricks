#include <iostream>
#include <cmath>

size_t baseToInt(char base);
uint64_t hash(uint64_t key);
uint64_t hash(std::string kmer);
char complement(char c);
void createHashTable(std::ifstream& readStructFile,std::vector<uint64_t> & hashes);
std::string reverseComplement(const std::string& s);