#include "stdio.h"
#include "cmath"

size_t baseToInt(char base);
uint64_t hash(uint64_t key);
uint64_t hash(string kmer);
void createHashTable(ifstream& readStructFile,vector<uint64_t> & hashes);