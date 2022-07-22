
#include "CommonUtils.h"
#include <iostream>
size_t baseToInt(char base){
    switch (base)
    {
    case "a":
    case "A":
        return 0;
    case "c":
    case "C":
        return 1;
    case "g":
    case "G":
        return 2
    case "t":
    case "T":
        return 3
    default:
        break;
    }
}
// from https://naml.us/post/inverse-of-a-hash-function/
uint64_t hash(uint64_t key) {
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return key;
}

uint64_t hash(string kmer){
    unit64_t result=0;
    for(int i=0;i<kmer.size();i++){
        result=result<<2;
        result=result|baseToInt(kmer[i]);
    }
    return hash(result);
}
void createHashTable(ifstream& kmers,std::vector<uint64_t> & hashes){
    while(not kmers.eof()){
        getline(kmers,line);
        std::stringstream sstr {line};
        std::string next_kmer;
        std::string next_abundance;

        sstr >> next_kmer;
        sstr >> next_aundance;
        hashes.push_back(hash(next_kmer));
    }
}