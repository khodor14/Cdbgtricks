
#include <CommonUtils.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
size_t baseToInt(char base){
    switch (base)
    {
    case 'a':
    case 'A':
        return 0;
    case 'c':
    case 'C':
        return 1;
    case 'g':
    case "G":
        return 2;
    case 't':
    case 'T':
        return 3;
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

uint64_t hash(std::string kmer){
    std::uint64_t result=0;
    for(int i=0;i<kmer.size();i++){
        result=result<<2;
        result=result|baseToInt(kmer[i]);
    }
    return hash(result);
}
void createHashTable(std::ifstream& kmers,std::vector<uint64_t> & hashes){
    std::string line;
    while(getline(kmers,line)){
        //getline(kmers,line);
        std::stringstream sstr {line};
        std::string next_kmer;
        std::string next_abundance;

        sstr >> next_kmer;
        sstr >> next_abundance;
        hashes.push_back(hash(next_kmer));
    }
}
char complement(char c){
 switch(c){
    case 'A':
    case 'a':
        return 'T';
    case 'T':
    case 't':
        return 'A';
    case 'C':
    case 'c':
        return 'G';
    case 'G':
    case 'g':
        return 'C';
 }
}
std::string reverseComplement(const std::string& s){
	std::string rc(s.size(),0);
	for (int i((int)s.length() - 1); i >= 0; i--){
		rc[s.size()-1-i] = complement(s[i]);
	}
	return rc;
}
std::string getCanonical(const std::string& s){
    return std::min(s,reverseComplement(s));
}