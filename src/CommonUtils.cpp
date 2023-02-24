#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <sstream>
#include "CommonUtils.h"
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
    case 'G':
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
std::unordered_map<std::string,bool> createHashTable(std::string file_name){
    /*
    input is a file name of kmers, the output of kmtricks pipeline
    output is a hash table kmers->boolean

    the boolean will be used later to indicate weither the kmer was used in a unitig or not
    */
    std::ifstream file{file_name, std::ios::in};
    std::string line;
    std::unordered_map<std::string,bool> result;// we store here out kmers
    while(getline(file,line)){
        //getline(kmers,line);
        std::stringstream sstr {line};
        std::string next_kmer;//the kmer in the line
        std::string next_abundance;//kmer output the abundance, for now we don't use it

        sstr >> next_kmer;
        sstr >> next_abundance;
        result[getCanonical(next_kmer)]=false;
    }
    return result;
}
void write_unitigs_to_fasta(std::unordered_map<int,std::string> unitigs,std::string filename){
    /*
    This function takes the vector of unitigs as input
    It writes this unitigs to a fasta file
    */
    std::string filenameout(filename);//the name of the file
	std::ofstream out(filenameout);//create the file
	for(auto unitig:unitigs){//loop over all unitigs
		out<<">sequence"+std::to_string(unitig.first)<<std::endl;//write the header in fasta format
		out<<unitig.second<<std::endl;//write the unitig
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

bool isCanonical(const std::string& seq){
	if (seq.size() > 1){
		char first(seq[0]);
		char last(complement(seq[seq.size()-1]));
		if (first < last){
			return true;
		} else {
			if (first == last){
				std::string seq2(seq.substr(1,seq.size()-2));
				return isCanonical(seq2);
			} else {
				return false;
			}
		}
	} else {
		if (seq.size() == 1){
			switch(seq[0]){
				case 'A': return true;
				case 'C': return true;
			}
			return false;
		} else {
			return true;
		}
	}
}