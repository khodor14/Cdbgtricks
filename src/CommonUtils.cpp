#include <fstream>
#include <cstring>
#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include "CommonUtils.h"
#include <cassert>
const char bToN[]={'A','C','G','T'};
const char revN[]={'T','G','C','A'};
/*
rev_comp is the 8 bit encoding of reverse complement
rev_comp[rev_comp[i]] is the reverse complement of 4-mer at i
for i=0 the 4-mer is 0xff(TTTT) just reading rev_comp[0xff] gives us 0x00 (000000) which is the reverse_complement of TTTT
for i=1(AAAC=00000001) the 4-mer is 0xbf (GTTTT=10111111)
*/
const uint8_t rev_comp[256] = {0xff,0xbf,0x7f,0x3f,0xef,0xaf,0x6f,0x2f,0xdf,0x9f,0x5f,0x1f,0xcf,0x8f,0x4f,0x0f,0xfb,0xbb,0x7b,0x3b,0xeb,
                               0xab,0x6b,0x2b,0xdb,0x9b,0x5b,0x1b,0xcb,0x8b,0x4b,0x0b,0xf7,0xb7,0x77,0x37,0xe7,0xa7,0x67,0x27,0xd7,0x97,
                               0x57,0x17,0xc7,0x87,0x47,0x07,0xf3,0xb3,0x73,0x33,0xe3,0xa3,0x63,0x23,0xd3,0x93,0x53,0x13,0xc3,0x83,0x43,
                               0x03,0xfe,0xbe,0x7e,0x3e,0xee,0xae,0x6e,0x2e,0xde,0x9e,0x5e,0x1e,0xce,0x8e,0x4e,0x0e,0xfa,0xba,0x7a,0x3a,
                               0xea,0xaa,0x6a,0x2a,0xda,0x9a,0x5a,0x1a,0xca,0x8a,0x4a,0x0a,0xf6,0xb6,0x76,0x36,0xe6,0xa6,0x66,0x26,0xd6,
                               0x96,0x56,0x16,0xc6,0x86,0x46,0x06,0xf2,0xb2,0x72,0x32,0xe2,0xa2,0x62,0x22,0xd2,0x92,0x52,0x12,0xc2,0x82,
                               0x42,0x02,0xfd,0xbd,0x7d,0x3d,0xed,0xad,0x6d,0x2d,0xdd,0x9d,0x5d,0x1d,0xcd,0x8d,0x4d,0x0d,0xf9,0xb9,0x79,
                               0x39,0xe9,0xa9,0x69,0x29,0xd9,0x99,0x59,0x19,0xc9,0x89,0x49,0x09,0xf5,0xb5,0x75,0x35,0xe5,0xa5,0x65,0x25,
                               0xd5,0x95,0x55,0x15,0xc5,0x85,0x45,0x05,0xf1,0xb1,0x71,0x31,0xe1,0xa1,0x61,0x21,0xd1,0x91,0x51,0x11,0xc1,
                               0x81,0x41,0x01,0xfc,0xbc,0x7c,0x3c,0xec,0xac,0x6c,0x2c,0xdc,0x9c,0x5c,0x1c,0xcc,0x8c,0x4c,0x0c,0xf8,0xb8,
                               0x78,0x38,0xe8,0xa8,0x68,0x28,0xd8,0x98,0x58,0x18,0xc8,0x88,0x48,0x08,0xf4,0xb4,0x74,0x34,0xe4,0xa4,0x64,
                               0x24,0xd4,0x94,0x54,0x14,0xc4,0x84,0x44,0x04,0xf0,0xb0,0x70,0x30,0xe0,0xa0,0x60,0x20,0xd0,0x90,0x50,0x10,
                               0xc0,0x80,0x40,0x00};
const uint8_t revB[]={3,2,1,0};
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
    /*
    convert kmer to 64 bits
    */
    std::uint64_t result=0;
    for(int i=0;i<kmer.size();i++){
        result=result<<2;
        result=result|baseToInt(kmer[i]);
    }
    return result;
}
std::tuple<uint64_t,bool> reverseComplementCanonical(uint64_t kmer_bits, int k){
    /*
        takes the bits representation of a k-mer 
        return the bits representation of its reverse complement
    */
   uint64_t rev=canonical_bits(kmer_bits,k);
   return std::tuple<uint64_t,bool>(rev,kmer_bits==rev);
}
std::string to_string(uint64_t kmer_bits,int k){
    /*
    convert 64 bits to kmer
    */
    int i;
    char kmer[k+1];
    uint64_t tmp=kmer_bits;
    for(i=k-1;i>=0;i--){
        kmer[i]=bToN[tmp&3];//take the right most 2 bits which represent the right most base in the k-mer
        tmp=tmp>>2; //shift right to get the next base in the next iteration
    }
    kmer[k]='\0';

    return kmer;
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
std::unordered_map<uint64_t,bool> createHashTable(std::string file_name){
    /*
    input is a file name of kmers, the output of kmtricks pipeline
    output is a hash table kmers->boolean

    the boolean will be used later to indicate weither the kmer was used in a unitig or not
    */
    std::ifstream file{file_name, std::ios::in};
    std::string line;
    std::unordered_map<uint64_t,bool> result;// we store here out kmers
    while(getline(file,line)){
        //getline(kmers,line);
        std::stringstream sstr {line};
        std::string next_kmer;//the kmer in the line
        std::string next_abundance;//kmer output the abundance, for now we don't use it

        sstr >> next_kmer;
        sstr >> next_abundance;
        result[std::get<0>(reverseComplementCanonical(hash(next_kmer),next_kmer.length()))]=false;
    }
    return result;
}
void write_unitigs_to_fasta(std::unordered_map<int,Unitig> unitigs,std::string filename){
    /*
    This function takes the vector of unitigs as input
    It writes this unitigs to a fasta file
    */
    std::string filenameout(filename);//the name of the file
	std::ofstream out(filenameout);//create the file
	for(std::pair<int, Unitig> unitig:unitigs){//loop over all unitigs
		out<<">sequence"+std::to_string(unitig.first)<<std::endl;//write the header in fasta format
		out<<unitig.second.to_string()<<std::endl;//write the unitig
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
uint8_t bit_encoding(std::string_view seq){
    assert(seq.length()<=4);
    uint8_t encode=0;
    int i;
    for(i=0;i<seq.length();i++){
        encode=encode|(baseToInt(seq[i])<<(8-2*(i+1)));
    }
    return encode;
}
std::string bits_to_seq_4(uint8_t encoding,int length){
    char seq[length+1];
    uint8_t tmp=encoding;
    for(int i=length-1;i>=0;i--){
        seq[i]=bToN[tmp&3];
        tmp=tmp>>2;
    }
    seq[length]='\0';
    return seq;
}
uint64_t reverse_complement(const uint64_t kmer,const int k){
    /*
        find reverse complement using the array rev_comp
    */
    uint64_t res=kmer;
    res=(res&0x00ffffffffffff00)|rev_comp[res>>56]|(((uint64_t)rev_comp[res&0x00000000000000ff])<<56);//reversing left most 4-mer and right most 4-mer
    res=(res&0xff00ffffffff00ff)|(((uint64_t)rev_comp[(res>>48)&0x00000000000000ff])<<8)|(((uint64_t)rev_comp[(res&0x000000000000ff00)>>8])<<48);//second left most 4-mer and second last most 4-mer
    res=(res&0xffff00ffff00ffff)|(((uint64_t)rev_comp[(res>>40)&0x00000000000000ff])<<16)|(((uint64_t)rev_comp[(res&0x0000000000ff0000)>>16])<<40);//same 
    res=(res&0xffffff0000ffffff)|(((uint64_t)rev_comp[(res>>32)&0x00000000000000ff])<<24)|(((uint64_t)rev_comp[(res&0x00000000ff000000)>>24])<<32);//reversing the 8-mer in the middle
    return res>>(64-2*k);//shift to encode only k-mer
}
uint64_t canonical_bits(const uint64_t kmer,const int k){
    uint64_t rev=reverse_complement(kmer,k);//reverse complement
    return rev<kmer? rev:kmer;
}
bool is_canonical(const uint64_t kmer,const int k){
    return kmer<reverse_complement(kmer,k);
}