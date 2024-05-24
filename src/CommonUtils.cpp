#include "CommonUtils.h"

const char bToN[] = {
  'A',
  'C',
  'G',
  'T'
};
const char revN[] = {
  'T',
  'G',
  'C',
  'A'
};
const uint8_t revB[] = {
  3,
  2,
  1,
  0
};
size_t baseToInt(char base) {
  switch (base) {
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
  return 0; //default value - TODO: throw exception
}
uint64_t hash(std::string kmer) {
  /*
  convert kmer to 64 bits
  */
  std::uint64_t result = 0;
  for (int i = 0; i < kmer.size(); i++) {
    result = result << 2;
    result = result | baseToInt(kmer[i]);
  }
  return result;
}
std::tuple < uint64_t, bool > reverseComplementCanonical(uint64_t kmer_bits, int k) {
  /*
      takes the bits representation of a k-mer 
      return the bits representation of its reverse complement
  */
  uint64_t rev = canonical_bits(kmer_bits, k);
  return std::tuple < uint64_t, bool > (rev, kmer_bits == rev);
}
std::string to_string(uint64_t kmer_bits, int k) {
  /*
  convert 64 bits to kmer
  */
  int i;
  char kmer[k + 1];
  uint64_t tmp = kmer_bits;
  for (i = k - 1; i >= 0; i--) {
    kmer[i] = bToN[tmp & 3]; //take the right most 2 bits which represent the right most base in the k-mer
    tmp = tmp >> 2; //shift right to get the next base in the next iteration
  }
  kmer[k] = '\0';

  return kmer;
}
void createHashTable(std::ifstream & kmers, std::vector < uint64_t > & hashes) {
  std::string line;
  while (getline(kmers, line)) {
    //getline(kmers,line);
    std::stringstream sstr {
      line
    };
    std::string next_kmer;
    std::string next_abundance;

    sstr >> next_kmer;
    sstr >> next_abundance;
    hashes.push_back(hash(next_kmer));
  }
}
std::unordered_map < uint64_t, bool > createHashTable(std::string file_name) {
  /*
  input is a file name of kmers, the output of kmtricks pipeline
  output is a hash table kmers->boolean

  the boolean will be used later to indicate weither the kmer was used in a unitig or not
  */
  std::ifstream file {
    file_name,
    std::ios::in
  };
  std::string line;
  std::unordered_map < uint64_t, bool > result; // we store here out kmers
  while (getline(file, line)) {
    //getline(kmers,line);
    std::stringstream sstr {
      line
    };
    std::string next_kmer; //the kmer in the line
    std::string next_abundance; //kmer output the abundance, for now we don't use it

    sstr >> next_kmer;
    sstr >> next_abundance;
    result[canonical_bits(hash(next_kmer), next_kmer.length())] = false;
  }
  return result;
}
void write_unitigs_to_fasta(std::unordered_map < int, Unitig > unitigs, std::string filename) {
  /*
  This function takes the vector of unitigs as input
  It writes this unitigs to a fasta file
  */
  std::string filenameout(filename); //the name of the file
  std::ofstream out(filenameout); //create the file
  for (std::pair < int, Unitig > unitig: unitigs) { //loop over all unitigs
    out << ">sequence" + std::to_string(unitig.first) << std::endl; //write the header in fasta format
    out << unitig.second.to_string() << std::endl; //write the unitig
  }
}
char complement(char c) {
  switch (c) {
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
  return 'N'; // default value todo: throw exception?
}
std::string reverseComplement(const std::string & s) {
  std::string rc(s.size(), 0);
  for (int i((int) s.length() - 1); i >= 0; i--) {
    rc[s.size() - 1 - i] = complement(s[i]);
  }
  return rc;
}
std::string getCanonical(const std::string & s) {
  return std::min(s, reverseComplement(s));
}

bool isCanonical(const std::string & seq) {
  if (seq.size() > 1) {
    char first(seq[0]);
    char last(complement(seq[seq.size() - 1]));
    if (first < last) {
      return true;
    } else {
      if (first == last) {
        std::string seq2(seq.substr(1, seq.size() - 2));
        return isCanonical(seq2);
      } else {
        return false;
      }
    }
  } else {
    if (seq.size() == 1) {
      switch (seq[0]) {
      case 'A':
        return true;
      case 'C':
        return true;
      }
      return false;
    } else {
      return true;
    }
  }
}
uint8_t bit_encoding(std::string_view seq) {
  uint8_t encode = 0;
  int i;
  for (i = 0; i < seq.length(); i++) {
    encode = encode | (baseToInt(seq[i]) << (8 - 2 * (i + 1)));
  }
  return encode;
}
std::string bits_to_seq_4(uint8_t encoding, int length) {
  char seq[length + 1];
  uint8_t tmp = encoding;
  for (int i = length - 1; i >= 0; i--) {
    seq[i] = bToN[tmp & 3];
    tmp = tmp >> 2;
  }
  seq[length] = '\0';
  return seq;
}
uint64_t reverse_complement(const uint64_t kmer,
  const int k) {
  /*
      find reverse complement using the array rev_comp
      Inspired from kmtricks implementation, changing the operand used for XOR just to adapt for {A,C,G,T} order
  */
  uint64_t res = kmer;
  res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
  res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
  res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
  res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
  res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);
  res = res ^ 0xffffffffffffffff;
  res = (res >> (2 * (32 - k)));
  return res;
}
uint64_t canonical_bits(const uint64_t kmer,
  const int k) {
  uint64_t rev = reverse_complement(kmer, k); //reverse complement
  return rev < kmer ? rev : kmer;
}
bool is_canonical(const uint64_t kmer,
  const int k) {
  return kmer < reverse_complement(kmer, k);
}
uint64_t get_next_kmer(uint64_t kmer, char c, int k) {
  /*
  takes a kmer in binary format and the next character c
  it returns the next kmer in binary
  */
  uint64_t next_kmer = kmer & (0xffffffffffffffff >> (66 - 2 * k));
  next_kmer = (next_kmer << 2) | baseToInt(c);
  return next_kmer;
}
uint64_t kmer_to_bits(std::string_view seq) {
  /*
  convert kmer to 64 bits
  */
  std::uint64_t result = 0;
  for (int i = 0; i < seq.size(); i++) {
    result = result << 2;
    result = result | baseToInt(seq[i]);
  }
  return result;
}
uint64_t revhash_min(uint64_t minimizer) {
  minimizer = ((minimizer >> 32) ^ minimizer) * 0xD6E8FEB86659FD93;
  minimizer = ((minimizer >> 32) ^ minimizer) * 0xD6E8FEB86659FD93;
  minimizer = ((minimizer >> 32) ^ minimizer);
  return minimizer;
}
uint64_t unrevhash_min(uint64_t key) {
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return key;
}