#include "CommonUtils.h"
#include "ParseGFA.h"
#include <iostream>
#include "index_kmer.h"
#include <unordered_map>
std::vector<char> possible_right_extension(std::string mer,std::unordered_map<std::string,bool> &kmers_to_add_to_graph);
std::string extend_right(std::string kmer, Index& unitigs_index,std::unordered_map<std::string,bool> &kmers_to_add_to_graph);
std::string unitig_from_kmers(std::string kmer,Index &unitig_index, std::unordered_map<std::string,bool> &kmers);
std::unordered_map<int,std::string> construct_unitigs_from_kmer(Index &unitig_index, std::unordered_map<std::string,bool> &kmers);
void split(std::unordered_map<int,std::string> &graph_unitigs,Index &ind,int id,int position);
