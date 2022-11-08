#include <CommonUtils.h>
#include <ParseGFA.h>
#include <iostream>
#include <index_kmer.h>
#include <unordered_map>
std::tuple<std::string,int,std::string,int> split_unitig(std::string unitig,int id,int new_id,int position,int k);
std::tuple<std::string,int> join(std::string unitig1,std::string unitig2,int id,int k);
std::string compact(std::string unitig,std::string kmer,bool conact_position);
void modify_unitigs(std::string kmer,Index & table,std::unordered_map<int,std::string> unitigs);
std::vector<char> possible_right_extension(std::string mer,std::unordered_map<std::string,bool> kmers_to_add_to_graph);
std::string extend_right(std::string kmer, Index& unitigs_index,std::unordered_map<std::string,bool> kmers_to_add_to_graph);
std::string unitig_from_kmers(std::string kmer,Index &unitig_index, std::unordered_map<std::string,bool> kmers);
std::vector<std::string> construct_unitigs_from_kmer(std::string kmer,Index &unitig_index, std::unordered_map<std::string,bool> kmers);

