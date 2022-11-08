#include <CommonUtils.h>
#include <ParseGFA.h>
#include <iostream>
#include <index_kmer.h>
#include <unordered_map>
std::tuple<std::string,int,std::string,int> split_unitig(std::string unitig,int id,int new_id,int position,int k);
std::tuple<std::string,int> join(std::string unitig1,std::string unitig2,int id,int k);
std::string compact(std::string unitig,std::string kmer,bool conact_position);
void modify_unitigs(std::string kmer,Index & table,std::unordered_map<int,std::string> unitigs);
std::vector<std::tuple<std::string,bool>> extend_right(std::string kmer,bool reverse, std::unordered_map<std::string,bool> kmers);
std::vector<std::tuple<std::string,bool>> extend_left(std::string kmer,bool reverse, std::unordered_map<std::string,bool> kmers);
std::vector<std::tuple<std::string,bool>> possible_extension_left(std::string kmer, std::unordered_map<std::string,bool> kmers);
std::vector<std::tuple<std::string,bool>> possible_extension_right(std::string kmer, std::unordered_map<std::string,bool> kmers);
std::string extend(std::string kmer,Index &unitig_index, std::unordered_map<std::string,bool> kmers);
std::string unitig_from_kmers(std::string kmer,Index &unitig_index, std::unordered_map<std::string,bool> kmers);
std::vector<std::string> construct_unitigs_from_kmer(std::string kmer,Index &unitig_index, std::unordered_map<std::string,bool> kmers);

