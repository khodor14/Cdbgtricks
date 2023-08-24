#include "CommonUtils.h"
#include "ParseGFA.h"
#include <iostream>
#include "index_kmer.h"
#include "unitig.h"
#include <unordered_map>
#include "bit_vector.hpp"
std::unordered_map<uint64_t,std::vector<std::tuple<int,int,bool>>> index_constructed_unitigs(const std::unordered_map<int,Unitig>& constructed_unitigs,int k);
std::vector<char> possible_right_extension(uint64_t mer,int k,std::unordered_map<uint64_t,bool> &kmers_to_add_to_graph);
std::string extend_right(std::string kmer,uint64_t mer, Index& unitigs_index,std::unordered_map<uint64_t,bool> &kmers_to_add_to_graph);
std::string unitig_from_kmers_string(std::string kmer,uint64_t mer,Index& unitig_index, std::unordered_map<uint64_t,bool>& kmers);
Unitig unitig_from_kmers(uint64_t kmer,Index &unitig_index, std::unordered_map<std::string,bool> &kmers);
std::unordered_map<int,Unitig> construct_unitigs_from_kmer(Index &unitig_index, std::unordered_map<uint64_t,bool> &kmers,int k);
size_t decision(const uint64_t k_1_mer,const std::vector<std::tuple<int,int,bool>> &occ_unitigs,Index& index_graph,GfaGraph& graph,int k);
void split(GfaGraph& graph,Index &ind,uint64_t kmer_data,int *max_node_id,int *num_split,float *time_split,bool verbose);
void checkAndMerge(uint64_t occurence_graph,std::tuple<int,int,bool> occurence_unitig,std::unordered_map<uint64_t,std::vector<std::tuple<int,int,bool>>> &unitig_index,Index& index_graph,GfaGraph &graph,std::unordered_map<int,Unitig> &constructed_unitigs,int *max_node_id,int *num_split,int *num_join,float *time_split,float * time_join,bool verbose);
std::tuple<std::string,bool> can_we_merge(int position_u,int position_g,bool orient_u,bool orient_g,Unitig unitig_constrct,Unitig unitig_graph,Index& index_table);
void merge_unitigs(std::unordered_map<uint64_t,std::vector<std::tuple<int,int,bool>>>& unitig_index,Index& graph_index,GfaGraph& graph,std::unordered_map<int,Unitig>& constructed_unitigs,int *max_node_id,int *num_split,int *num_join,float *time_split,float * time_join,float * time_update,bool verbose,bool update_index);
std::vector<char> possible_right_extension(uint64_t mer,int k,std::unordered_map<uint64_t,std::pair<bool,BitVector>> &kmers_to_add_to_graph);
std::string extend_right(std::string kmer,uint64_t mer, Index& unitigs_index,std::unordered_map<uint64_t,std::pair<bool,BitVector>> &kmers_to_add_to_graph);
std::string unitig_from_kmers_string(std::string kmer,uint64_t mer,Index& unitig_index, std::unordered_map<uint64_t,std::pair<bool,BitVector>>& kmers);
Unitig unitig_from_kmers(uint64_t kmer,Index &unitig_index, std::unordered_map<uint64_t,std::pair<bool,BitVector>> &kmers);
std::unordered_map<uint64_t,bool> kmers_for_update(std::string filename,GfaGraph& graph,Index &ind);
void construct_graph_from_kmers(std::unordered_map<uint64_t,std::pair<bool,BitVector>>& kmers,GfaGraph & graph,Index& ind);