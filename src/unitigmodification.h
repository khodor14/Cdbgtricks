#include "CommonUtils.h"
#include "ParseGFA.h"
#include <iostream>
#include "index_kmer.h"
#include "unitig.h"
#include <unordered_map>
std::unordered_map<uint64_t,std::vector<std::tuple<int,int,bool>>> index_constructed_unitigs(const std::unordered_map<int,Unitig>& constructed_unitigs,int k);
std::vector<char> possible_right_extension(std::string mer,std::unordered_map<std::string,bool> &kmers_to_add_to_graph);
std::string extend_right(std::string kmer, Index& unitigs_index,std::unordered_map<std::string,bool> &kmers_to_add_to_graph);
std::string unitig_from_kmers(std::string kmer,Index& unitig_index, std::unordered_map<std::string,bool>& kmers);
Unitig unitig_from_kmers(uint64_t kmer,Index &unitig_index, std::unordered_map<std::string,bool> &kmers);
std::unordered_map<int,Unitig> construct_unitigs_from_kmer(Index &unitig_index, std::unordered_map<uint64_t,bool> &kmers,int k);
int decision(const std::vector<std::tuple<int,int,bool>> &occ_graph,const std::vector<std::tuple<int,int,bool>> &occ_unitigs,const std::unordered_map<int,Unitig> &unitigs,int k);
void split(std::unordered_map<int,Unitig> &graph_unitigs,Index &ind,int id,int position,int k,int *max_node_id,int *num_split,float *time_split,bool verbose);
void checkAndMerge(std::tuple<int,int,bool> occurence_graph,std::tuple<int,int,bool> occurence_unitig,std::unordered_map<uint64_t,std::vector<std::tuple<int,int,bool>>> &unitig_index,Index& index_graph,std::unordered_map<int,Unitig> &graph_unitigs,std::unordered_map<int,Unitig> &constructed_unitigs,int *max_node_id,int *num_split,int *num_join,float *time_split,float * time_join,bool verbose);
std::tuple<std::string,bool> can_we_merge(int position_u,int position_g,bool orient_u,bool orient_g,Unitig unitig_constrct,Unitig unitig_graph,Index& index_table);
void merge_unitigs(std::unordered_map<uint64_t,std::vector<std::tuple<int,int,bool>>>& unitig_index,Index& graph_index,std::unordered_map<int,Unitig>& graph_unitigs,std::unordered_map<int,Unitig>& constructed_unitigs,int *max_node_id,int *num_split,int *num_join,float *time_split,float * time_join,float * time_update,bool verbose,bool update_index);