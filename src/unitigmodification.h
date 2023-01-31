#include "CommonUtils.h"
#include "ParseGFA.h"
#include <iostream>
#include "index_kmer.h"
#include <unordered_map>
#include <sparsehash/sparse_hash_map>
std::vector<char> possible_right_extension(std::string mer,google::sparse_hash_map<std::string,bool> &kmers_to_add_to_graph);
std::string extend_right(std::string kmer, Index& unitigs_index,google::sparse_hash_map<std::string,bool> &kmers_to_add_to_graph);
std::string unitig_from_kmers(std::string kmer,Index &unitig_index, google::sparse_hash_map<std::string,bool> &kmers);
google::sparse_hash_map<int,std::string> construct_unitigs_from_kmer(Index &unitig_index, google::sparse_hash_map<std::string,bool> &kmers);
void split(google::sparse_hash_map<int,std::string> &graph_unitigs,Index &ind,int id,int position,int k,int *max_node_id,int *num_split,float *time_split,bool verbose);
int decision(const std::vector<std::tuple<int,int,bool>> &occ_graph,const std::vector<std::tuple<int,int,bool>> &occ_unitigs,const google::sparse_hash_map<int,std::string> &unitigs,int k);
void checkAndMerge(std::tuple<int,int,bool> occurence_graph,std::tuple<int,int,bool> occurence_unitig,google::sparse_hash_map<std::string,std::vector<std::tuple<int,int,bool>>> &unitig_index,Index& index_graph,google::sparse_hash_map<int,std::string> &graph_unitigs,google::sparse_hash_map<int,std::string> &constructed_unitigs,int *max_node_id,int *num_split,int *num_join,float *time_split,float * time_join,bool verbose);
std::tuple<std::string,bool> can_we_merge(int position_u,int position_g,bool orient_u,bool orient_g,std::string unitig_constrct,std::string unitig_graph,Index& index_table);
void merge_unitigs(google::sparse_hash_map<std::string,std::vector<std::tuple<int,int,bool>>>& unitig_index,Index& graph_index,google::sparse_hash_map<int,std::string>& graph_unitigs,google::sparse_hash_map<int,std::string>& constructed_unitigs,int *max_node_id,int *num_split,int *num_join,float *time_split,float * time_join,float * time_update,bool verbose,bool update_index);
google::sparse_hash_map<std::string,std::vector<std::tuple<int,int,bool>>> index_constructed_unitigs(const google::sparse_hash_map<int,std::string>& constructed_unitigs,int k);