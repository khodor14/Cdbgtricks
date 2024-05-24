#include "CommonUtils.h"
#include "ParseGFA.h"
#include <iostream>
#include "unitig.h"
#include <unordered_map>
#include "index.hpp"
#include <vector>
#include <algorithm>
#include <utility>
#include <iterator> 
std::vector<char> possible_right_extension(uint64_t mer,int k,std::unordered_map<uint64_t,bool> &kmers_to_add_to_graph,Index_mphf& ind,GfaGraph& graph,std::unordered_map<uint64_t,bool>& candidate_splits,uint64_t *possible_join);
std::string extend_right(std::string kmer,uint64_t mer, Index_mphf& unitigs_index,std::unordered_map<uint64_t,bool> &kmers_to_add_to_graph,GfaGraph& graph,std::unordered_map<uint64_t,bool>& candidate_splits,std::unordered_map<int,std::pair<int,int>> &candidate_join,int id,bool extended_from_right,std::unordered_map<int,std::pair<int,int>> &funitigs_joins);
std::unordered_map<int,Unitig> construct_unitigs_from_kmer(Index_mphf &unitig_index, std::unordered_map<uint64_t,bool> &kmers,int k,
    GfaGraph& graph,std::unordered_map<uint64_t,bool> &candidate_splits,
    std::unordered_map<int,std::pair<int,int>> &candidate_join,std::unordered_map<int,std::pair<int,int>> &funitigs_join);
Unitig unitig_from_kmers(uint64_t kmer,Index_mphf &unitig_index, std::unordered_map<uint64_t,bool> &kmers,GfaGraph& graph,int id,std::unordered_map<uint64_t,bool> &candidate_splits,std::unordered_map<int,std::pair<int,int>> &candidate_join,std::unordered_map<int,std::pair<int,int>> &funitigs_joins);
void update_CdBG(Index_mphf & graph_index, GfaGraph & graph, std::unordered_map < uint64_t, bool > & candidate_splits, 
                std::unordered_map < int, std::pair < int, int >> & candidate_join, std::unordered_map<int,int> &map_on_split_if_join,bool update_index);
void insert_funitigs(GfaGraph& graph,std::unordered_map<int,Unitig>& constructed_unitigs, std::vector< int > & unused_ids);
void merge_all(GfaGraph& graph,std::unordered_map<int,Unitig>& funitigs,Index_mphf &unitig_index,std::unordered_map<int,std::pair<int,int>> &funitig_joins,
            std::unordered_map<int,std::pair<int,int>> &unitig_joins,std::unordered_map<int,int> &right_after_split,
            std::unordered_map<int,std::vector<std::pair<int,int>>> &track_funitig_offsets,std::vector<int> &unused_id);
void merge_one_funitig(uint64_t kmmer,uint64_t position,std::unordered_map<uint64_t,uint64_t> &suffix_prefix_funitigs,
          std::unordered_map<uint64_t,uint64_t> &suffix_prefix_unitigs,std::unordered_map<int,Unitig>& funitigs,GfaGraph& graph,Index_mphf &unitig_index,
          std::unordered_map<int,std::vector<std::pair<int,int>>> &track_funitig_offsets,std::vector<int> &unused_id);
void update_all_with_skips(const std::unordered_map< int, std::vector<std::pair< int, int >>> &track_funitig_offsets, GfaGraph &graph, Index_mphf &unitig_index);
void index_concerned_funitigs(std::unordered_map<uint64_t,uint64_t> &suffix_prefix_funitigs,std::unordered_map<int,Unitig>& funitigs,std::unordered_map<int,std::pair<int,int>> &funitig_joins,int k);
void index_concerned_unitigs(std::unordered_map<uint64_t,uint64_t> &suffix_prefix_unitigs,GfaGraph& graph,std::unordered_map<int,std::pair<int,int>> &unitig_joins,
                            std::unordered_map<int,int> right_after_split,int k);
void index_suffix_prefix(const int id,const int left,const int right,Unitig &u,std::unordered_map<uint64_t,uint64_t> &suffix_prefix_funitigs,int k);
std::string concat(Unitig funitig, Unitig uni,int offset_funitig,int offset_unitig,int k,bool &order, bool &reversed);
void shift_and_merge(std::vector<std::pair<int, int>> &destination, const std::vector<std::pair<int, int>> &source, int shift);
void merge_and_track(std::string merged,std::string merged2,int unitig1_id,int unitig2_id,
                      int unitig1_length,int unitig2_length,int funitig_length,
                      Index_mphf &unitig_index,GfaGraph& graph,std::unordered_map<int,std::vector<std::pair<int,int>>> &track_funitig_offsets,
                      std::vector<int> &unused_id,std::unordered_map<uint64_t,uint64_t> &suffix_prefix_unitigs);
void insert_suff_pref(uint64_t kmmer,int id,int offset,int k,std::unordered_map<uint64_t,uint64_t> &suffix_prefix_index);
void update_graph(GfaGraph &graph, Index_mphf &unitig_index,std::unordered_map<int,Unitig> &funitigs, std::unordered_map<int , std::pair<int, int>> &unitig_joins,
                  std::unordered_map<int , std::pair<int, int>> &funitig_joins,std::unordered_map < uint64_t, bool > & candidate_splits, bool update_index);
void update_suff_pref_index(uint64_t kmmer,int id,int offset,int k,std::unordered_map<uint64_t,uint64_t> &suffix_prefix_unitigs);
void update_all_with_skips(const std::unordered_map< int, std::vector<std::pair< int, int >>> &track_funitig_offsets, GfaGraph &graph, Index_mphf &unitig_index);
std::vector<std::pair<int,int>> reverseOffsets(std::vector<std::pair<int,int>> funiOffsets,int shift);