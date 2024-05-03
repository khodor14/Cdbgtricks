#include "CommonUtils.h"

#include "ParseGFA.h"

#include <iostream>

#include "unitig.h"

#include <unordered_map>

#include "index.hpp"

#include <strings.h>

#include <algorithm>

#include <vector>

#include <chrono>

#include <cmath>

#include <bitset>

#include <map>

std::vector < char > possible_right_extension(uint64_t mer, int k, std::unordered_map < uint64_t, bool > & kmers_to_add_to_graph, Index_mphf & ind, GfaGraph & graph, std::unordered_map < uint64_t, bool > & candidate_splits, uint64_t * possible_join);
std::string extend_right(std::string kmer, uint64_t mer, Index_mphf & unitigs_index, std::unordered_map < uint64_t, bool > & kmers_to_add_to_graph, GfaGraph & graph, std::unordered_map < uint64_t, bool > & candidate_splits, std::unordered_map < int, std::pair < int, int >> & candidate_join, int id, bool extended_from_right, std::unordered_map < int, std::pair < int, int >> & funitigs_joins);
std::unordered_map < int, Unitig > construct_unitigs_from_kmer(Index_mphf & unitig_index, std::unordered_map < uint64_t, bool > & kmers, int k,
  GfaGraph & graph, std::unordered_map < uint64_t, bool > & candidate_splits,
  std::unordered_map < int, std::pair < int, int >> & candidate_join, std::unordered_map < int, std::pair < int, int >> funitigs_join);
Unitig unitig_from_kmers(uint64_t kmer, Index_mphf & unitig_index, std::unordered_map < uint64_t, bool > & kmers, GfaGraph & graph, int id, std::unordered_map < uint64_t, bool > & candidate_splits, std::unordered_map < int, std::pair < int, int >> & candidate_join, std::unordered_map < int, std::pair < int, int >> & funitigs_joins);
//std::unordered_map<int,Unitig> construct_unitigs_from_kmer(Index_mphf &unitig_index, std::unordered_map<uint64_t,bool> &kmers,int k,GfaGraph& graph);
void update_CdBG(Index_mphf & graph_index, GfaGraph & graph, std::unordered_map < int, Unitig > & constructed_unitigs, std::unordered_map < uint64_t, bool > & candidate_splits, std::unordered_map < int, std::pair < int, int >> & candidate_join, bool update_index);
void insert_funitigs(GfaGraph & graph, std::unordered_map < int, Unitig > & constructed_unitigs);