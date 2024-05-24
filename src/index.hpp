#pragma once

#include <string>

#include "CommonUtils.h"

#include "ParseGFA.h"

#include <vector>

#include <pthash.hpp>

#include "zstr.hpp"

#include <unordered_set>

struct minimizer_info {
  uint64_t mphf_size;
  uint64_t bucket_id;
  bool small;
  bool empty;
  bool created;
  template < typename Visitor >
    void visit(Visitor & visitor) {
      visitor.visit(mphf_size);
      visitor.visit(bucket_id);
      visitor.visit(empty);
      visitor.visit(created);
    }
};
class Index_mphf {
  private: int k; //the k-mer length
  int m; //the minimizer length
  int log2_super_bucket; //log2 of the number of super_buckets
  size_t small_bucket_size;
  int multiplier_bucketing; //
  uint64_t super_bucket_max_size;
  uint64_t minimizer_number;
  uint64_t bucket_per_super_bucket;
  uint64_t number_of_super_buckets;
  uint64_t num_sup_buckets; //number of groups
  uint64_t smallest_super_bucket_id = 0xffffffffffffffff;
  uint64_t smallest_super_bucket_size = 0xffffffffffffffff;
  uint64_t smallest_bucket_id = 0xffffffffffffffff;
  uint64_t smallest_bucket_size = 0xffffffffffffffff;
  std::vector < minimizer_info > mphfs_info;
  std::vector < pthash::single_phf < pthash::murmurhash2_64,
  pthash::skew_bucketer,
  pthash::dictionary_dictionary,
  true,
  pthash::add_displacement >> all_mphfs;
  std::vector < std::vector < uint64_t >> position_kmers;
  template < typename T >
  void create_mphf_per_super_bucket(std::vector < uint64_t > & kmers, std::vector < uint64_t > & positions,
    const T & track_minimizer, uint64_t bucket_id);
  template < typename T >
  void update_super_bucket(std::vector < uint64_t > & kmers, std::vector < uint64_t > & positions,
    const T & track_minimizer, uint64_t bucket_id);
  void create_mphfs();
  void rearrange_positions(pthash::single_phf < pthash::murmurhash2_64, pthash::skew_bucketer, pthash::dictionary_dictionary, true, pthash::add_displacement > mphf_ref, std::vector < uint64_t > kmers, std::vector < uint64_t > positions, uint64_t bucket_id);
  void prepare_super_buckets(GfaGraph & graph);
  void read_super_file(std::string filename, std::unordered_map < uint64_t, std::vector < std::tuple < uint64_t, uint64_t, uint64_t >>> & super_bucket_data);
  void read_super_buckets(GfaGraph & graph, std::vector < uint64_t > & kmers, std::vector < uint64_t > & positions,
    std::vector < uint64_t > & track_minimizers_in_super_bucket);
  void read_super_bucket_update(GfaGraph & graph, uint64_t & num_new_kmers_new_supb, uint64_t & last_super_bucket_id,
    std::unordered_map < uint64_t, std::vector < std::tuple < uint64_t, uint64_t, uint64_t >>> & kmers_new_super,
    std::unordered_map < uint64_t, std::vector < std::tuple < uint64_t, uint64_t, uint64_t >>> & kmers_super_b_updates);
  void create_mphf(std::vector < uint64_t > & kmers, std::vector < uint64_t > & positions, GfaGraph & graph,
    std::unordered_map < uint64_t, std::vector < std::tuple < uint64_t, uint64_t, uint64_t >>> superkeys,
    std::vector < uint64_t > & track_minimizers_in_super_bucket, uint64_t & super_bucket_created);
  void update_mphfs(GfaGraph & graph, uint64_t & num_new_kmers_new_supb, uint64_t & last_super_bucket_id,
    std::unordered_map < uint64_t, std::vector < std::tuple < uint64_t, uint64_t, uint64_t >>> & kmers_new_super,
    std::unordered_map < uint64_t, std::vector < std::tuple < uint64_t, uint64_t, uint64_t >>> & kmers_super_b_updates,
    std::unordered_map < uint64_t, std::vector < std::tuple < uint64_t, uint64_t, uint64_t >>> superkeys);
  void update_super_bucket(GfaGraph & graph, uint64_t super_bucket_id, std::vector < std::tuple < uint64_t, uint64_t, uint64_t >> & new_super_keys);
  void update_all_super_buckets(GfaGraph & graph, std::unordered_map < uint64_t, std::vector < std::tuple < uint64_t, uint64_t, uint64_t >>> kmers_super_b_updates, uint64_t last_super_bucket_id);
  public: Index_mphf() =
    default;
  Index_mphf(size_t k_size, size_t m_size, size_t log_super_bucket, size_t small_b_size, size_t multiplier_bucket);
  void build(GfaGraph & graph);
  void update_unitig(Unitig seq, int id, int previous_id, int starting_position, int ending_position, bool keep_orient);
  void update_unitig_with_skips(Unitig seq, int id, int previous_id, std::vector<std::pair<int,int>> skips);
  void extract_kmers_from_unitig(Unitig &u,const int id,int start,int end,std::vector < std::ostream * > &minimizer_outs);
  void extract_kmers_from_funitigs_in_unitigs(GfaGraph & graph, std::unordered_map<int,std::vector<std::pair<int, int>>> funitig_offsets, std::vector < std::ostream * > &minimizer_outs);
  void extract_kmers_from_funitigs(std::unordered_map < int, Unitig > & constructed_unitigs, GfaGraph & graph, std::vector< int > & unused_ids, std::unordered_map < int, std::vector< std::pair< int, int > > > & track_funitig_offsets);
  uint64_t kmer_position(uint64_t kmer);
  uint64_t kmer_position_minimizer(uint64_t kmer, uint64_t minimizer);
  int get_k_length();
  int get_m_length();
  uint64_t compute_minimizer_position(uint64_t kmer, uint64_t & position);
  uint64_t compute_minimizer(uint64_t kmer);
  void update_index(std::unordered_map < int, Unitig > & constructed_funitigs, GfaGraph & graph, std::vector< int > & unused_ids, std::unordered_map < int, std::vector< std::pair< int, int > > > & track_funitig_offsets);
  template < typename Visitor >
  void visit(Visitor & visitor) {
    visitor.visit(k);
    visitor.visit(m);
    visitor.visit(small_bucket_size);
    visitor.visit(log2_super_bucket);
    visitor.visit(multiplier_bucketing);
    visitor.visit(bucket_per_super_bucket);
    visitor.visit(number_of_super_buckets);
    visitor.visit(minimizer_number);
    visitor.visit(super_bucket_max_size);
    visitor.visit(num_sup_buckets);
    visitor.visit(smallest_bucket_id);
    visitor.visit(smallest_bucket_size);
    visitor.visit(smallest_super_bucket_id);
    visitor.visit(smallest_super_bucket_size);
    visitor.visit(mphfs_info);
    visitor.visit(all_mphfs);
    visitor.visit(position_kmers);
  }
  ~Index_mphf() =
  default;
};
//#endif // !index_mphf_H