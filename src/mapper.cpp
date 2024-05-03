#include "mapper.hpp"

Mapper::Mapper(float r) {
  ratio = r;
}
std::vector < int > Mapper::map_back(const std::string_view seq, Index_mphf & graph_ind, GfaGraph & graph) {
  int k = graph_ind.get_k_length();
  //as the length is less than the actual kmer, then it cannot be mapped
  if (seq.length() < k) {
    return std::vector < int > ();
  }
  int mapped = 0;
  std::vector < int > res; //store response here
  uint64_t kmer = kmer_to_bits(seq.substr(0, k)); //first kmer
  uint64_t canonical_kmer = canonical_bits(kmer, k); //its canonical form
  uint64_t kmer_pos = graph_ind.kmer_position(canonical_kmer); //its possible position in the graph
  uint64_t pref_kmer = canonical_bits(kmer >> 2, k - 1); //its prefix
  int id_unitig_graph = (int)(kmer_pos >> 32);
  if (kmer_pos != 0) { //non zero values are valid position as id cannot be zero
    if (graph.test_kmer_presence(kmer, pref_kmer, kmer_pos, k) != -1) {
      res.push_back(id_unitig_graph);
      mapped++;
    } else {
      res.push_back(0);
    }
  } else {
    //the kmer is not in the index so it is not in the graph
    res.push_back(0);
  }
  for (int i = k; i < seq.length(); i++) {
    canonical_kmer = canonical_bits(kmer, k);
    kmer_pos = graph_ind.kmer_position(canonical_kmer);
    pref_kmer = canonical_bits(kmer >> 2, k - 1);
    if (kmer_pos != 0) {
      if (graph.test_kmer_presence(kmer, pref_kmer, kmer_pos, k) != -1) {
        id_unitig_graph = (int)(kmer_pos >> 32);
        res.push_back(id_unitig_graph);
        mapped++;
      } else {
        res.push_back(0);
      }
    } else {
      res.push_back(0);
    }
  }
  if (mapped / res.size() < ratio) {
    return std::vector < int > ();
  }
  return res;
}