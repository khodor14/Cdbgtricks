#ifndef index_mphf_H
#define index_mphf_H
#include <string>
#include "CommonUtils.h"
#include "ParseGFA.h"
#include <vector>
#include "../external/pthash/include/pthash.hpp"
#include <ankerl/unordered_dense.h>
#include "zstr.hpp"
struct mphf_data {
	uint64_t mphf_size;
	pthash::single_phf<pthash::murmurhash2_64,pthash::dictionary_dictionary,true> kmer_MPHF;
	bool empty;
    bool created;
    bool small;
    template <typename Visitor>
    void visit(Visitor& visitor){
        visitor.visit(mphf_size);
        visitor.visit(kmer_MPHF);
        visitor.visit(empty);
        visitor.visit(small);
        visitor.visit(created);
    }
};
class Index_mphf
{
private:
    int k;//the k-mer length
    int m;//the minimizer length
    size_t small_bucket_size;
    uint64_t minimizer_number;
    std::vector<mphf_data> all_mphfs;
    std::vector<std::vector<uint64_t>> position_kmers;
    std::vector<uint64_t> position_kmers_small_buckets;
    pthash::single_phf<pthash::murmurhash2_64,pthash::dictionary_dictionary,true> mphf_grouped_buckets;
    void create_mphfs();
    void rearrange_positions(pthash::single_phf<pthash::murmurhash2_64,pthash::dictionary_dictionary,true> mphf_ref,std::vector<uint64_t>& kmers,std::vector<uint64_t>& positions,uint64_t min_i,bool grouped);
public:
    Index_mphf()=default;
    Index_mphf(size_t k_size,size_t m_size,size_t small_b_size);
    void build(GfaGraph& graph);
    void update(GfaGraph& graph);
    void update_unitig(Unitig seq,int id,int previous_id,int starting_position,int ending_position,bool keep_orient);
    void extract_kmers_from_funitigs(std::unordered_map<int,Unitig>& constructed_unitigs,GfaGraph& graph);
    uint64_t kmer_position(uint64_t kmer);
    int get_k_length();
    uint64_t compute_minimizer_position(uint64_t kmer,uint64_t &position);
    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(k);
        visitor.visit(m);
        visitor.visit(small_bucket_size);
        visitor.visit(all_mphfs);
        visitor.visit(position_kmers);
        visitor.visit(mphf_grouped_buckets);
        visitor.visit(position_kmers_small_buckets);
    }
    ~Index_mphf()=default;
};
#endif // !index_mphf_H