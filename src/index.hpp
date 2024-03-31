#ifndef index_mphf_H
#define index_mphf_H
#include <string>
#include "CommonUtils.h"
#include "ParseGFA.h"
#include <vector>
#include <ankerl/unordered_dense.h>
#include "zstr.hpp"
#include "BooPHF.h"
typedef boomphf::SingleHashFunctor<uint64_t> hasher_t;
typedef boomphf::mphf<uint64_t, hasher_t> MPHF;
struct minimizer_info {
	uint64_t mphf_size;
    uint64_t bucket_id;
    bool small;
    bool empty;
    bool created;
};
class Index_mphf
{
private:
    int k;//the k-mer length
    int m;//the minimizer length
    int log2_super_bucket;//log2 of the number of super_buckets
    size_t small_bucket_size;
    int multiplier_bucketing;//
    uint64_t super_bucket_max_size;
    uint64_t minimizer_number;
    uint64_t bucket_per_super_bucket;
    uint64_t number_of_super_buckets;
    uint64_t num_sup_buckets;//number of groups
    uint64_t smallest_super_bucket_id=0xffffffffffffffff;
    uint64_t smallest_super_bucket_size=0xffffffffffffffff;
    uint64_t smallest_bucket_id=0xffffffffffffffff;
    uint64_t smallest_bucket_size=0xffffffffffffffff;
    std::vector<minimizer_info> mphfs_info;
    std::vector<MPHF> all_mphfs;
    std::vector<std::vector<uint64_t>> position_kmers;
    template <typename T>
    void create_mphf_per_super_bucket(std::vector<uint64_t>& kmers,std::vector<uint64_t>& positions,const T& track_minimizer,uint64_t bucket_id);
    template <typename T>
    void update_super_bucket(std::vector<uint64_t>& kmers,std::vector<uint64_t>& positions,const T& track_minimizer,uint64_t bucket_id);
    void create_mphfs();
    void rearrange_positions(MPHF mphf_ref,std::vector<uint64_t> kmers,std::vector<uint64_t> positions,uint64_t bucket_id);
    void prepare_super_buckets(GfaGraph& graph);
    void read_super_file(std::string filename,std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> &super_bucket_data);
    void read_super_buckets(GfaGraph& graph,std::vector<uint64_t> &kmers,std::vector<uint64_t>& positions,
                                    std::vector<uint64_t>& track_minimizers_in_super_bucket);
    void read_super_bucket_update(GfaGraph& graph, uint64_t &num_new_kmers_new_supb,uint64_t &last_super_bucket_id,
        std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> &kmers_new_super,
        std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> &kmers_super_b_updates);
    void create_mphf(std::vector<uint64_t>& kmers,std::vector<uint64_t>& positions,GfaGraph& graph,
                    std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> superkeys,
                    std::vector<uint64_t>& track_minimizers_in_super_bucket,uint64_t &super_bucket_created);
    void update_mphfs(GfaGraph& graph, uint64_t &num_new_kmers_new_supb,uint64_t &last_super_bucket_id,
                                std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> &kmers_new_super,
                                std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> &kmers_super_b_updates,
                                std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> superkeys);
    void update_super_bucket(GfaGraph& graph,uint64_t super_bucket_id,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>> &new_super_keys);
    void update_all_super_buckets(GfaGraph& graph,std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> kmers_super_b_updates,uint64_t last_super_bucket_id);
    void compute_mphf_per_bucket(uint64_t bucket_id,std::vector<uint64_t> kmers,bool recompute);
public:
    Index_mphf()=default;
    Index_mphf(size_t k_size,size_t m_size,size_t log_super_bucket,size_t small_b_size,size_t multiplier_bucket);
    void build(GfaGraph& graph);
    void update_unitig(Unitig seq,int id,int previous_id,int starting_position,int ending_position,bool keep_orient);
    void extract_kmers_from_funitigs(std::unordered_map<int,Unitig>& constructed_unitigs,GfaGraph& graph);
    uint64_t kmer_position(uint64_t kmer);
    int get_k_length();
    uint64_t compute_minimizer_position(uint64_t kmer,uint64_t &position);
    uint64_t compute_minimizer(uint64_t kmer);
    void update_index(std::unordered_map<int,Unitig>& constructed_funitigs,GfaGraph& graph);
    void save(std::string filename);
    void load(std::string filename);
    ~Index_mphf()=default;
};
#endif // !index_mphf_H