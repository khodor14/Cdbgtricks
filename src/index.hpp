#ifndef index_H
#define index_H
#include <string>
#include "CommonUtils.h"
#include "ParseGFA.h"
#include <vector>
#include "../external/pthash/include/pthash.hpp"
#include "../external/essentials/include/essentials.hpp"
#include "zstr.hpp"
struct mphf_data {
	uint64_t mphf_size;
	pthash::single_phf<pthash::murmurhash2_64,pthash::dictionary_dictionary,true> *kmer_MPHF;
	bool empty;
    template <typename Visitor>
    void visit(Visitor& visitor){
        visitor.visit(mphf_size);
        visitor.visit(kmer_MPHF);
        visitor.visit(empty);
    }
};
class Index
{
private:
    size_t k;//the k-mer length
    size_t m;//the minimizer length
    uint64_t minimizer_number=0xffffffffffffffff>>(64-2*m);
    mphf_data *all_mphfs;
    std::vector<uint64_t> *position_kmers;
    void create_mphfs();
public:
    Index(size_t k_size,size_t m_size);
    void build(GfaGraph& graph);
    uint64_t kmer_position(uint64_t kmer);
    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(k);
        visitor.visit(m);
        visitor.visit(all_mphfs);
        visitor.visit(position_kmers);
    }
    ~Index()=default;
};
#endif // !index_H