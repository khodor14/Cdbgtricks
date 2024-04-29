#ifndef Query_S
#define Query_S
#include "index.hpp"
#include "ParseGFA.h"
#include "unitig.h"
#include "zstr.hpp"
#include "kseq.h"
#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
KSEQ_INIT(gzFile, gzread);
#endif
class QStreamer{
    public:
        QStreamer(GfaGraph& g,Index_mphf& ind);
        bool query(uint64_t kmer,uint64_t minimizer);
        float stream(const std::string &sequence);
        void stream_all(const std::string & filename_in,const std::string & filename_out,const double &ratio);
    private:
        Index_mphf& graph_index;
        GfaGraph& graph;
        Unitig last_uni;
        uint64_t last_kmer;
        uint64_t minimizer_number;
        int k;
        int m;
        int last_id;
        int last_position;
};

#endif // !Query_S