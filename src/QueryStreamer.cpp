#include "QueryStreamer.hpp"
#include "CommonUtils.h"
QStreamer::QStreamer(GfaGraph& g, Index_mphf& index) : graph(g), graph_index(index) {
    last_id = 0;
    last_position = -1;
    last_kmer = 0;
    k = index.get_k_length();
    m = index.get_m_length();
    minimizer_number = 1ULL << (2 * m);
}
bool QStreamer::query(uint64_t kmer,uint64_t minimizer){
    kmer=canonical_bits(kmer,k);
    uint64_t position_tuple=graph_index.kmer_position_minimizer(kmer,minimizer);
    int id_unitig_graph=(int)(position_tuple>>32);//the id of the unitig having the first occurence
    int position_unitig=(int)((position_tuple>>1)&0x7FFFFFFF);//the position of k-mer in this unitig
    if(id_unitig_graph==last_id){
        if (position_unitig==last_position+1)
        {
            last_position=last_position+1;
            last_kmer=last_uni.get_next_mer(last_kmer,last_position,k);
        }
        else if (position_unitig==last_position-1)
        {
            last_position=last_position-1;
            last_kmer=last_uni.get_prev_mer(last_kmer,last_position,k);
        }
        else if (position_unitig!=last_position)
        {
            last_position=position_unitig;
            last_kmer=last_uni.get_ith_mer(last_position,k);
        }
        return kmer==canonical_bits(last_kmer,k);
    }
    else{
        last_id=id_unitig_graph;
        last_uni=graph.get_unitig(last_id);
        last_position=position_unitig;
        last_kmer=last_uni.get_ith_mer(last_position,k);
        return kmer==canonical_bits(last_kmer,k);
    }
}
float QStreamer::stream(const std::string &sequence){
    uint64_t first_kmer=hash(sequence.substr(0,k));
    uint64_t count_present=0;
    uint64_t all_count=sequence.length()-k+1;//total number of k-mers
    uint64_t pos_mini=0;
    uint64_t minimizer=graph_index.compute_minimizer_position(first_kmer,pos_mini);
    uint64_t min_seq=hash(sequence.substr(k-m,m));
    uint64_t hash_min=unrevhash_min(minimizer);
    count_present=query(first_kmer,revhash_min(minimizer)%minimizer_number)? count_present+1:count_present;

    for(int i=1;i<sequence.length()-k+1;i++){
        first_kmer=get_next_kmer(first_kmer,sequence[i+k-1],k);
        min_seq=get_next_kmer(min_seq,sequence[i+k-1],m);
        uint64_t canon_min_seq=canonical_bits(min_seq,m);
        uint64_t new_hash=unrevhash_min(canon_min_seq);
        if(new_hash<hash_min){
            minimizer=canon_min_seq;
            hash_min=new_hash;
            pos_mini=i+k-m;
        }
        else{
            if(i>pos_mini){//the minimizer is outdated
                minimizer=graph_index.compute_minimizer_position(first_kmer,pos_mini);
                hash_min=unrevhash_min(minimizer);
                pos_mini=pos_mini+i;
            }
        }
        count_present=query(first_kmer,revhash_min(minimizer)%minimizer_number)? count_present+1:count_present;
    }
    return (1.0*count_present)/all_count;
}
void QStreamer::stream_all(const std::string & filename_in,const std::string & filename_out,const double &ratio){
    FILE* fp = fopen(filename_in.c_str(), "r");
	if(fp==0){
		std::cerr<<"Couldn't open the file "<<filename_in<<std::endl;
	}
    std::ofstream res_map_file;
	std::ostream file_out(0);
	res_map_file.open(filename_out.c_str());
	file_out.rdbuf(res_map_file.rdbuf());
	kseq_t* kseq;
    kseq = kseq_init(gzopen(filename_in.c_str(),"r"));
    file_out<<"Query ID\t"<<"Found\t"<<"ratio shared\n";
    while(kseq_read(kseq)>=0){
        float shared=stream(kseq->seq.s);
        bool found=shared>=ratio;
        file_out<<kseq->comment.s<<"\t"<<std::to_string(found)<<"\t"<<std::to_string(shared)<<std::endl;
    }
    res_map_file.close();
}