#include "index.hpp"
Index::Index(size_t k_size,size_t m_size){
    k=k_size;
    m=m_size;
    all_mphfs=new mphf_data[minimizer_number];
    position_kmers=new std::vector<uint64_t>[minimizer_number];
}
uint64_t Index::kmer_position(uint64_t kmer){
    uint64_t minimizer=compute_canonical_minimizer(kmer,k,m);//compute the minimizer of the k-mer
    return kmer_position[minimizer][all_mphfs[minimizer](canonical_bits(kmer,k))];// return the position associated to this k-mer
}
void Index::build(GfaGraph& graph){
    std::vector<std::ostream*> minimizer_outs;//write k-mers per file
	for (uint64_t i(0); i < minimizer_number; ++i) {
        all_mphfs[i].empty=true;
        all_mphfs[i].mphf_size=0;
		auto out = new zstr::ofstream("ccdbgupdater_out" + to_string(i) + ".gz");
		if (not out->good()) {
			cout << "Problem with files opening" << endl;
			exit(1);
		}
		out_files.push_back(out);
	}
    //go over unitigs
    for(auto node:graph.get_nodes()){
        uint64_t kmer=node.second.get_ith_mer(0,k);//get the i-th mer
        std::tuple<uint64_t,bool> seq_data=reverseComplementCanonical(kmer,k);//compute canonical k-mer and its orientation in graph
        uint64_t position=(uint64_t)node.first;//assign unitig id to position
        position=(position<<32)|std::get<1>(seq_data);//assign 0 as position in id and the computed orientation
        uint64_t minimizer=compute_canonical_minimizer(kmer,k,m);//compute minimizer
        *(minimizer_outs[minimizer])<<to_string(std::get<0>(seq_data))<<"\n";//write k-mer to temporary output file
        *(minimizer_outs[minimizer])<<to_string(position)<<"\n";//write k-mer position to output file
        all_mphfs[minimizer].mphf_size++;
        all_mphfs[minimizer].empty=false;
       for(int i=1;i<=node.second.unitig_length()-k+1;i++){
            kmer=node.second.get_next_mer(i_th_mer,i,k);//compute next k-mer
            seq_data=reverseComplementCanonical(i_th_mer,k);//compute canonical form of k-mer and its orientation
            position=(uint64_t)node.first;//assign unitig id to position
            position=(position<<32)|(i<<1)|std::get<1>(seq_data);//assign i as position in unitig and the computed orientation
            minimizer=compute_canonical_minimizer(kmer,k,m);//compute minimizer
            *(minimizer_outs[minimizer])<<to_string(std::get<0>(seq_data))<<"\n";//write k-mer to temporary output file
            *(minimizer_outs[minimizer])<<to_string(position)<<"\n";//write k-mer position to output file
            all_mphfs[minimizer].mphf_size++;
            all_mphfs[minimizer].empty=false;
       }

    }
	for (uint64_t i(0); i < minimizer_number; ++i) {
		*minimizer_outs[i]<<std::flush;
		delete (minimizer_outs[i]);
	}
}
void Index::create_mphfs(){
    //configuration to be used for all mphf
    pthash::build_configuration config;
    config.c = 6.0;
    config.alpha = 0.94;
    config.minimal_output = true;  // mphf
    config.verbose_output = true;
    for (uint64_t i(0); i < minimizer_number; ++i) {
        if(!all_mphfs[i].empty){
            //read the kmer file
            zstr::ifstream in("ccdbgupdater_out" + to_string(i) + ".gz");
            in.peek();
            std::vector<uint64_t> kmers;
            kmers.reserve(all_mphfs[i].mphf_size);
            std::vector<uint64_t> kmer_data;
            kmer_data.reserve(all_mphfs[i].mphf_size);
            size_t j=0;
            while(not in.eof() and in.good()){
                std::string line;
                std::getline(in,line);
                uint64_t kmer=std::stoi(line);
                std::getline(in,line);
                uint64_t kmer_position=std::stoi(line);
                kmers[j]=kmer;
                kmer_data[j]=kmer_position;             
                j++;
            }
            all_mphfs[i].kmer_MPHF->build_in_internal_memory(kmers.begin(),kmers.size(),config);//build the mphf
            /*
            we need to rearrange in place the data wrt to the values computed by the mphf
            */
           uint64_t current=0;
           size_t computed=0;
           while(computed!=kmers.size()){
                uint64_t true_position=all_mphfs[i].kmer_MPHF->(kmers[current]);
                if(true_position==current){
                    current++;
                }
                else{
                    uint64_t temp=kmers[current];
                    //swap kmers
                    kmers[current]=kmers[true_position];
                    kmers[true_position]=temp;
                    //swap kmer positions
                    temp=kmer_data[current];
                    kmer_data[current]=kmer_data[true_position];
                    kmer_data[true_position]=temp;
                }
                computed++;
           }
           //assign the rearranged k-mer positions to the minimizer
           position_kmers[i]=kmer_data;
        }
 
	}
}