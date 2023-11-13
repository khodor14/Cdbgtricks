#include "index.hpp"
#include <unordered_set>
Index_mphf::Index_mphf(size_t k_size,size_t m_size,size_t log_super_buckets,size_t small_b_size){
    k=k_size;//k-mer size
    m=m_size;//minimizer size
    log2_super_bucket=log_super_buckets;
    number_of_super_buckets=1ULL<<(2*log2_super_bucket);
    bucket_per_super_bucket=1ULL<<(2*(m-log2_super_bucket));
    minimizer_number=1ULL<<(2*m);
    all_mphfs.resize(minimizer_number);
    position_kmers.resize(minimizer_number);
    small_bucket_size=small_b_size;
}
void Index_mphf::prepare_super_buckets(GfaGraph& graph){
    std::vector<std::ostream*> minimizer_outs;//write k-mers per file
	for (uint64_t i(0); i < number_of_super_buckets; ++i) {
		auto out = new zstr::ofstream("ccdbgupdater_out" + std::to_string(i) + ".gz");
		if (not out->good()) {
			std::cout << "Problem with files opening" <<std::endl;
			exit(1);
		}
		minimizer_outs.push_back(out);
	}

    for(auto node:graph.get_nodes()){
        uint64_t pos_min=0;
        uint64_t kmer=node.second.get_ith_mer(0,k);//get the i-th mer
        uint64_t prev_kmer=kmer;
        uint64_t min_seq=kmer&((0xffffffffffffffff<<(64-2*m))>>(64-2*m));//
        uint64_t canon_min_seq;
        int pos_min_seq=k-m;
        std::tuple<uint64_t,bool> seq_data=reverseComplementCanonical(kmer,k);//compute canonical k-mer and its orientation in graph
        uint64_t position=(uint64_t)node.first;//assign unitig id to position
        position=position<<32;//assign 0 as position in id and the computed orientation
        uint64_t minimizer=compute_minimizer_position(std::get<0>(seq_data),pos_min);//minimizer of kmer
        uint64_t old_minimizer=minimizer;//old minimizer for comparison only
        int counter=1;
        if(!std::get<1>(seq_data)){
            pos_min=(k-m-pos_min)%(k-m+1);
        }
        //when we have one kmer in the seq, write the data to the super-bucket
        if(node.second.unitig_length()==k){
            *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(old_minimizer)<<"\n";//write minimizer to super bucket file
            *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(prev_kmer)<<"\n";//write k-mer position to output file
            *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(position)<<"\n";//write k-mer position to output file
            *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(counter)<<"\n";//write k-mer position to output file
            all_mphfs[old_minimizer].mphf_size+=counter;
            all_mphfs[old_minimizer].empty=false;   
        }
        for(int i=1;i<node.second.unitig_length()-k+1;i++){
            pos_min_seq++;
            kmer=node.second.get_next_mer(kmer,i,k);//compute next k-mer
            seq_data=reverseComplementCanonical(kmer,k);//compute canonical form of k-mer and its orientation
            min_seq=node.second.get_next_mer(min_seq,pos_min_seq,m);
            canon_min_seq=canonical_bits(min_seq,m);
            if(canon_min_seq<minimizer){
                minimizer=canon_min_seq;
                pos_min=i+k-m;
            }
            else{
                if(i>pos_min){
                    minimizer=compute_minimizer_position(std::get<0>(seq_data),pos_min);
                    if(!std::get<1>(seq_data)){
                        pos_min=(k-m-pos_min)%(k-m+1);
                    }
                    pos_min+=i;
                }
            }
            if(minimizer==old_minimizer){
                counter++;
                if(i==node.second.unitig_length()-k){//when same minimizer and last kmer write everything to the superbucket
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(old_minimizer)<<"\n";//write minimizer to super bucket file
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(prev_kmer)<<"\n";//write k-mer position to output file
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(position)<<"\n";//write k-mer position to output file
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(counter)<<"\n";//write k-mer position to output file
                    all_mphfs[old_minimizer].mphf_size+=counter;
                    all_mphfs[old_minimizer].empty=false;
                }
            }
            else{
                *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(old_minimizer)<<"\n";//write minimizer to super bucket file
                *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(prev_kmer)<<"\n";//write k-mer position to output file
                *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(position)<<"\n";//write k-mer position to output file
                *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(counter)<<"\n";//write k-mer position to output file
                all_mphfs[old_minimizer].mphf_size+=counter;
                all_mphfs[old_minimizer].empty=false;
                counter=1;
                prev_kmer=kmer;
                position=(uint64_t)node.first;//assign unitig id to position
                position=(position<<32)|i;//assign i as position in unitig and the computed orientation
                old_minimizer=minimizer;
                if(i==node.second.unitig_length()-k){//when same minimizer and last kmer write everything to the superbucket
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(minimizer)<<"\n";//write minimizer to super bucket file
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(kmer)<<"\n";//write k-mer position to output file
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(position)<<"\n";//write k-mer position to output file
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(counter)<<"\n";//write k-mer position to output file
                    all_mphfs[old_minimizer].mphf_size+=counter;
                    all_mphfs[old_minimizer].empty=false;
                }
            }
       }
    }
	for (uint64_t i(0); i < number_of_super_buckets; ++i) {
		*minimizer_outs[i]<<std::flush;
		delete (minimizer_outs[i]);
	}
}
uint64_t Index_mphf::kmer_position(uint64_t kmer){
    uint64_t pos=0;
    uint64_t minimizer=compute_minimizer_position(kmer,pos);//compute the minimizer of the k-mer
    if(all_mphfs[minimizer].empty){
        //0 is undefined k-mer position because the position is a triplet (unitig id,pos in unitig,orientation) represented in 64 bits
        // as id is not zero then the defined position cannot be zero
        return 0;//the mphf related to the minimizer of the k-mer is not yet created so the k-mer does not exist in the graph
    }
    else if(all_mphfs[minimizer].small)
    {   
        return position_kmers_small_buckets[mphf_grouped_buckets(kmer)];//this is a singleton mphf
    }
    return position_kmers[minimizer][all_mphfs[minimizer].kmer_MPHF(kmer)];//return the position associated to this k-mer
}
void Index_mphf::build(GfaGraph& graph){
    for(uint64_t i=0;i<minimizer_number;i++){
        all_mphfs[i].mphf_size=0;
        all_mphfs[i].empty=true;
    }
    std::vector<uint64_t> kmers;
    std::vector<uint64_t> positions;
    //create files for each super bucket
    prepare_super_buckets(graph);
    //read and create mphf if the bucket is not small
    //along the way group ungrouped mphf
    read_super_buckets(graph,kmers,positions,false);
    //build the grouped mphf
    pthash::build_configuration config;
    config.c = 6.0;
    config.alpha = 0.94;
    config.minimal_output = true;  // mphf
    config.verbose_output = false;
    config.seed=1234567890;
    mphf_grouped_buckets.build_in_internal_memory(kmers.begin(),kmers.size(),config);
    position_kmers_small_buckets.resize(kmers.size());
    //rearrange the positions of kmers in grouped mphf
    for(uint64_t i(0);i<kmers.size();i++){
        position_kmers_small_buckets[mphf_grouped_buckets(kmers[i])]=positions[i];//
    }
}
void Index_mphf::read_super_buckets(GfaGraph& graph,std::vector<uint64_t> &kmers,std::vector<uint64_t>& positions,bool update_phase){
    /*
    takes the graph, an empty vector of kmers and empty vector of their position
    a boolean update phase: this function is used for construction and update
    */
    for(uint64_t super_bucket=0;super_bucket<number_of_super_buckets;super_bucket++){
        zstr::ifstream in("ccdbgupdater_out" + std::to_string(super_bucket) + ".gz");
        in.peek();
        std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> super_bucket_data;
        while(not in.eof() and in.good()){
            std::string line;
            std::getline(in,line);
            if(line.length()==0){
                break;
            }
            uint64_t minimizer=std::stoll(line);
            std::getline(in,line);
            uint64_t kmer=std::stoll(line);
            std::getline(in,line);
            uint64_t position=std::stoll(line);
            std::getline(in,line);
            int counter=std::stoi(line);
            if(super_bucket_data.count(minimizer)==0){
                super_bucket_data[minimizer]=std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>();
            }
            super_bucket_data[minimizer].push_back(std::tuple<uint64_t,uint64_t,uint64_t>(kmer,position,counter));
        }
        std::remove(("ccdbgupdater_out" + std::to_string(super_bucket) + ".gz").c_str());
        if(!update_phase){
            //construction phase
            create_mphf(kmers,positions,graph,super_bucket_data);
        }
        else{
            //update the mphfs
            update_mphfs(kmers,positions,graph,super_bucket_data);
        }
    }
}
void Index_mphf::update_mphfs(std::vector<uint64_t>& kmers,std::vector<uint64_t>& positions,GfaGraph& graph,
                                std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> superkeys){
    pthash::build_configuration config;
    config.c = 6.0;
    config.alpha = 0.94;
    config.minimal_output = true;  // mphf
    config.verbose_output = false;
    for(auto minimizer_data:superkeys){
        if(all_mphfs[minimizer_data.first].small||all_mphfs[minimizer_data.first].mphf_size<small_bucket_size){
            all_mphfs[minimizer_data.first].small=true;
            all_mphfs[minimizer_data.first].empty=false;
            for(auto super_key_data:minimizer_data.second){
                uint64_t k_mer=std::get<0>(super_key_data);
                uint64_t pos_kmer=std::get<1>(super_key_data);
                uint64_t counter=std::get<2>(super_key_data);
                std::tuple<uint64_t,bool> seq_data=reverseComplementCanonical(k_mer,k);
                uint64_t id=pos_kmer>>32;
                uint64_t pos_in_unitig=pos_kmer&0xffffffff;
                uint64_t full_position=(id<<32)|(pos_in_unitig<<1)|std::get<1>(seq_data);
                Unitig seq=graph.get_unitig(id);
                kmers.push_back(std::get<0>(seq_data));
                positions.push_back(full_position);
                for(uint64_t i=1;i<counter;i++){
                    k_mer=seq.get_next_mer(k_mer,pos_in_unitig+i,k);
                    seq_data=reverseComplementCanonical(k_mer,k);
                    full_position=(id<<32)|((pos_in_unitig+i)<<1)|std::get<1>(seq_data);
                    kmers.push_back(std::get<0>(seq_data));
                    positions.push_back(full_position);
                }
            }
        }
        else{
            std::vector<uint64_t> minimizer_kmers;
            std::vector<uint64_t> kmer_positions;
            for(auto super_key_data:minimizer_data.second){
                uint64_t k_mer=std::get<0>(super_key_data);
                uint64_t pos_kmer=std::get<1>(super_key_data);
                uint64_t counter=std::get<2>(super_key_data);
                std::tuple<uint64_t,bool> seq_data=reverseComplementCanonical(k_mer,k);
                uint64_t id=pos_kmer>>32;
                uint64_t pos_in_unitig=pos_kmer&0xffffffff;
                uint64_t full_position=(id<<32)|(pos_in_unitig<<1)|std::get<1>(seq_data);
                Unitig seq=graph.get_unitig(id);
                minimizer_kmers.push_back(std::get<0>(seq_data));
                kmer_positions.push_back(full_position);
                for(uint64_t i=1;i<counter;i++){
                    k_mer=seq.get_next_mer(k_mer,pos_in_unitig+i,k);
                    seq_data=reverseComplementCanonical(k_mer,k);
                    full_position=(id<<32)|((pos_in_unitig+i)<<1)|std::get<1>(seq_data);
                    minimizer_kmers.push_back(std::get<0>(seq_data));
                    kmer_positions.push_back(full_position);
                }
            }
            if(!all_mphfs[minimizer_data.first].empty){
                for(uint64_t pos_kmer:position_kmers[minimizer_data.first]){
                    kmer_positions.push_back(pos_kmer);
                    uint64_t km=graph.get_kmer(pos_kmer,k);
                    minimizer_kmers.push_back(canonical_bits(km,k));
                }
            }
            config.num_buckets=4*minimizer_kmers.size();
            all_mphfs[minimizer_data.first].empty=false;
            all_mphfs[minimizer_data.first].kmer_MPHF.build_in_internal_memory(minimizer_kmers.begin(),minimizer_kmers.size(),config);
            rearrange_positions(all_mphfs[minimizer_data.first].kmer_MPHF,minimizer_kmers,kmer_positions,minimizer_data.first,false);
        }
    }
}
void Index_mphf::create_mphf(std::vector<uint64_t>& kmers,std::vector<uint64_t>& positions,GfaGraph& graph,
                                std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> superkeys){
    /*
        takes a vector of kmers to which the kmers from small buckets will be added
              a vector of positions to which the positions of kmers from small buckets will be added
              A graph to retrieve the kmers associated with the super keys
              a map whose keys are the minimizers and values is a vector of tuples (kmer, position,counter) where kmer is the first kmer of the superkmer
                                                                    position is the position of this kmer, and counter is the length of super kmer (number of kmer)
        if the bucket is small => create the mphf
            else: group the kmers from small buckets
    */
    //pthash configuration
    pthash::build_configuration config;
    config.c = 6.0;
    config.alpha = 0.94;
    config.minimal_output = true;  // mphf
    config.verbose_output = false;
    config.seed=1234567890;
    for(auto minimizer_data:superkeys){
        if(all_mphfs[minimizer_data.first].mphf_size<small_bucket_size){
            all_mphfs[minimizer_data.first].small=true;
            for(auto super_key_data:minimizer_data.second){
                uint64_t k_mer=std::get<0>(super_key_data);
                uint64_t pos_kmer=std::get<1>(super_key_data);
                uint64_t counter=std::get<2>(super_key_data);
                std::tuple<uint64_t,bool> seq_data=reverseComplementCanonical(k_mer,k);
                uint64_t id=pos_kmer>>32;
                uint64_t pos_in_unitig=pos_kmer&0xffffffff;
                uint64_t full_position=(id<<32)|(pos_in_unitig<<1)|std::get<1>(seq_data);
                Unitig seq=graph.get_unitig(id);
                kmers.push_back(std::get<0>(seq_data));
                positions.push_back(full_position);
                for(uint64_t i=1;i<counter;i++){
                    k_mer=seq.get_next_mer(k_mer,pos_in_unitig+i,k);
                    seq_data=reverseComplementCanonical(k_mer,k);
                    full_position=(id<<32)|((pos_in_unitig+i)<<1)|std::get<1>(seq_data);
                    kmers.push_back(std::get<0>(seq_data));
                    positions.push_back(full_position);
                }
            }
        }
        else{
            std::vector<uint64_t> minimizer_kmers;
            std::vector<uint64_t> kmer_positions;
            for(auto super_key_data:minimizer_data.second){
                uint64_t k_mer=std::get<0>(super_key_data);
                uint64_t pos_kmer=std::get<1>(super_key_data);
                uint64_t counter=std::get<2>(super_key_data);
                std::tuple<uint64_t,bool> seq_data=reverseComplementCanonical(k_mer,k);
                uint64_t id=pos_kmer>>32;
                uint64_t pos_in_unitig=pos_kmer&0xffffffff;
                uint64_t full_position=(id<<32)|(pos_in_unitig<<1)|std::get<1>(seq_data);
                Unitig seq=graph.get_unitig(id);
                minimizer_kmers.push_back(std::get<0>(seq_data));
                kmer_positions.push_back(full_position);
                for(uint64_t i=1;i<counter;i++){
                    k_mer=seq.get_next_mer(k_mer,pos_in_unitig+i,k);
                    seq_data=reverseComplementCanonical(k_mer,k);
                    full_position=(id<<32)|((pos_in_unitig+i)<<1)|std::get<1>(seq_data);
                    minimizer_kmers.push_back(std::get<0>(seq_data));
                    kmer_positions.push_back(full_position);
                }
            }
            all_mphfs[minimizer_data.first].kmer_MPHF.build_in_internal_memory(minimizer_kmers.begin(),minimizer_kmers.size(),config);
            rearrange_positions(all_mphfs[minimizer_data.first].kmer_MPHF,minimizer_kmers,kmer_positions,minimizer_data.first,false);
        }
    }
}
void Index_mphf::update_index(std::unordered_map<int,Unitig>& constructed_funitigs,GfaGraph& graph){
    std::vector<uint64_t> kmers;
    std::vector<uint64_t> positions;
    //create files for each super bucket and write kmers from funitigs
    //
    extract_kmers_from_funitigs(constructed_funitigs,graph);
    //read and create mphf if the bucket is not small
    //along the way group ungrouped mphf
    read_super_buckets(graph,kmers,positions,true);
    //build the grouped mphf
    if(kmers.size()>0){
        pthash::build_configuration config;
        config.c = 6.0;
        config.alpha = 0.94;
        config.minimal_output = true;  // mphf
        config.verbose_output = false;
        config.seed=1234567890;
        for(uint64_t pos_kmer:position_kmers_small_buckets){
            positions.push_back(pos_kmer);
            kmers.push_back(graph.get_kmer(pos_kmer,k));
        }
        mphf_grouped_buckets.build_in_internal_memory(kmers.begin(),kmers.size(),config);
        position_kmers_small_buckets.resize(kmers.size());
        //rearrange the positions of kmers in grouped mphf
        for(uint64_t i(0);i<kmers.size();i++){
            position_kmers_small_buckets[mphf_grouped_buckets(kmers[i])]=positions[i];//
        } 
    }
}
void Index_mphf::extract_kmers_from_funitigs(std::unordered_map<int,Unitig>& constructed_unitigs,GfaGraph& graph){
    std::vector<std::ostream*> minimizer_outs(number_of_super_buckets);//write k-mers per file
	for (uint64_t i(0); i < number_of_super_buckets; ++i) {
		auto out = new zstr::ofstream("ccdbgupdater_out" + std::to_string(i) + ".gz");
		if (not out->good()) {
			std::cout << "Problem with files opening" <<std::endl;
			exit(1);
		}
		minimizer_outs[i]=out;
	}
    int id=graph.get_max_node_id();
    //go over unitigs
    for(auto node:constructed_unitigs){
        id++;
        uint64_t pos_min=0;
        uint64_t kmer=node.second.get_ith_mer(0,k);//get the i-th mer
        uint64_t prev_kmer=kmer;
        uint64_t min_seq=kmer&((0xffffffffffffffff<<(64-2*m))>>(64-2*m));//
        uint64_t canon_min_seq;
        int pos_min_seq=k-m;
        std::tuple<uint64_t,bool> seq_data=reverseComplementCanonical(kmer,k);//compute canonical k-mer and its orientation in graph
        uint64_t position=(uint64_t)id;//assign unitig id to position
        position=position<<32;//assign 0 as position in id and the computed orientation
        uint64_t minimizer=compute_minimizer_position(std::get<0>(seq_data),pos_min);//minimizer of kmer
        uint64_t old_minimizer=minimizer;//old minimizer for comparison only
        int counter=1;
        if(!std::get<1>(seq_data)){
            pos_min=(k-m-pos_min)%(k-m+1);
        }
        //when we have one kmer in the seq, write the data to the super-bucket
        if(node.second.unitig_length()==k){
            *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(old_minimizer)<<"\n";//write minimizer to super bucket file
            *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(prev_kmer)<<"\n";//write k-mer position to output file
            *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(position)<<"\n";//write k-mer position to output file
            *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(counter)<<"\n";//write k-mer position to output file
            all_mphfs[old_minimizer].mphf_size+=counter;
        }
        for(int i=1;i<node.second.unitig_length()-k+1;i++){
            pos_min_seq++;
            kmer=node.second.get_next_mer(kmer,i,k);//compute next k-mer
            seq_data=reverseComplementCanonical(kmer,k);//compute canonical form of k-mer and its orientation
            min_seq=node.second.get_next_mer(min_seq,pos_min_seq,m);
            canon_min_seq=canonical_bits(min_seq,m);
            if(canon_min_seq<minimizer){
                minimizer=canon_min_seq;
                pos_min=i+k-m;
            }
            else{
                if(i>pos_min){
                    minimizer=compute_minimizer_position(std::get<0>(seq_data),pos_min);
                    if(!std::get<1>(seq_data)){
                        pos_min=(k-m-pos_min)%(k-m+1);
                    }
                    pos_min+=i;
                }
            }
            if(minimizer==old_minimizer){
                counter++;
                if(i==node.second.unitig_length()-k){//when same minimizer and last kmer write everything to the superbucket
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(old_minimizer)<<"\n";//write minimizer to super bucket file
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(prev_kmer)<<"\n";//write k-mer position to output file
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(position)<<"\n";//write k-mer position to output file
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(counter)<<"\n";//write k-mer position to output file
                    all_mphfs[old_minimizer].mphf_size+=counter;
                }
            }
            else{
                *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(old_minimizer)<<"\n";//write minimizer to super bucket file
                *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(prev_kmer)<<"\n";//write k-mer position to output file
                *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(position)<<"\n";//write k-mer position to output file
                *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(counter)<<"\n";//write k-mer position to output file
                all_mphfs[old_minimizer].mphf_size+=counter;
                counter=1;
                prev_kmer=kmer;
                position=(uint64_t)id;//assign unitig id to position
                position=(position<<32)|i;//assign i as position in unitig and the computed orientation
                old_minimizer=minimizer;
                if(i==node.second.unitig_length()-k){//when same minimizer and last kmer write everything to the superbucket
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(minimizer)<<"\n";//write minimizer to super bucket file
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(kmer)<<"\n";//write k-mer position to output file
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(position)<<"\n";//write k-mer position to output file
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(counter)<<"\n";//write k-mer position to output file
                    all_mphfs[old_minimizer].mphf_size+=counter;
                }
            }
       }
        graph.insert_unitig(id,node.second);
    }
    graph.set_max_node_id(id);
	for (uint64_t i(0); i < number_of_super_buckets; ++i) {
		*minimizer_outs[i]<<std::flush;
		delete (minimizer_outs[i]);
	}
}
void Index_mphf::update_unitig(Unitig seq,int id,int previous_id,int starting_position,int ending_position,bool keep_orient){
        uint64_t kmer=seq.get_ith_mer(starting_position,k);//get the i-th mer
        std::tuple<uint64_t,bool> seq_data=reverseComplementCanonical(kmer,k);//compute canonical k-mer and its orientation in graph
        uint64_t position=(uint64_t)id;//assign unitig id to position
        position=(position<<32)|std::get<1>(seq_data);//assign 0 as position in id and the computed orientation
        uint64_t minimizer=compute_canonical_minimizer(std::get<0>(seq_data),k,m);//compute minimizer
        uint64_t index_by_mphf;
        //check first if the mphf is small (i.e grouped with other)
        if(all_mphfs[minimizer].small){
            index_by_mphf=mphf_grouped_buckets(std::get<0>(seq_data));
            position_kmers_small_buckets[index_by_mphf]=position;
        }
        else{
            index_by_mphf=all_mphfs[minimizer].kmer_MPHF(std::get<0>(seq_data));
            position_kmers[minimizer][index_by_mphf]=position;
        }
        uint64_t j=1;
       for(int i=starting_position+1;i<=ending_position;i++){
            kmer=seq.get_next_mer(kmer,i,k);//compute next k-mer
            seq_data=reverseComplementCanonical(kmer,k);//compute canonical form of k-mer and its orientation
            position=(uint64_t)id;//assign unitig id to position
            position=(position<<32)|(j<<1)|std::get<1>(seq_data);//assign i as position in unitig and the computed orientation
            //position=(position<<32)|(j<<1)|std::get<1>(seq_data);//assign i as position in unitig and the computed orientation
            minimizer=compute_canonical_minimizer(std::get<0>(seq_data),k,m);//compute minimizer
            //check first if the mphf is small (i.e grouped with other)
            if(all_mphfs[minimizer].small){
                index_by_mphf=mphf_grouped_buckets(std::get<0>(seq_data));
                position_kmers_small_buckets[index_by_mphf]=position;
            }
            else{
                index_by_mphf=all_mphfs[minimizer].kmer_MPHF(std::get<0>(seq_data));
                position_kmers[minimizer][index_by_mphf]=position;
            }
            j++;
       }
}
int Index_mphf::get_k_length(){
    return k;
}
void Index_mphf::rearrange_positions(pthash::single_phf<pthash::murmurhash2_64,pthash::dictionary_dictionary,true> mphf_ref,std::vector<uint64_t> kmers,std::vector<uint64_t> positions,uint64_t min_i,bool grouped){
    std::vector<uint64_t> kmer_pos(kmers.size());
    for(size_t i=0;i<kmers.size();i++){
        kmer_pos[mphf_ref(kmers[i])]=positions[i];
    }
    if(!grouped){
        position_kmers[min_i]=kmer_pos;
    }
}
uint64_t Index_mphf::compute_minimizer_position(uint64_t kmer,uint64_t& position){
    uint64_t res=0xffffffffffffffff;
    position=0;
    for(int i=0;i<k-m+1;i++){
        uint64_t mini=(kmer&((0xffffffffffffffff<<(64-2*k+2*i))>>(64-2*k+2*i)))>>(2*k-2*m-2*i);
        mini=canonical_bits(mini,m);
        if(mini<=res){
            res=mini;
            position=i;
        }
    }
    return res;
}
