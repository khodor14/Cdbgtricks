#include "index.hpp"
#include <unordered_set>
Index_mphf::Index_mphf(size_t k_size,size_t m_size,size_t log_super_buckets,size_t small_b_size,size_t multiplier_bucket){
    k=k_size;//k-mer size
    m=m_size;//minimizer size
    multiplier_bucketing=multiplier_bucket;
    log2_super_bucket=log_super_buckets;
    number_of_super_buckets=1ULL<<(2*log2_super_bucket);
    bucket_per_super_bucket=1ULL<<(2*(m-log2_super_bucket));
    minimizer_number=1ULL<<(2*m);
    super_bucket_max_size=multiplier_bucketing*small_b_size;
    mphfs_info.resize(minimizer_number);
    small_bucket_size=small_b_size;
}
void Index_mphf::build(GfaGraph& graph){
    for(uint64_t i=0;i<minimizer_number;i++){
        mphfs_info[i].mphf_size=0;
        mphfs_info[i].empty=true;
    }
    std::vector<uint64_t> kmers;
    std::vector<uint64_t> positions;
    std::vector<uint64_t> track_minimizers;
    //create files for each super bucket
    prepare_super_buckets(graph);
    //read and create mphf if the bucket is not small
    //along the way group ungrouped mphf
    read_super_buckets(graph,kmers,positions,track_minimizers);
    //the last group
    if(kmers.size()>0){
        create_mphf_per_super_bucket(kmers,positions,track_minimizers,all_mphfs.size());
    }
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
        uint64_t position=(uint64_t)node.first;//assign unitig id to position
        position=position<<32;//assign 0 as position in id and the computed orientation
        uint64_t minimizer=compute_minimizer_position(kmer,pos_min);//minimizer of kmer
        uint64_t old_minimizer=minimizer;//old minimizer for comparison only
        uint64_t hash_min=unrevhash_min(minimizer);
        int counter=1;
        //when we have one kmer in the seq, write the data to the super-bucket
        if(node.second.unitig_length()==k){
            old_minimizer=revhash_min(old_minimizer)%minimizer_number;
            *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(old_minimizer)<<"\n";//write minimizer to super bucket file
            *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(prev_kmer)<<"\n";//write k-mer position to output file
            *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(position)<<"\n";//write k-mer position to output file
            *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(counter)<<"\n";//write k-mer position to output file
            mphfs_info[old_minimizer].mphf_size+=counter;
            mphfs_info[old_minimizer].empty=false;   
        }
        for(int i=1;i<node.second.unitig_length()-k+1;i++){
            pos_min_seq++;
            kmer=node.second.get_next_mer(kmer,i,k);//compute next k-mer
            min_seq=node.second.get_next_mer(min_seq,pos_min_seq,m);
            canon_min_seq=canonical_bits(min_seq,m);
            uint64_t new_hash=unrevhash_min(canon_min_seq);
            if(new_hash<hash_min){//new minimizer is found
                minimizer=canon_min_seq;
                hash_min=new_hash;
                pos_min=i+k-m;
            }
            else{
                if(i>pos_min){//the minimizer is outdated
                    minimizer=compute_minimizer_position(kmer,pos_min);
                    hash_min=unrevhash_min(minimizer);
                    pos_min=pos_min+i;
                }
            }
            if(revhash_min(minimizer)%minimizer_number==revhash_min(old_minimizer)%minimizer_number){
                counter++;
                if(i==node.second.unitig_length()-k){//when same minimizer and last kmer write everything to the superbucket
                    old_minimizer=revhash_min(old_minimizer)%minimizer_number;
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(old_minimizer)<<"\n";//write minimizer to super bucket file
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(prev_kmer)<<"\n";//write k-mer position to output file
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(position)<<"\n";//write k-mer position to output file
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(counter)<<"\n";//write k-mer position to output file
                    mphfs_info[old_minimizer].mphf_size+=counter;
                    mphfs_info[old_minimizer].empty=false;
                }
            }
            else{
                old_minimizer=revhash_min(old_minimizer)%minimizer_number;
                *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(old_minimizer)<<"\n";//write minimizer to super bucket file
                *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(prev_kmer)<<"\n";//write k-mer position to output file
                *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(position)<<"\n";//write k-mer position to output file
                *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(counter)<<"\n";//write k-mer position to output file
                mphfs_info[old_minimizer].mphf_size+=counter;
                mphfs_info[old_minimizer].empty=false;
                counter=1;
                prev_kmer=kmer;
                position=(uint64_t)node.first;//assign unitig id to position
                position=(position<<32)|i;//assign i as position in unitig and the computed orientation
                old_minimizer=minimizer;
                if(i==node.second.unitig_length()-k){//when same minimizer and last kmer write everything to the superbucket
                    minimizer=revhash_min(minimizer)%minimizer_number;
                    *(minimizer_outs[minimizer/bucket_per_super_bucket])<<std::to_string(minimizer)<<"\n";//write minimizer to super bucket file
                    *(minimizer_outs[minimizer/bucket_per_super_bucket])<<std::to_string(kmer)<<"\n";//write k-mer position to output file
                    *(minimizer_outs[minimizer/bucket_per_super_bucket])<<std::to_string(position)<<"\n";//write k-mer position to output file
                    *(minimizer_outs[minimizer/bucket_per_super_bucket])<<std::to_string(counter)<<"\n";//write k-mer position to output file
                    mphfs_info[minimizer].mphf_size+=counter;
                    mphfs_info[minimizer].empty=false;
                }
            }
       }
    }
	for (uint64_t i(0); i < number_of_super_buckets; ++i) {
		*minimizer_outs[i]<<std::flush;
		delete (minimizer_outs[i]);
	}
    uint64_t total_kmer_in_super_buckets=0;
    for (uint64_t i(0); i < minimizer_number; ++i) {
        if(!mphfs_info[i].empty && mphfs_info[i].mphf_size<small_bucket_size){
            total_kmer_in_super_buckets=total_kmer_in_super_buckets+mphfs_info[i].mphf_size;
        }
	}
    num_sup_buckets=total_kmer_in_super_buckets/super_bucket_max_size;
}
void Index_mphf::read_super_buckets(GfaGraph& graph,std::vector<uint64_t> &kmers,std::vector<uint64_t>& positions,
                                    std::vector<uint64_t>& track_minimizers_in_super_bucket){
    /*
    takes the graph, an empty vector of kmers and empty vector of their position
    a boolean update phase: this function is used for construction and update
    */
    uint64_t super_bucket_created=0;
    for(uint64_t super_bucket=0;super_bucket<number_of_super_buckets;super_bucket++){
        std::string filename("ccdbgupdater_out" + std::to_string(super_bucket) + ".gz");
        std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> super_bucket_data;
        read_super_file(filename,super_bucket_data);
        create_mphf(kmers,positions,graph,super_bucket_data,track_minimizers_in_super_bucket,super_bucket_created);
    }
}
uint64_t Index_mphf::kmer_position(uint64_t kmer){
    uint64_t pos=0;
    uint64_t minimizer=compute_minimizer_position(kmer,pos);//compute the minimizer of the k-mer
    minimizer=revhash_min(minimizer)%minimizer_number;
    if(mphfs_info[minimizer].empty){
        //0 is undefined k-mer position because the position is a triplet (unitig id,pos in unitig,orientation) represented in 64 bits
        // as id is not zero then the defined position cannot be zero
        return 0;//the mphf related to the minimizer of the k-mer is not yet created so the k-mer does not exist in the graph
    }
    uint64_t bucket_id=mphfs_info[minimizer].bucket_id;
    uint64_t h_val=all_mphfs[bucket_id].lookup(kmer)%position_kmers[bucket_id].size();//compute hash value and use module to capture out of range returned by bbhash
    return position_kmers[bucket_id][h_val];//postion of this k-mer in the graph
}
void Index_mphf::read_super_file(std::string filename,std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> &super_bucket_data){
    zstr::ifstream in(filename);
    in.peek();
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
    std::remove(filename.c_str());
}
void Index_mphf::create_mphf(std::vector<uint64_t>& kmers,std::vector<uint64_t>& positions,GfaGraph& graph,
                                std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> superkeys,
                                std::vector<uint64_t>& track_minimizers_in_super_bucket,uint64_t &super_bucket_created){
    /*
        takes a vector of kmers to which the kmers from small buckets will be added
              a vector of positions to which the positions of kmers from small buckets will be added
              A graph to retrieve the kmers associated with the super keys
              a map whose keys are the minimizers and values is a vector of tuples (kmer, position,counter) where kmer is the first kmer of the superkmer
                                                                    position is the position of this kmer, and counter is the length of super kmer (number of kmer)
        if the bucket is small => create the mphf
            else: group the kmers from small buckets
    */
    for(auto minimizer_data:superkeys){
        if(mphfs_info[minimizer_data.first].mphf_size<small_bucket_size){
            mphfs_info[minimizer_data.first].small=true;
            if(kmers.size()>=super_bucket_max_size && super_bucket_created<num_sup_buckets-1){
                super_bucket_created+=1;
                uint64_t bucket_id=all_mphfs.size();
                create_mphf_per_super_bucket(kmers,positions,track_minimizers_in_super_bucket,bucket_id);
                //finding smallest super-bucket id
                if(kmers.size()<smallest_super_bucket_size){
                    smallest_super_bucket_size=kmers.size();
                    smallest_bucket_id=bucket_id;
                }
                kmers.clear();//empty the buckets
                positions.clear();//empty the buckets
                track_minimizers_in_super_bucket.clear();//empty the buckets
            }
            track_minimizers_in_super_bucket.push_back(minimizer_data.first);
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
            mphfs_info[minimizer_data.first].bucket_id=all_mphfs.size();
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
            MPHF kmer_MPHF=boomphf::mphf<uint64_t,hasher_t>(minimizer_kmers.size(),minimizer_kmers);
            uint64_t bucket_id=all_mphfs.size();
            all_mphfs.push_back(kmer_MPHF);
            mphfs_info[minimizer_data.first].bucket_id=bucket_id;
            //finding smallest bucket id
            if(kmer_positions.size()<smallest_bucket_size){
                smallest_bucket_size=kmer_positions.size();
                smallest_bucket_id=bucket_id;
            }
            position_kmers.push_back(std::vector<uint64_t>());//push back an empty vector
            rearrange_positions(kmer_MPHF,minimizer_kmers,kmer_positions,bucket_id);
        }
    }
}
void Index_mphf::read_super_bucket_update(GfaGraph& graph, uint64_t &num_new_kmers_new_supb,uint64_t &last_super_bucket_id,
    std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> &kmers_new_super,
    std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> &kmers_super_b_updates){
    for(uint64_t super_bucket=0;super_bucket<number_of_super_buckets;super_bucket++){
        std::string filename("ccdbgupdater_out" + std::to_string(super_bucket) + ".gz");
        std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> super_bucket_data;
        read_super_file(filename,super_bucket_data);
        update_mphfs(graph,num_new_kmers_new_supb,last_super_bucket_id,kmers_new_super,kmers_super_b_updates,super_bucket_data);
    }

}
void Index_mphf::update_mphfs(GfaGraph& graph, uint64_t &num_new_kmers_new_supb,uint64_t &last_super_bucket_id,
                                std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> &kmers_new_super,
                                std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> &kmers_super_b_updates,
                                std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> superkeys){
    for(auto minimizer_data:superkeys){
        if(mphfs_info[minimizer_data.first].small){
            uint64_t super_bucket_id=mphfs_info[minimizer_data.first].bucket_id;
            last_super_bucket_id=super_bucket_id;
            if(kmers_super_b_updates.count(super_bucket_id)==0){
                kmers_super_b_updates[super_bucket_id]=std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>();
            }
            for(auto super_key:minimizer_data.second){
                kmers_super_b_updates[super_bucket_id].push_back(super_key);//insert super-keys to the vector of super-bucket
            }
        }
        else if(mphfs_info[minimizer_data.first].mphf_size<small_bucket_size){
            num_new_kmers_new_supb+=mphfs_info[minimizer_data.first].mphf_size;
            mphfs_info[minimizer_data.first].empty=false;
            mphfs_info[minimizer_data.first].small=true;
            kmers_new_super[minimizer_data.first]=minimizer_data.second;
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
            uint64_t bucket_id=0xffffffffffffffff;//undefined bucket id
            if(!mphfs_info[minimizer_data.first].empty){
                bucket_id=mphfs_info[minimizer_data.first].bucket_id;
                for(uint64_t pos_kmer:position_kmers[bucket_id]){
                    kmer_positions.push_back(pos_kmer);
                    uint64_t km=graph.get_kmer(pos_kmer,k);
                    minimizer_kmers.push_back(canonical_bits(km,k));
                }
            }
            if(bucket_id==0xffffffffffffffff){//the minimizer is not seen yet and need to be added 
                bucket_id=mphfs_info.size();
                MPHF kmer_MPHF=boomphf::mphf<uint64_t,hasher_t>(minimizer_kmers.size(),minimizer_kmers);
                all_mphfs.push_back(kmer_MPHF);
                mphfs_info[minimizer_data.first].bucket_id=bucket_id;
                position_kmers.push_back(std::vector<uint64_t>());
            }
            else{
                //we already have this minimizer in the index, so we need to update the mphf of its bucket
                all_mphfs[mphfs_info[minimizer_data.first].bucket_id]=boomphf::mphf<uint64_t,hasher_t>(minimizer_kmers.size(),minimizer_kmers);
            }
            rearrange_positions(all_mphfs[bucket_id],minimizer_kmers,kmer_positions,bucket_id);
        }
    }
}
template <typename T>
void Index_mphf::create_mphf_per_super_bucket(std::vector<uint64_t>& kmers,std::vector<uint64_t>& positions,
        const T& track_minimizer,uint64_t bucket_id){
    //takes the k-mers, their positions and the minimizers of this super-bucket
    MPHF mphf_super_bucket=boomphf::mphf<uint64_t,hasher_t>(kmers.size(),kmers);
    std::vector<uint64_t> rearranged_positions(kmers.size());
    for(uint64_t i=0;i<kmers.size();i++){
        rearranged_positions[mphf_super_bucket.lookup(kmers[i])]=positions[i];
    }
    position_kmers.push_back(rearranged_positions);
    all_mphfs.push_back(mphf_super_bucket);
    for(uint64_t minimizer:track_minimizer){
        mphfs_info[minimizer].bucket_id=bucket_id;
    }
}
template <typename T>
void Index_mphf::update_super_bucket(std::vector<uint64_t>& kmers,std::vector<uint64_t>& positions,
        const T& track_minimizer,uint64_t bucket_id){
    //takes the k-mers, their positions and the minimizers of this super-bucket
    std::vector<uint64_t> rearranged_positions(kmers.size());
    all_mphfs[bucket_id]=boomphf::mphf<uint64_t,hasher_t>(kmers.size(),kmers);
    MPHF mphf_super_bucket=all_mphfs[bucket_id];
    for(uint64_t i=0;i<kmers.size();i++){
        rearranged_positions[mphf_super_bucket.lookup(kmers[i])]=positions[i];
    }
    position_kmers[bucket_id]=rearranged_positions;
    for(uint64_t minimizer:track_minimizer){
        mphfs_info[minimizer].bucket_id=bucket_id;
    }
}
void Index_mphf::update_index(std::unordered_map<int,Unitig>& constructed_funitigs,GfaGraph& graph){
    //create files for each super bucket and write kmers from funitigs
    //
    extract_kmers_from_funitigs(constructed_funitigs,graph);
    //read and create mphf if the bucket is not small
    //along the way group ungrouped mphf
    uint64_t num_new_kmers_new_supb=0;
    uint64_t last_super_bucket_id=0xffffffffffffffff;//id of last super bucket created before update to be used for adding new 
    std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> kmers_new_super;
    std::unordered_map<uint64_t,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> kmers_super_b_updates;
    read_super_bucket_update(graph,num_new_kmers_new_supb,last_super_bucket_id,kmers_new_super,kmers_super_b_updates);
    update_all_super_buckets(graph,kmers_super_b_updates,last_super_bucket_id);
    //create new super bucket for new k-mers with new minimizers
    bool merged_with_super_b=false;
    if(num_new_kmers_new_supb>=super_bucket_max_size){
        std::vector<uint64_t> track_minimizers;
        std::vector<uint64_t> kmers;
        std::vector<uint64_t> positions;
        uint64_t num_new_super_bucket=num_new_kmers_new_supb/super_bucket_max_size;
        uint64_t new_created_super_bucket=1;
        for(auto minimizerData:kmers_new_super){
            if(kmers.size()>super_bucket_max_size && new_created_super_bucket<num_new_super_bucket){
                new_created_super_bucket++;
                create_mphf_per_super_bucket(kmers,positions,track_minimizers,all_mphfs.size());
                kmers.clear();
                positions.clear();
                track_minimizers.clear();
            }
            track_minimizers.push_back(minimizerData.first);
            for(auto super_key_data:minimizerData.second){
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
        //create the new last super bucket
        create_mphf_per_super_bucket(kmers,positions,track_minimizers,all_mphfs.size());
    }
    //we have few kmers we add them to the last super bucket that should be updated
    else if(last_super_bucket_id!=0xffffffffffffffff){
        merged_with_super_b=true;
        //add new super-keys to last super bucket
        for(auto minimizerData:kmers_new_super){
            mphfs_info[minimizerData.first].bucket_id=last_super_bucket_id;
            for(auto superkey:minimizerData.second){
                kmers_super_b_updates[last_super_bucket_id].push_back(superkey);
            }
        }
        //update last super bucket
        update_super_bucket(graph,last_super_bucket_id,kmers_super_b_updates[last_super_bucket_id]);
    }
    //we have few k-mers but we don't have super bucket that needs update
    else if(num_new_kmers_new_supb>0){
        
        std::vector<std::tuple<uint64_t,uint64_t,uint64_t>> new_super_keys;
        //we add these k-mers to the smallest super bucket created
        if(smallest_super_bucket_id!=0xffffffffffffffff){
            for(auto minimizerData:kmers_new_super){
                mphfs_info[minimizerData.first].bucket_id=smallest_super_bucket_id;
                for(auto superkey:minimizerData.second){
                    new_super_keys.push_back(superkey);
                }
            }
            update_super_bucket(graph,smallest_super_bucket_id,new_super_keys);
        }
        //we don't have super bucket in the index, so we add these k-mers to the smallest bucket and we create the first super bucket
        else
        {   
            for(auto minimizerData:kmers_new_super){
                mphfs_info[minimizerData.first].bucket_id=smallest_bucket_id;
                for(auto superkey:minimizerData.second){
                    new_super_keys.push_back(superkey);
                }
            }
            update_super_bucket(graph,smallest_bucket_id,new_super_keys);
        }
    }
    if(!merged_with_super_b && last_super_bucket_id!=0xffffffffffffffff){
        update_super_bucket(graph,last_super_bucket_id,kmers_super_b_updates[last_super_bucket_id]);
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
        uint64_t position=(uint64_t)id;//assign unitig id to position
        position=position<<32;//assign 0 as position in id and the computed orientation
        uint64_t minimizer=compute_minimizer_position(kmer,pos_min);//minimizer of kmer
        uint64_t old_minimizer=minimizer;//old minimizer for comparison only
        uint64_t hash_min=unrevhash_min(minimizer);
        int counter=1;
        //when we have one kmer in the seq, write the data to the super-bucket
        if(node.second.unitig_length()==k){
            old_minimizer=revhash_min(old_minimizer)%minimizer_number;
            *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(old_minimizer)<<"\n";//write minimizer to super bucket file
            *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(prev_kmer)<<"\n";//write k-mer position to output file
            *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(position)<<"\n";//write k-mer position to output file
            *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(counter)<<"\n";//write k-mer position to output file
            mphfs_info[old_minimizer].mphf_size+=counter;
        }
        for(int i=1;i<node.second.unitig_length()-k+1;i++){
            pos_min_seq++;
            kmer=node.second.get_next_mer(kmer,i,k);//compute next k-mer
            min_seq=node.second.get_next_mer(min_seq,pos_min_seq,m);
            canon_min_seq=canonical_bits(min_seq,m);
            uint64_t new_hash=unrevhash_min(canon_min_seq);
            if(new_hash<hash_min){//new minimizer is found
                minimizer=canon_min_seq;
                hash_min=new_hash;
                pos_min=i+k-m;
            }
            else{
                if(i>pos_min){//the minimizer is outdated
                    minimizer=compute_minimizer_position(kmer,pos_min);
                    hash_min=unrevhash_min(minimizer);
                    pos_min=pos_min+i;
                }
            }
            if(revhash_min(minimizer)%minimizer_number==revhash_min(old_minimizer)%minimizer_number){
                counter++;
                if(i==node.second.unitig_length()-k){//when same minimizer and last kmer write everything to the superbucket
                    old_minimizer=revhash_min(old_minimizer)%minimizer_number;
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(old_minimizer)<<"\n";//write minimizer to super bucket file
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(prev_kmer)<<"\n";//write k-mer position to output file
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(position)<<"\n";//write k-mer position to output file
                    *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(counter)<<"\n";//write k-mer position to output file
                    mphfs_info[old_minimizer].mphf_size+=counter;
                }
            }
            else{
                old_minimizer=revhash_min(old_minimizer)%minimizer_number;
                *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(old_minimizer)<<"\n";//write minimizer to super bucket file
                *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(prev_kmer)<<"\n";//write k-mer position to output file
                *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(position)<<"\n";//write k-mer position to output file
                *(minimizer_outs[old_minimizer/bucket_per_super_bucket])<<std::to_string(counter)<<"\n";//write k-mer position to output file
                mphfs_info[old_minimizer].mphf_size+=counter;
                counter=1;
                prev_kmer=kmer;
                position=(uint64_t)id;//assign unitig id to position
                position=(position<<32)|i;//assign i as position in unitig and the computed orientation
                old_minimizer=minimizer;
                if(i==node.second.unitig_length()-k){//when same minimizer and last kmer write everything to the superbucket
                   minimizer=revhash_min(minimizer)%minimizer_number;
                    *(minimizer_outs[minimizer/bucket_per_super_bucket])<<std::to_string(minimizer)<<"\n";//write minimizer to super bucket file
                    *(minimizer_outs[minimizer/bucket_per_super_bucket])<<std::to_string(kmer)<<"\n";//write k-mer position to output file
                    *(minimizer_outs[minimizer/bucket_per_super_bucket])<<std::to_string(position)<<"\n";//write k-mer position to output file
                    *(minimizer_outs[minimizer/bucket_per_super_bucket])<<std::to_string(counter)<<"\n";//write k-mer position to output file
                    mphfs_info[minimizer].mphf_size+=counter;
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
        uint64_t minimizer=compute_minimizer(kmer);//compute minimizer
        minimizer=revhash_min(minimizer)%minimizer_number;
        uint64_t bucket_id=mphfs_info[minimizer].bucket_id;
        uint64_t index_by_mphf=all_mphfs[bucket_id].lookup(std::get<0>(seq_data));
        position_kmers[bucket_id][index_by_mphf]=position;
        uint64_t j=1;
       for(int i=starting_position+1;i<=ending_position;i++){
            kmer=seq.get_next_mer(kmer,i,k);//compute next k-mer
            seq_data=reverseComplementCanonical(kmer,k);//compute canonical form of k-mer and its orientation
            position=(uint64_t)id;//assign unitig id to position
            position=(position<<32)|(j<<1)|std::get<1>(seq_data);//assign i as position in unitig and the computed orientation
            //position=(position<<32)|(j<<1)|std::get<1>(seq_data);//assign i as position in unitig and the computed orientation
            minimizer=compute_minimizer(kmer);//compute minimizer
            minimizer=revhash_min(minimizer)%minimizer_number;
            bucket_id=mphfs_info[minimizer].bucket_id;
            index_by_mphf=all_mphfs[bucket_id].lookup(std::get<0>(seq_data));
            position_kmers[bucket_id][index_by_mphf]=position;
            j++;
       }
}
int Index_mphf::get_k_length(){
    return k;
}
void Index_mphf::rearrange_positions(MPHF mphf_ref,std::vector<uint64_t> kmers,std::vector<uint64_t> positions,uint64_t bucket_id){
    std::vector<uint64_t> kmer_pos(kmers.size());
    for(size_t i=0;i<kmers.size();i++){
        kmer_pos[mphf_ref.lookup(kmers[i])]=positions[i];
    }
    position_kmers[bucket_id]=kmer_pos;
}
uint64_t Index_mphf::compute_minimizer_position(uint64_t kmer,uint64_t& position){
    uint64_t mini, mmer;
	mmer = kmer%minimizer_number;//m-suffix of the k-mer
	mini = mmer = canonical_bits(mmer, m);
	uint64_t hash_mini=(unrevhash_min(mmer));//hash the minimizer
	position = k-m;
	for (uint64_t i(1); i<=k-m; i++) {
		kmer>>=2;
		mmer=kmer%minimizer_number;
		mmer=canonical_bits(mmer, m);
		uint64_t hash=(unrevhash_min(mmer));
		if (hash_mini>hash) {
			position=k-m-i;
			mini=mmer;
			hash_mini=hash;
		}
	}
	return mini;
}
uint64_t Index_mphf::compute_minimizer(uint64_t kmer){
    uint64_t mini, mmer;
	mmer = kmer%minimizer_number;//m-suffix of the k-mer
	mini = mmer = canonical_bits(mmer, m);
	uint64_t hash_mini=(unrevhash_min(mmer));//hash the minimizer
	for (uint64_t i(1); i<=k-m; i++) {
		kmer>>=2;
		mmer=kmer%minimizer_number;
		mmer=canonical_bits(mmer, m);
		uint64_t hash=(unrevhash_min(mmer));
		if (hash_mini>hash) {
			mini=mmer;
			hash_mini=hash;
		}
	}
	return mini;    
}
void Index_mphf::update_all_super_buckets(GfaGraph& graph,std::unordered_map<uint64_t,
                                std::vector<std::tuple<uint64_t,uint64_t,uint64_t>>> kmers_super_b_updates,
                                uint64_t last_super_bucket_id){
    for(auto elem:kmers_super_b_updates){
        if(elem.first!=last_super_bucket_id){//skip last super-bucket
            update_super_bucket(graph,elem.first,elem.second);
        }
    }
}
void Index_mphf::update_super_bucket(GfaGraph& graph,uint64_t super_bucket_id,std::vector<std::tuple<uint64_t,uint64_t,uint64_t>> &new_super_keys){
    //take super bucket id with new set of k-mers
    //update or split the super bucket
    uint64_t super_bucket_new_size=position_kmers[super_bucket_id].size();
    for(auto super_key:new_super_keys){
        super_bucket_new_size+=std::get<2>(super_key);
    }
    //we recreate one mphf
    if(super_bucket_new_size<2*super_bucket_max_size){
        std::vector<uint64_t> kmers;
        std::vector<uint64_t> positions;
        //parse new kmers
        for(auto super_key:new_super_keys){
            uint64_t k_mer=std::get<0>(super_key);
            uint64_t pos_kmer=std::get<1>(super_key);
            uint64_t counter=std::get<2>(super_key);
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
        for(auto position:position_kmers[super_bucket_id]){
            positions.push_back(position);
            kmers.push_back(canonical_bits(graph.get_kmer(position,k),k));
        }
        all_mphfs[super_bucket_id]=boomphf::mphf<uint64_t,hasher_t>(kmers.size(),kmers);
        rearrange_positions(all_mphfs[super_bucket_id],kmers,positions,super_bucket_id);
    }
    //we do greedy split into two super-buckets
    else{
        uint64_t counter_b_1=0;//number of k-mer added to first bucket
        uint64_t counter_b_2=0;//number of k-mers added to second bucket
        std::unordered_set<uint64_t> minimizer_b_1;//track minimizer first super bucket
        std::unordered_set<uint64_t> minimizer_b_2;//track minimizer second super bucket
        std::vector<uint64_t> kmer_b_1;//k-mers in first bucket
        std::vector<uint64_t> kmer_b_2;//k-mers in second bucket
        std::vector<uint64_t> positions_b_1;//positions in first bucket
        std::vector<uint64_t> positions_b_2;//positions in second bucket
        /*
         first check the minimizer is associated to which bucket and the insert the k-mer into it
         if no bucket is found, then check the number of k-mer inserted so far to the bucket and insert the k-mer to the smallest bucket
         
        */
        for(auto pos_kmer:position_kmers[super_bucket_id]){
            uint64_t kmer=canonical_bits(graph.get_kmer(pos_kmer,k),k);
            uint64_t minimizer=compute_minimizer(kmer);
            minimizer=revhash_min(minimizer)%minimizer_number;
            if(minimizer_b_1.count(minimizer)>0){
                kmer_b_1.push_back(kmer);
                positions_b_1.push_back(pos_kmer);
                counter_b_1++;
            }
            else if (minimizer_b_2.count(minimizer)>0)
            {
                kmer_b_2.push_back(kmer);
                positions_b_2.push_back(pos_kmer);
                counter_b_2++;
            }
            else if(counter_b_1<=counter_b_2){
                minimizer_b_1.emplace(minimizer);
                kmer_b_1.push_back(kmer);
                positions_b_1.push_back(pos_kmer);
                counter_b_1++;
            }
            else{
                minimizer_b_2.emplace(minimizer);
                kmer_b_2.push_back(kmer);
                positions_b_2.push_back(pos_kmer);
                counter_b_2++;
            }

        }
        for(auto super_key:new_super_keys){
            uint64_t k_mer=std::get<0>(super_key);
            uint64_t minimizer=compute_minimizer(k_mer);
            minimizer=revhash_min(minimizer)%minimizer_number;
            uint64_t pos_kmer=std::get<1>(super_key);
            uint64_t counter=std::get<2>(super_key);
            std::tuple<uint64_t,bool> seq_data=reverseComplementCanonical(k_mer,k);
            uint64_t id=pos_kmer>>32;
            uint64_t pos_in_unitig=pos_kmer&0xffffffff;
            uint64_t full_position=(id<<32)|(pos_in_unitig<<1)|std::get<1>(seq_data);
            Unitig seq=graph.get_unitig(id);
            if(minimizer_b_1.count(minimizer)>0){
                kmer_b_1.push_back(std::get<0>(seq_data));
                positions_b_1.push_back(full_position);
                counter_b_1++;
            }
            else if (minimizer_b_2.count(minimizer)>0)
            {
                kmer_b_2.push_back(std::get<0>(seq_data));
                positions_b_2.push_back(full_position);
                counter_b_2++;
            }
            else if(counter_b_1<=counter_b_2){
                minimizer_b_1.emplace(minimizer);
                kmer_b_1.push_back(std::get<0>(seq_data));
                positions_b_1.push_back(full_position);
                counter_b_1++;
            }
            else{
                minimizer_b_2.emplace(minimizer);
                kmer_b_2.push_back(std::get<0>(seq_data));
                positions_b_2.push_back(full_position);
                counter_b_2++;
            }
            for(uint64_t i=1;i<counter;i++){
                k_mer=seq.get_next_mer(k_mer,pos_in_unitig+i,k);
                seq_data=reverseComplementCanonical(k_mer,k);
                full_position=(id<<32)|((pos_in_unitig+i)<<1)|std::get<1>(seq_data);
                if(minimizer_b_1.count(minimizer)>0){
                    kmer_b_1.push_back(std::get<0>(seq_data));
                    positions_b_1.push_back(full_position);
                    counter_b_1++;
                }
                else if (minimizer_b_2.count(minimizer)>0)
                {
                    kmer_b_2.push_back(std::get<0>(seq_data));
                    positions_b_2.push_back(full_position);
                    counter_b_2++;
                }
                else if(counter_b_1<=counter_b_2){
                    minimizer_b_1.emplace(minimizer);
                    kmer_b_1.push_back(std::get<0>(seq_data));
                    positions_b_1.push_back(full_position);
                    counter_b_1++;
                }
                else{
                    minimizer_b_2.emplace(minimizer);
                    kmer_b_2.push_back(std::get<0>(seq_data));
                    positions_b_2.push_back(full_position);
                    counter_b_2++;
                }
            }
        }
        bool singleton_b1=false;// we may get only one minimizer in b1 so we need to update the info in the index
        bool singleton_b2=false;// we may get only one minimizer in b2 so we need to update the info in the index
        //to create a mphf for b1 and for b2, their sizes must respect the threshold
        if(kmer_b_1.size()>small_bucket_size && kmer_b_2.size()>small_bucket_size){
            update_super_bucket(kmer_b_1,positions_b_1,minimizer_b_1,super_bucket_id);//create mphf for first super-bucket
            create_mphf_per_super_bucket(kmer_b_2,positions_b_2,minimizer_b_2,all_mphfs.size());//create mphf for second super-bucket
            if(minimizer_b_1.size()==1){
                singleton_b1=true;
            }
            if(minimizer_b_2.size()==1){
                singleton_b2=true;
            }
        }
        // In case b1 does not respect the threshold strategy, we have to merge them again
        else if(kmer_b_1.size()<kmer_b_2.size()){
            for(int i=0;i<kmer_b_1.size();i++){
                kmer_b_2.push_back(kmer_b_1[i]);
                positions_b_2.push_back(positions_b_1[i]);
            }
            //
            if(kmer_b_1.size()==0 && minimizer_b_2.size()==1){
                singleton_b2=true;
            }
            update_super_bucket(kmer_b_2,positions_b_2,minimizer_b_2,super_bucket_id);//create mphf for first super-bucket
        }
        else{
            for(int i=0;i<kmer_b_2.size();i++){
                kmer_b_1.push_back(kmer_b_2[i]);
                positions_b_1.push_back(positions_b_2[i]);
            }
            if(kmer_b_2.size()==0 && minimizer_b_1.size()==1){
                singleton_b1=true;
            }
            update_super_bucket(kmer_b_1,positions_b_1,minimizer_b_1,super_bucket_id);//create mphf for second super-bucket
        }
        if(singleton_b1){
            uint64_t single_minimizer_b1=*(minimizer_b_1.begin());
            mphfs_info[single_minimizer_b1].small=false;
        }
        if(singleton_b2){
            uint64_t single_minimizer_b2=*(minimizer_b_2.begin());
            mphfs_info[single_minimizer_b2].small=false;
        }
    }
}
void Index_mphf::save(std::string filename){
    std::filebuf fb;
	std::remove(filename.c_str());
	fb.open(filename, std::ios::out | std::ios::binary | std::ios::trunc);
	zstr::ostream out(&fb);
    out.write(reinterpret_cast<const char*>(&k), sizeof(k));
    out.write(reinterpret_cast<const char*>(&m), sizeof(m));
    out.write(reinterpret_cast<const char*>(&log2_super_bucket), sizeof(log2_super_bucket));
    out.write(reinterpret_cast<const char*>(&small_bucket_size), sizeof(small_bucket_size));
    out.write(reinterpret_cast<const char*>(&multiplier_bucketing), sizeof(multiplier_bucketing));
    out.write(reinterpret_cast<const char*>(&super_bucket_max_size), sizeof(super_bucket_max_size));
    out.write(reinterpret_cast<const char*>(&minimizer_number), sizeof(minimizer_number));
    out.write(reinterpret_cast<const char*>(&bucket_per_super_bucket), sizeof(bucket_per_super_bucket));
    out.write(reinterpret_cast<const char*>(&number_of_super_buckets), sizeof(number_of_super_buckets));
    out.write(reinterpret_cast<const char*>(&num_sup_buckets), sizeof(num_sup_buckets));
    out.write(reinterpret_cast<const char*>(&smallest_super_bucket_id), sizeof(smallest_super_bucket_id));
    out.write(reinterpret_cast<const char*>(&smallest_super_bucket_size), sizeof(smallest_super_bucket_size));
    out.write(reinterpret_cast<const char*>(&smallest_bucket_id), sizeof(smallest_bucket_id));
    out.write(reinterpret_cast<const char*>(&smallest_bucket_size), sizeof(smallest_bucket_size));
    for(uint64_t i=0;i<mphfs_info.size();i++){
        out.write(reinterpret_cast<const char*>(&mphfs_info[i].bucket_id), sizeof(mphfs_info[i].bucket_id));
        out.write(reinterpret_cast<const char*>(&mphfs_info[i].created), sizeof(mphfs_info[i].created));
        out.write(reinterpret_cast<const char*>(&mphfs_info[i].empty), sizeof(mphfs_info[i].empty));
        out.write(reinterpret_cast<const char*>(&mphfs_info[i].mphf_size), sizeof(mphfs_info[i].mphf_size));
        out.write(reinterpret_cast<const char*>(&mphfs_info[i].small), sizeof(mphfs_info[i].small));
    }
    uint64_t num_mphf=all_mphfs.size();
    out.write(reinterpret_cast<const char*>(&num_mphf), sizeof(num_mphf));
    for(uint64_t i=0;i<all_mphfs.size();i++){
        uint64_t num_positions=position_kmers[i].size();
        out.write(reinterpret_cast<const char*>(&num_positions), sizeof(num_positions));
        out.write(reinterpret_cast<char const*>(position_kmers[i].data()), (std::streamsize)(sizeof(position_kmers[i][0]) * num_positions));
        all_mphfs[i].save(out);
    }
    out <<std::flush;
	fb.close();

}
void Index_mphf::load(std::string filename){
	zstr::ifstream in(filename);
    in.read(reinterpret_cast<char*>(&k), sizeof(k));
   in.read(reinterpret_cast<char*>(&m), sizeof(m));
    in.read(reinterpret_cast<char*>(&log2_super_bucket), sizeof(log2_super_bucket));
    in.read(reinterpret_cast<char*>(&small_bucket_size), sizeof(small_bucket_size));
    in.read(reinterpret_cast<char*>(&multiplier_bucketing), sizeof(multiplier_bucketing));
    in.read(reinterpret_cast<char*>(&super_bucket_max_size), sizeof(super_bucket_max_size));
    in.read(reinterpret_cast<char*>(&minimizer_number), sizeof(minimizer_number));
    in.read(reinterpret_cast<char*>(&bucket_per_super_bucket), sizeof(bucket_per_super_bucket));
    in.read(reinterpret_cast<char*>(&number_of_super_buckets), sizeof(number_of_super_buckets));
    in.read(reinterpret_cast<char*>(&num_sup_buckets), sizeof(num_sup_buckets));
    in.read(reinterpret_cast<char*>(&smallest_super_bucket_id), sizeof(smallest_super_bucket_id));
    in.read(reinterpret_cast<char*>(&smallest_super_bucket_size), sizeof(smallest_super_bucket_size));
    in.read(reinterpret_cast<char*>(&smallest_bucket_id), sizeof(smallest_bucket_id));
    in.read(reinterpret_cast<char*>(&smallest_bucket_size), sizeof(smallest_bucket_size));
    mphfs_info.resize(minimizer_number);
    for(uint64_t i=0;i<mphfs_info.size();i++){
        in.read(reinterpret_cast<char*>(&mphfs_info[i].bucket_id), sizeof(mphfs_info[i].bucket_id));
        in.read(reinterpret_cast<char*>(&mphfs_info[i].created), sizeof(mphfs_info[i].created));
        in.read(reinterpret_cast<char*>(&mphfs_info[i].empty), sizeof(mphfs_info[i].empty));
        in.read(reinterpret_cast<char*>(&mphfs_info[i].mphf_size), sizeof(mphfs_info[i].mphf_size));
        in.read(reinterpret_cast<char*>(&mphfs_info[i].small), sizeof(mphfs_info[i].small));
    }
    uint64_t num_mphf;
    in.read(reinterpret_cast<char*>(&num_mphf), sizeof(num_mphf));
    all_mphfs.resize(num_mphf);
    position_kmers.resize(num_mphf);
    for(uint64_t i=0;i<all_mphfs.size();i++){
        uint64_t num_positions=position_kmers[i].size();
		uint64_t sizer;
		in.read(reinterpret_cast<char *>(&sizer),  sizeof(size_t));
		position_kmers[i].resize(sizer);
		in.read(reinterpret_cast<char*>(position_kmers[i].data()), (std::streamsize)(sizeof(position_kmers[i][0]) * sizer));
        all_mphfs[i].load(in);
    }
}
