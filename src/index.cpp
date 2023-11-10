#include "index.hpp"
#include <unordered_set>
Index_mphf::Index_mphf(size_t k_size,size_t m_size,size_t small_b_size){
    k=k_size;
    m=m_size;
    minimizer_number=1ULL<<(2*m);
    all_mphfs.resize(minimizer_number);
    position_kmers.resize(minimizer_number);
    small_bucket_size=small_b_size;
}
uint64_t Index_mphf::kmer_position(uint64_t kmer){
    uint64_t minimizer=compute_canonical_minimizer(kmer,k,m);//compute the minimizer of the k-mer
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
    std::vector<std::ostream*> minimizer_outs;//write k-mers per file
	for (uint64_t i(0); i < minimizer_number; ++i) {
        all_mphfs[i].empty=true;
        all_mphfs[i].mphf_size=0;
		auto out = new zstr::ofstream("ccdbgupdater_out" + std::to_string(i) + ".gz");
		if (not out->good()) {
			std::cout << "Problem with files opening" <<std::endl;
			exit(1);
		}
		minimizer_outs.push_back(out);
	}
    //go over unitigs
    for(auto node:graph.get_nodes()){
        uint64_t pos_min=0;
        uint64_t kmer=node.second.get_ith_mer(0,k);//get the i-th mer
        uint64_t min_seq=kmer&((0xffffffffffffffff<<(64-2*m))>>(64-2*m));//
        uint64_t canon_min_seq;
        int pos_min_seq=k-m;
        std::tuple<uint64_t,bool> seq_data=reverseComplementCanonical(kmer,k);//compute canonical k-mer and its orientation in graph
        uint64_t position=(uint64_t)node.first;//assign unitig id to position
        position=(position<<32)|std::get<1>(seq_data);//assign 0 as position in id and the computed orientation
        uint64_t minimizer=compute_minimizer_position(std::get<0>(seq_data),pos_min);
        *(minimizer_outs[minimizer])<<std::to_string(std::get<0>(seq_data))<<"\n";//write k-mer to temporary output file
        *(minimizer_outs[minimizer])<<std::to_string(position)<<"\n";//write k-mer position to output file
        all_mphfs[minimizer].mphf_size++;
        all_mphfs[minimizer].empty=false;
        if(!std::get<1>(seq_data)){
            pos_min=(k-m-pos_min)%(k-m);
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
                    pos_min=(k-m-pos_min)%(k-m);
                }
                    pos_min+=i;
                }
            }
            seq_data=reverseComplementCanonical(kmer,k);//compute canonical form of k-mer and its orientation
            position=(uint64_t)node.first;//assign unitig id to position
            position=(position<<32)|(i<<1)|std::get<1>(seq_data);//assign i as position in unitig and the computed orientation
            *(minimizer_outs[minimizer])<<std::to_string(std::get<0>(seq_data))<<"\n";//write k-mer to temporary output file
            *(minimizer_outs[minimizer])<<std::to_string(position)<<"\n";//write k-mer position to output file
            all_mphfs[minimizer].mphf_size++;
            all_mphfs[minimizer].empty=false;
       }
    }
	for (uint64_t i(0); i < minimizer_number; ++i) {
        if(all_mphfs[i].mphf_size<small_bucket_size){
            all_mphfs[i].small=true;
        }
		*minimizer_outs[i]<<std::flush;
		delete (minimizer_outs[i]);
	}
    create_mphfs();
}
void Index_mphf::create_mphfs(){
    //configuration to be used for all mphf
    pthash::build_configuration config;
    config.c = 6.0;
    config.alpha = 0.94;
    config.minimal_output = true;  // mphf
    config.verbose_output = false;
    config.seed=1234567890;
    std::vector<uint64_t> kmers_from_small_buckets;
    std::vector<uint64_t> positions_kmers;
    int grouped=0;
    int separated=0;
    for (uint64_t i(0); i < minimizer_number; ++i) {
        if(!all_mphfs[i].empty){
            //read the kmer file
            zstr::ifstream in("ccdbgupdater_out" + std::to_string(i) + ".gz");
            in.peek();
            if(all_mphfs[i].small){
                grouped++;
                while(not in.eof() and in.good()){
                    std::string line;
                    std::getline(in,line);
                    if(line.length()==0){
                        break;
                    }
                    uint64_t kmer=std::stoll(line);
                    std::getline(in,line);
                    uint64_t kmer_position=std::stoll(line);
                    kmers_from_small_buckets.push_back(kmer);
                    positions_kmers.push_back(kmer_position);             
                }               
            }
            else{
                separated++;
            uint64_t size=all_mphfs[i].mphf_size;
            std::vector<uint64_t> kmers(all_mphfs[i].mphf_size);
            std::vector<uint64_t> kmer_data(all_mphfs[i].mphf_size);
            size_t j=0;
            //check mphf size to decide weither to construct a mphf or to insert the k-mers to the hash table
            while(not in.eof() and in.good()){
                std::string line;
                std::getline(in,line);
                if(line.length()==0){
                    break;
                }
                uint64_t kmer=std::stoll(line);
                std::getline(in,line);
                
                uint64_t kmer_position=std::stoll(line);
                kmers[j]=kmer;
                kmer_data[j]=kmer_position;             
                j++;
            }
            all_mphfs[i].kmer_MPHF.build_in_internal_memory(kmers.begin(),kmers.size(),config);//build the mphf
           //rearrange the positions with respect to the mphf
           rearrange_positions(all_mphfs[i].kmer_MPHF,kmers,kmer_data,i,false);
            }
           
        }
        std::remove(("ccdbgupdater_out" + std::to_string(i) + ".gz").c_str());
	}
    mphf_grouped_buckets.build_in_internal_memory(kmers_from_small_buckets.begin(),kmers_from_small_buckets.size(),config);
    position_kmers_small_buckets.resize(kmers_from_small_buckets.size());
    for(uint64_t i(0);i<kmers_from_small_buckets.size();i++){
        position_kmers_small_buckets[mphf_grouped_buckets(kmers_from_small_buckets[i])]=positions_kmers[i];//
    }
}
void Index_mphf::update(GfaGraph& graph){
    //configuration to be used for all mphf
    pthash::build_configuration config;
    config.c = 6.0;
    config.alpha = 0.94;
    config.minimal_output = true;  // mphf
    config.verbose_output = true;
    std::vector<uint64_t> kmers_from_small_buckets;
    std::vector<uint64_t> positions_kmers;
    for (uint64_t i(0); i < minimizer_number; ++i) {
        std::vector<uint64_t> kmers;
        std::vector<uint64_t> kmer_data;
        if(!all_mphfs[i].empty && all_mphfs[i].created){
            std::unordered_map<uint64_t,bool> track_error;
            zstr::ifstream in("ccdbgupdater_out" + std::to_string(i) + ".gz");
            in.peek();
            if(all_mphfs[i].small){
                    while(not in.eof() and in.good()){
                        std::string line;
                        std::getline(in,line);
                        if(line.length()==0){
                            break;
                        }
                        uint64_t kmer=std::stoll(line);
                        std::getline(in,line);
                        uint64_t kmer_position=std::stoll(line);
                        kmers_from_small_buckets.push_back(kmer);
                        positions_kmers.push_back(kmer_position);             
                    }

            }
            else{

                int size=position_kmers[i].size();
                for(int j=0;j<size;j++){
                    uint64_t k_track=canonical_bits(graph.get_kmer(position_kmers[i][j],k),k);
                    kmers.push_back(k_track);
                    kmer_data.push_back(position_kmers[i][j]);
                }
                size_t j=size;
                while(not in.eof() and in.good()){
                    std::string line;
                    std::getline(in,line);
                    if(line.length()==0){
                        break;
                    }
                    uint64_t kmer=std::stoll(line);
                    std::getline(in,line);
                    uint64_t kmer_position=std::stoll(line);
                    kmers.push_back(kmer);
                    kmer_data.push_back(kmer_position);             
                }      
                all_mphfs[i].kmer_MPHF.build_in_internal_memory(kmers.begin(),kmers.size(),config);//build the mphf
                rearrange_positions(all_mphfs[i].kmer_MPHF,kmers,kmer_data,i,false);//arrange the position of k-mers in unitig based on the indices computed by the mphf
            }

        }
            //add one more boolean to say it is created
        else if(all_mphfs[i].created){
                //read the kmer file
                all_mphfs[i].empty=false;
                zstr::ifstream in("ccdbgupdater_out" + std::to_string(i) + ".gz");
                in.peek();
                if(all_mphfs[i].mphf_size<small_bucket_size){
                    //the mphf related to this bucket was not created and has a number of keys less than the threshold=> we group keys in one mphf
                    all_mphfs[i].small=true;
                    while(not in.eof() and in.good()){
                        std::string line;
                        std::getline(in,line);
                        if(line.length()==0){
                            break;
                        }
                        uint64_t kmer=std::stoll(line);
                        std::getline(in,line);
                        uint64_t kmer_position=std::stoll(line);
                        kmers_from_small_buckets.push_back(kmer);
                        positions_kmers.push_back(kmer_position);
                    }
                }
                else{
                    size_t j=0;
                    while(not in.eof() and in.good()){
                        std::string line;
                        std::getline(in,line);
                        if(line.length()==0){
                            break;
                        }
                        uint64_t kmer=std::stoll(line);
                        std::getline(in,line);
                        uint64_t kmer_position=std::stoll(line);
                        kmers[j]=kmer;
                        kmer_data[j]=kmer_position;             
                        j++;
                    }
                    all_mphfs[i].kmer_MPHF.build_in_internal_memory(kmers.begin(),kmers.size(),config);//build the mphf
                    rearrange_positions(all_mphfs[i].kmer_MPHF,kmers,kmer_data,i,false);
                }
            }
           std::remove(("ccdbgupdater_out" + std::to_string(i) + ".gz").c_str());
           all_mphfs[i].created=false;
        }
        if(kmers_from_small_buckets.size()>0){
            for(auto val:position_kmers_small_buckets){
                positions_kmers.push_back(val);
                kmers_from_small_buckets.push_back(graph.get_kmer(val,k));
            }
            position_kmers_small_buckets.resize(kmers_from_small_buckets.size());
            mphf_grouped_buckets.build_in_internal_memory(kmers_from_small_buckets.begin(),kmers_from_small_buckets.size(),config);
            for(uint64_t i(0);i<kmers_from_small_buckets.size();i++){
                position_kmers_small_buckets[mphf_grouped_buckets(kmers_from_small_buckets[i])]=positions_kmers[i];
            }
        }
}
void Index_mphf::extract_kmers_from_funitigs(std::unordered_map<int,Unitig>& constructed_unitigs,GfaGraph& graph){
    std::vector<std::ostream*> minimizer_outs(minimizer_number);//write k-mers per file
	for (uint64_t i(0); i < minimizer_number; ++i) {
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
        uint64_t min_seq=kmer&((0xffffffffffffffff<<(64-2*m))>>(64-2*m));//
        uint64_t canon_min_seq;
        int pos_min_seq=k-m;
        std::tuple<uint64_t,bool> seq_data=reverseComplementCanonical(kmer,k);//compute canonical k-mer and its orientation in graph
        uint64_t position=(uint64_t)id;//assign unitig id to position
        position=(position<<32)|std::get<1>(seq_data);//assign 0 as position in id and the computed orientation
        uint64_t minimizer=compute_minimizer_position(std::get<0>(seq_data),pos_min);
        *(minimizer_outs[minimizer])<<std::to_string(std::get<0>(seq_data))<<"\n";//write k-mer to temporary output file
        *(minimizer_outs[minimizer])<<std::to_string(position)<<"\n";//write k-mer position to output file
        all_mphfs[minimizer].mphf_size++;
        all_mphfs[minimizer].created=true;
        if(!std::get<1>(seq_data)){
            pos_min=(k-m-pos_min)%(k-m);
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
                    pos_min=(k-m-pos_min)%(k-m);
                }
                    pos_min+=i;
                }
            }
            seq_data=reverseComplementCanonical(kmer,k);//compute canonical form of k-mer and its orientation
            position=(uint64_t)id;//assign unitig id to position
            position=(position<<32)|(i<<1)|std::get<1>(seq_data);//assign i as position in unitig and the computed orientation
            *(minimizer_outs[minimizer])<<std::to_string(std::get<0>(seq_data))<<"\n";//write k-mer to temporary output file
            *(minimizer_outs[minimizer])<<std::to_string(position)<<"\n";//write k-mer position to output file
            all_mphfs[minimizer].mphf_size++;
            all_mphfs[minimizer].created=true;//to differentiate between mphf that should be created from those already created
       }
       graph.insert_unitig(id,node.second);
    }
    graph.set_max_node_id(id);
	for (uint64_t i(0); i < minimizer_number; ++i) {
		*minimizer_outs[i]<<std::flush;
		delete (minimizer_outs[i]);
	}
}
void Index_mphf::update_unitig(Unitig seq,int id,int previous_id,int starting_position,int ending_position,bool keep_orient){
        uint64_t kmer=seq.get_ith_mer(starting_position,k);//get the i-th mer
        uint64_t pos_min=0;
        uint64_t min_seq=kmer&((0xffffffffffffffff<<(64-2*m))>>(64-2*m));//
        uint64_t canon_min_seq;
        std::tuple<uint64_t,bool> seq_data=reverseComplementCanonical(kmer,k);//compute canonical k-mer and its orientation in graph
        uint64_t position=(uint64_t)id;//assign unitig id to position
        position=(position<<32)|std::get<1>(seq_data);//assign 0 as position in id and the computed orientation
        uint64_t minimizer=compute_minimizer_position(kmer,pos_min);//compute minimizer
        uint64_t index_by_mphf;
        int pos_min_seq=starting_position+k-m;
        //check first if the mphf is small (i.e grouped with other)
        /*if(all_mphfs[minimizer].small){
            index_by_mphf=mphf_grouped_buckets(std::get<0>(seq_data));
            position_kmers_small_buckets[index_by_mphf]=position;
        }
        else{
            index_by_mphf=all_mphfs[minimizer].kmer_MPHF(std::get<0>(seq_data));
            position_kmers[minimizer][index_by_mphf]=position;
        }*/
        std::cout<<pos_min<<" posmin\n";
        if(!std::get<1>(seq_data)){
            pos_min=(k-m-pos_min)%m;
        }
        std::cout<<to_string(kmer,k)<<"    "<<to_string(minimizer,m)<<"    "<<to_string(min_seq,m)<<"    0"<<"    "<<pos_min<<"\n";
        uint64_t j=1;
       for(int i=starting_position+1;i<=ending_position;i++){
            pos_min_seq++;
            kmer=seq.get_next_mer(kmer,i,k);//compute next k-mer
            seq_data=reverseComplementCanonical(kmer,k);//compute canonical form of k-mer and its orientation
            min_seq=seq.get_next_mer(min_seq,pos_min_seq,m);
            canon_min_seq=canonical_bits(min_seq,m);
            if(canon_min_seq<minimizer){
                std::cout<<"Found new minimizer\n";
                minimizer=canon_min_seq;
                pos_min=j+k-m;
            }
            else{
                if(j>pos_min){
                    std::cout<<"Recompute minimizer  "<<j<<"   "<<pos_min<<"\n";
                    minimizer=compute_minimizer_position(std::get<0>(seq_data),pos_min);
                    if(!std::get<1>(seq_data)){
                        pos_min=(k-m-pos_min)%m;
                    }
                    pos_min+=j;
                }
            }
            position=(uint64_t)id;//assign unitig id to position
            position=(position<<32)|(j<<1)|std::get<1>(seq_data);//assign i as position in unitig and the computed orientation
            std::cout<<to_string(kmer,k)<<"    "<<to_string(minimizer,m)<<"    "<<to_string(min_seq,m)<<"    "<<j<<"    "<<pos_min<<"\n";
            //check first if the mphf is small (i.e grouped with other)
            /*if(all_mphfs[minimizer].small){
                index_by_mphf=mphf_grouped_buckets(std::get<0>(seq_data));
                position_kmers_small_buckets[index_by_mphf]=position;
            }
            else{
                index_by_mphf=all_mphfs[minimizer].kmer_MPHF(std::get<0>(seq_data));
                position_kmers[minimizer][index_by_mphf]=position;
            }*/
            j++;
       }
}
int Index_mphf::get_k_length(){
    return k;
}
void Index_mphf::rearrange_positions(pthash::single_phf<pthash::murmurhash2_64,pthash::dictionary_dictionary,true> mphf_ref,std::vector<uint64_t>& kmers,std::vector<uint64_t>& positions,uint64_t min_i,bool grouped){
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
