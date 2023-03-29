#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/unordered_map.hpp>
#include <boost/serialization/utility.hpp>
#include "index_kmer.h"
#include "CommonUtils.h"
#include <unordered_map>
#include <algorithm>
#include <cstring>
#include <vector>
#include <string>
#include <tuple>
#include <iostream>
#include <fstream>
#include <cstdio>
#include "FileSerializer.hpp"
Index::Index(int buckets,int k_size){
    k=k_size;
    number_of_buckets=buckets;
    index_table.rehash(number_of_buckets);
}
int Index::get_k(){
    return k;
}
void Index::create(GfaGraph& graph){
   //index_table.set_empty_key(NULL);
   for (std::pair<int, Unitig> node : graph.get_nodes()){
        int id=node.first;
        int i;
        uint64_t i_th_mer=node.second.get_ith_mer(0,k-1);
        std::tuple<uint64_t,bool> seq_data=reverseComplementCanonical(i_th_mer,k-1);
        if(index_table.count(std::get<0>(seq_data))==0){
            std::vector<std::tuple<int,int,bool>> data;
            data.push_back(std::tuple<int,int,bool>(node.first,0,std::get<0>(seq_data)));
            index_table[std::get<0>(seq_data)]=data;
        }
        else{
            std::vector<std::tuple<int,int,bool>> data=index_table[std::get<0>(seq_data)];
            data.push_back(std::tuple<int,int,bool>(id,i,std::get<1>(seq_data)));
            index_table[std::get<0>(seq_data)]=data;
        }
        for(i=1;i<=node.second.unitig_length()-k+1;++i){
            i_th_mer=node.second.get_next_mer(i_th_mer,i,k-1);
            seq_data=reverseComplementCanonical(i_th_mer,k-1);
            if(index_table.count(std::get<0>(seq_data))==0){
                    std::vector<std::tuple<int,int,bool>> data_i;
                    data_i.push_back(std::tuple<int,int,bool>(id,i,std::get<1>(seq_data)));
                    index_table[std::get<0>(seq_data)]=data_i;
            }
            else{
                std::vector<std::tuple<int,int,bool>> data_i=index_table[std::get<0>(seq_data)];
                data_i.push_back(std::tuple<int,int,bool>(id,i,std::get<1>(seq_data)));
                index_table[std::get<0>(seq_data)]=data_i;
            }
        }
	}
}
std::vector<std::tuple<int,int,bool>> Index::find(uint64_t kmer){
    uint64_t canonical_kmer=canonical_bits(kmer,k-1);
    if(index_table.count(canonical_kmer)==0){
        return std::vector<std::tuple<int,int,bool>>();
    }
    return index_table[canonical_kmer];
}
void Index::update_k_1_mer(uint64_t k_1_mer,int prev_id,int current_id,int position,bool keep_orient){
    /*
        the input are: a string of size (k-1)
        the previous id where this string occurs
        the current id where it occurs after the updte performed
        the new position where it occurs
    */
    std::vector<std::tuple<int,int,bool>> data=find(k_1_mer);//get the occurences of the string from the index
    //std::tuple<uint64_t,bool> kmer_data=reverseComplementCanonical(k_1_mer,k-1);
    //go over the occurences
    for (int i=0;i<data.size();i++){
        //which occurens has the previous id, the old occurence
        if(std::get<0>(data[i])==prev_id){
            if(keep_orient){
                data[i]=std::tuple<int,int,bool>(current_id,position,std::get<2>(data[i]));//update the values 
            }
            else{
                std::tuple<uint64_t,bool> kmer_data=reverseComplementCanonical(k_1_mer,k-1);
                data[i]=std::tuple<int,int,bool>(current_id,position,std::get<1>(kmer_data));//update the values 
            }
            break;//break since the (k-1)-mer occurs once in a unitig except for palindromic k-mers
        }
    }
    index_table[canonical_bits(k_1_mer,k-1)]=data;
}
void Index::insert(uint64_t k_1_mer,int id,int position){
    /*
        the input is a string of size (k-1)
        the id of the unitig where this string occurs
        the position

        The k-1-mer gets inserted into the index
    */
   std::vector<std::tuple<int,int,bool>> occurences=find(k_1_mer);//get previous occurences
   std::tuple<uint64_t,bool> kmer_data=reverseComplementCanonical(k_1_mer,k-1);
   occurences.push_back(std::tuple<int,int,bool>(id,position,std::get<1>(kmer_data)));//insert the new occurence into vector
   index_table[std::get<0>(kmer_data)]=occurences;

}
void Index::insertSubUnitig(Unitig unitig,int id,int starting_position,int ending_position){
    /*
        the inputs are:
                     a string unitig
                     its id
                     the starting position (int)
                     the ending position (int)
        Actions performed: all the (k-1)-mers occuring between starting_position and ending position gets inserted
    */
   uint64_t i_th_mer=unitig.get_ith_mer(starting_position,k-1);
   insert(i_th_mer,id,starting_position);
   for(int position=starting_position+1;position<=ending_position;position++){
        //call the function insert to insert the (k-1)-mer at position i
        i_th_mer=unitig.get_next_mer(i_th_mer,position,k-1);
        insert(i_th_mer,id,position);
   }
}
void Index::update_unitig(Unitig seq,int current_id,int previous_id,int starting_position,int ending_position,bool keep_orient){
    /*
        the inputs are:
                     a string unitig
                     its current id
                     its previous id
                     the starting position (int)
                     the ending position (int)
        Actions performed: the info all the (k-1)-mers occuring between starting_position and ending position get updated
    */
    uint64_t i_th_mer=seq.get_ith_mer(starting_position,k-1);
    update_k_1_mer(i_th_mer,previous_id,current_id,starting_position,keep_orient);
    for(int position=starting_position+1;position<=ending_position;position++){
        i_th_mer=seq.get_next_mer(i_th_mer,position,k-1);
        update_k_1_mer(i_th_mer,previous_id,current_id,position,keep_orient);
    }
}
void Index::serialize(const std::string filename){
  std::ofstream ofs(filename, std::ios::binary);
  boost::archive::binary_oarchive oa(ofs);
  oa << index_table;
  ofs.close();
}
void Index::deserialize(const std::string filename){
  std::ifstream ifs(filename, std::ios::binary);
  boost::archive::binary_iarchive ia(ifs);
  ia >> index_table;
  ifs.close();
}
