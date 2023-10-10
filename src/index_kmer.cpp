#include <zstd.h>
#include "index_kmer.h"
#include "CommonUtils.h"
#include <algorithm>
#include <cstring>
#include <vector>
#include <string>
#include <tuple>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <limits>
#include "FileSerializer.hpp"
template <typename T>
class CustomAllocator {
public:
    using value_type = T;

    T* allocate(std::size_t n) {
        if (n > std::numeric_limits<std::size_t>::max() / sizeof(T))
            throw std::bad_alloc();
        return static_cast<T*>(std::malloc(n * sizeof(T)));
    }

    void deallocate(T* p, std::size_t n) noexcept {
        std::free(p);
    }
};
using LargeVector = std::vector<char, CustomAllocator<char>>;
Index::Index(int k_size){
    k=k_size;
}
int Index::get_k(){
    return k;
}
void Index::create(GfaGraph& graph){
    // add description of function
   //index_table.set_empty_key(NULL);
    for(std::pair<int,Unitig> node:graph.get_nodes()){
       uint64_t i_th_mer=node.second.get_ith_mer(0,k-1);
       std::tuple<uint64_t,bool> seq_data=reverseComplementCanonical(i_th_mer,k-1);
       uint64_t position=(uint64_t)node.first;
       position=(position<<32)|std::get<1>(seq_data);//add clarification
       kmer_occurences[find_which_table(std::get<0>(seq_data))][std::get<0>(seq_data)]=position;
       for(int i=1;i<=node.second.unitig_length()-k+1;i++){
            i_th_mer=node.second.get_next_mer(i_th_mer,i,k-1);
            seq_data=reverseComplementCanonical(i_th_mer,k-1);
            position=(uint64_t)node.first;
            position=(position<<32)|(i<<1)|std::get<1>(seq_data);
            kmer_occurences[find_which_table(std::get<0>(seq_data))][std::get<0>(seq_data)]=position;
       }
    }
}
void Index::update_k_1_mer(uint64_t k_1_mer,int prev_id,int current_id,int position,bool keep_orient){
    /*
        the input are: a string of size (k-1)
        the previous id where this string occurs
        the current id where it occurs after the updte performed
        the new position where it occurs
    */
   std::tuple<uint64_t,bool> kmer_data=reverseComplementCanonical(k_1_mer,k-1);
   uint64_t kmer_data_bits=(uint64_t)current_id;
   kmer_data_bits=(kmer_data_bits<<32)|(position<<1)|std::get<1>(kmer_data);
   kmer_occurences[find_where_update(std::get<0>(kmer_data),prev_id)][std::get<0>(kmer_data)]=kmer_data_bits;//find where to update and update it
}
void Index::insert(uint64_t k_1_mer,int id,int position){
    /*
        the input is a string of size (k-1)
        the id of the unitig where this string occurs
        the position

        The k-1-mer gets inserted into the index
    */
   std::tuple<uint64_t,bool> kmer_data=reverseComplementCanonical(k_1_mer,k-1);
   uint64_t kmer_data_bits=(uint64_t)id;
   kmer_data_bits=(kmer_data_bits<<32)|(position<<1)|std::get<1>(kmer_data);
   kmer_occurences[find_which_table(std::get<0>(kmer_data))][std::get<0>(kmer_data)]=kmer_data_bits;//find where to insert and insert it

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
    std::ofstream outFile(filename, std::ios::binary);
    for (const auto& vec : kmer_occurences) {
        std::vector<char> buffer;
        size_t vSize=vec.size();
        outFile.write(reinterpret_cast<char*>(&vSize), sizeof(vSize));
        const char* vecData = reinterpret_cast<const char*>(vec.values().data());
        size_t vecSize = vec.size() * sizeof(std::pair<uint64_t, uint64_t>);
        buffer.insert(buffer.end(), vecData, vecData + vecSize);
        std::vector<char> compressedData(ZSTD_compressBound(buffer.size()));
        size_t compressedSize = ZSTD_compress(compressedData.data(), compressedData.size(),
                                  buffer.data(), buffer.size(), 1);
        compressedData.resize(compressedSize);
        outFile.write(reinterpret_cast<char*>(&compressedSize), sizeof(compressedSize));
        // Write the compressed data to a file
        outFile.write(compressedData.data(), compressedData.size());
    }
    
    outFile.close();
}
void Index::deserialize(const std::string filename){
    std::ifstream inFile(filename, std::ios::binary);
    // Read the sizes of each unordered_map from the input file
    size_t size[8];
    for (int i = 0; i < 8; i++) {
        size_t vecSize;
        inFile.read(reinterpret_cast<char*>(&vecSize), sizeof(vecSize));
        size_t compressedSize;
        inFile.read(reinterpret_cast<char*>(&compressedSize), sizeof(compressedSize));
        LargeVector compressedData(static_cast<size_t>(compressedSize));
        inFile.read(compressedData.data(), compressedSize);

    // Decompress the compressed buffer using Zstd
        size_t decompressedSize = ZSTD_getFrameContentSize(compressedData.data(), compressedData.size());
        LargeVector decompressedData(decompressedSize);
        decompressedSize = ZSTD_decompress(decompressedData.data(), decompressedData.size(), compressedData.data(), compressedData.size());
        //kmer_occurences[i].reserve(vecSize);
        const std::pair<uint64_t, uint64_t>* mapData = reinterpret_cast<const std::pair<uint64_t, uint64_t>*>(decompressedData.data());
        for(int j=0;j<vecSize;j++){
            kmer_occurences[i].emplace(mapData->first, mapData->second);
            mapData++;
        }
    }
    inFile.close();
}
size_t Index::find_which_table(uint64_t kmer){
    /*
    takes a (k-1)-mer as input then return where it should be stored
    */
    size_t where=0;

    while(where!=8){
        if(kmer_occurences[where].count(kmer)==0){
            return where;
        }
        where++;
    }
    return where;
}
size_t Index::how_many(uint64_t k_1_mer){
    /*
    takes a (k-1)-mer and return how many occurences of it
    it is the same analogy as find which table.
    If it should be inserted at the table whose index is 0, then we have zero occurences.
    the maximum value returned by the function find_which_table is 8
    */
    return find_which_table(canonical_bits(k_1_mer,k-1));
}
uint64_t Index::find_data(uint64_t k_1_mer,size_t i){
    /*
    This is useful for split and join
    in such a case the (k-1)-mer exist only once so it's in the
    */
   uint64_t canonical_kmer=canonical_bits(k_1_mer,k-1);
    if(kmer_occurences[i].count(canonical_kmer)==0){
        return 0;
    }
    return kmer_occurences[i].at(canonical_kmer);
}
size_t Index::find_where_update(uint64_t kmer,int id){
    /*
    takes a k-1 mer and the id of the unitig
    it finds where the update should occur
    */
    size_t where=0;
    while(where<8){
        uint64_t data=find_data(kmer,where);
        data=data>>32;
        if(data==id){
            return where;
        }
        where++;
    }
    return where;
}