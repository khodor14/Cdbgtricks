#include <string>
#include "CommonUtils.h"
#include "ParseGFA.h"
#include <vector>
#include <unordered_map>
#include <sparsehash/dense_hash_map>
#include <tuple>
#ifndef index_kmer_H
#define index_kmer_H
class Index
{
private:
    /* data */
    int k;//the size of the k-mer
    /*
    we use a map<(k-1)-mer,a vector tuple<unitig id,position orientation>> to index all the unitigs
    using a tuple to store the index <unitig id,position,orientation> of each (k-1)-mer

    */
   /*
    mask to get orientation: &1
    mask to get position: (>>1)&0x7FFFFFFF
    mask to get identity: >>32
   */
    std::unordered_map<uint64_t, uint64_t> kmer_occurences[8];//eight tables where table i contains the ith occurence of the kmer
    size_t find_which_table(const uint64_t kmer);//useful to insert in a table
    size_t find_where_update(const uint64_t kmer,int id);//where update should occur
    //std::unordered_map<std::string,std::vector<std::tuple<int,int,bool>>> index_table;
    
public:
    //initialize the index
    Index(int k);

    //create it from the unitigs of the graph
    void create(GfaGraph & graph);
    int get_k();
    //update the table for one (k-1)-mer i.e change the id and position of the unitig
    void update_k_1_mer(uint64_t k_1_mer,int prev_id,int current_id,int position,bool keep_orient);
    //insert a (k-1)-mer in unitig u of id id at the given position into the index
    void insert(uint64_t k_1_mer,int id,int position);
    //insert all (k-1)-mers in unitig from starting till ending
    void insertSubUnitig(Unitig unitig,int id,int starting_position,int ending_position);
    //updating for all the (k-1)-mers of the unitig
    void update_unitig(Unitig seq,int id,int previous_id,int starting_position,int ending_position,bool keep_orient);
    void serialize(const std::string filename);
    void deserialize(const std::string filename);
    size_t how_many(const uint64_t k_1_mer);//the number of occurrences of a (k-1)-mer
    uint64_t find_data(uint64_t k_1_mer,size_t i);//find the data in the i-th table
    ~Index()=default;
};
#endif // !index_kmer_H