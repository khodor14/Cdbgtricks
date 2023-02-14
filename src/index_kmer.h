#include <string>
#include "CommonUtils.h"
#include "ParseGFA.h"
#include <vector>
#include <unordered_map>
#include <sparsehash/sparse_hash_map>
#include <tuple>
#ifndef index_kmer_H
#define index_kmer_H
class Index
{
private:
    /* data */
    int number_of_buckets;//the maximum number of buckets
    int k;//the size of the k-mer
    /*
    we use a map<(k-1)-mer,a vector tuple<unitig id,position orientation>> to index all the unitigs
    using a tuple to store the index <unitig id,position,orientation> of each (k-1)-mer

    */
    google::sparse_hash_map<std::string,std::vector<std::tuple<int,int,bool>>> index_table;

    //std::unordered_map<std::string,std::vector<std::tuple<int,int,bool>>> index_table;
    
public:
    //initialize the index
    Index(int number_of_buckets,int k);

    //create it from the unitigs of the graph
    void create(GfaGraph & graph);
    int get_k();
    //this function returns the index of a (k-1)-mer, its all occurences which is equal to 8 at max
    std::vector<std::tuple<int,int,bool>> find(std::string k_1_mer);

    //update the table for one (k-1)-mer i.e change the id and position of the unitig
    void update_k_1_mer(std::string k_1_mer,int prev_id,int current_id,int position);

    //insert a (k-1)-mer in unitig u of id id at the given position into the index
    void insert(std::string k_1_mer,int id,int position);
    //insert all (k-1)-mers in unitig from starting till ending
    void insertSubUnitig(std::string unitig,int id,int starting_position,int ending_position);
    //updating for all the (k-1)-mers of the unitig
    void update_unitig(std::string seq,int id,int previous_id,int starting_position,int ending_position);
    void serialize(const std::string filename);
    void deserialize(const std::string filename);
    ~Index()=default;
};
#endif // !index_kmer_H