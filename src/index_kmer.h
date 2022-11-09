#include <string>
#include "CommonUtils.h"
#include "ParseGFA.h"
#include <vector>
#include <unordered_map>

#ifndef index_kmer_H
#define index_kmer_H
class Index
{
private:
    /* data */
    int number_of_buckets;//the maximum number of buckets

    /*
    we use a map<(k-1)-mer,a vector tuple<unitig id,position orientation>> to index all the unitigs
    using a tuple to store the index <unitig id,position,orientation> of each (k-1)-mer

    */
    std::unordered_map<std::string,std::vector<std::tuple<int,int,bool>>> index_table;
    
public:
    //initialize the index
    Index(int number_of_buckets);

    //create it from the unitigs of the graph
    void create(GfaGraph & graph,int k);

    //this function returns the index of a (k-1)-mer, its all occurences which is equal to 8 at max
    std::vector<std::tuple<int,int,bool>> find(std::string k_1_mer);

    //update the table for one (k-1)-mer i.e change the id and position of the unitig
    void update_k_1_mer(std::string k_1_mer,int prev_id,int current_id,int position);

    //updating for all the (k-1)-mers of the unitig
    void update_unitig(std::string seq,int id,int previous_id,int k);

    //adding a new created unitig
    void add_unitig(std::string unitig,int id,int k);
    ~Index()=default;
};
#endif // !index_kmer_H