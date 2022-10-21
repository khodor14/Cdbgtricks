#include <string>
#include <ParseGFA.h>
#include <vector>
#include <unordered_map>
class Index
{
private:
    /* data */
    int number_of_buckets;
    std::unordered_map<std::string,std::vector<std::tuple<int,int,bool>>> index_table;
    
public:
    Index(int number_of_buckets);
    void create(GfaGraph & graph,int k);
    std::vector<std::tuple<int,int,bool>> find(std::string k_1_mer);
    void update_k_1_mer(std::string k_1_mer,int prev_id,int current_id,int position);
    void update_unitig(Node & unitig,int previous_id,int k);
    ~Index();
};

