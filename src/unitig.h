#ifndef unitig_H
#define unitig_H
#include <vector>
class Unitig
{
private:
    /* data */
    uint8_t left_unused_bits;//the leftmost unused bits in the first bit encoding element of the vector
    uint8_t right_unused_bits;//the righmost unsed bits in the first bit encoding
    /*
    storing every 4 bases in uint8_t (8 bits)
    */
    std::vector<uint8_t> unitig_encoding;

    //std::unordered_map<std::string,std::vector<std::tuple<int,int,bool>>> index_table;
    
public:
    //create it from a string, here the left unused  bits is 0
    Unitig()= default;
    explicit Unitig(std::string unitig);
    //create it from 
    Unitig(uint8_t left_unused,uint8_t right_unused, std::vector<uint8_t> encoding);
    bool operator==(const Unitig& other) const;
    uint8_t get_left_unused();
    uint8_t get_right_unused();
    uint64_t get_ith_mer(int i,int k);
    uint64_t get_next_mer(uint64_t mer,int i,int k);
    void insert_front(char base);
    void insert_back(char base);
    int unitig_length();
    std::string to_string();
    std::vector<uint8_t> get_encoding();
    ~Unitig()=default;
};
#endif // !index_kmer_H