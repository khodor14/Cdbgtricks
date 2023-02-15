#include <string>
#include "CommonUtils.h"
#include <vector>
#ifndef unitig_H
#define unitig_H
class Unitig
{
private:
    /* data */
    uint8_t left_unused_bits;//the leftmost unused bits in the first bit encoding element of the vector
    uint8_t right_uused_bits;//the righmost uused bits in the first bit encoding
    /*
    storing every 4 bases in uint8_t (8 bits)
    */
    std::vector<uint8_t> unitig_encoding

    //std::unordered_map<std::string,std::vector<std::tuple<int,int,bool>>> index_table;
    
public:
    //create it from a string, here the left unused  bits is 0
    Unitig(std::string unitig);

    //create it from 
    Unitig(uint8_t left_unused,uint8_t right_unused, std::vector<uint8_t> encoding);
    uint8_t get_left_unused();
    uint8_t get_right_unused();
    int unitig_length();
    std::string to_string();
    ~Unitig=default;
};
#endif // !index_kmer_H