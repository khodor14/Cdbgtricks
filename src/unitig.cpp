#include "unitig.h"
#include "CommonUtils.h"
#include <string_view>
#include<algorithm> 
#include<iterator>
Unitig::Unitig(std::string unitig){
    left_unused_bits=0;//when storing 8 bits for every other 4 bases, the left unused bits is 0 since we store from left to right
    right_uused_bits=(uint8_t)((unitig.length()%4)*2);//the right unused bits is the number of remaining characters times 2
    std::string_view copy{unitig};
    //encode every 4 bases
    for(int i=0;i<(int)(copy.length()/4);i+=4){
        unitig_encoding.push_back(bit_ecoding(copy.substr(i,4)));
    }
    //encoding the remaining bases, a substring of 0<length<4
    if(copy.length()%4!=0){
        unitig_encoding.push_back(bit_ecoding(copy.substr(copy.length()-copy.length()%4)));
    }


}
Unitig::Unitig(uint8_t left_unused,uint8_t right_unused, std::vector<uint8_t> encoding){
    /*
    constructor to create a unitig from an already computed vector of bits
    useful when we need to split a unitig at a position i
    */
    left_unused_bits=left_unused;//set left unused bits
    right_uused_bits=right_unused;//set right unused bits
    std::copy(encoding.begin(),encoding.end(),std::back_inserter(unitig_encoding));//copy element from vector encoding to unitig_encoding
}
uint8_t Unitig::get_left_unused(){
    return left_unused_bits;
}
uint8_t Unitig::get_right_unused(){
    return right_uused_bits;
}
int Unitig::unitig_length(){
    return 4*unitig_encoding.size()-(int)(left_unused_bits+right_uused_bits)/2;//number of elements in the bit vector * 4 bases - number of unused bits/2
}
std::string Unitig::to_string(){
    /*
    recompute the string from the vector of bits
    useful to output the graph
    */
    return "";
}