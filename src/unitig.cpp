#include "unitig.h"
#include "CommonUtils.h"
#include <string_view>
#include<algorithm> 
#include<iterator>
#include <cstdint>
#include <string>
#include <cmath>
Unitig::Unitig(std::string_view unitig){
    left_unused_bits=0;//when storing 8 bits for every other 4 bases, the left unused bits is 0 since we store from left to right
    //encode every 4 bases
    for(int i=0;i+4<=unitig.length();i+=4){
        unitig_encoding.push_back(bit_encoding(unitig.substr(i,4)));
    }
    //encoding the remaining bases, a substring of 0<length<4
    if(unitig.length()%4!=0){
        unitig_encoding.push_back(bit_encoding(unitig.substr(unitig.length()-unitig.length()%4)));
        right_unused_bits=(4-unitig.length()%4)*2;//the right unused bits is the number of remaining characters times 2
    }
    else{
        right_unused_bits=0;//0 since we don't have remaining characters
    }


}
Unitig::Unitig(uint8_t left_unused,uint8_t right_unused, std::vector<uint8_t> encoding){
    /*
    constructor to create a unitig from an already computed vector of bits
    useful when we need to split a unitig at a position i
    */
    left_unused_bits=left_unused;//set left unused bits
    right_unused_bits=right_unused;//set right unused bits
    std::copy(encoding.begin(),encoding.end(),std::back_inserter(unitig_encoding));//copy element from vector encoding to unitig_encoding
}
uint8_t Unitig::get_left_unused(){
    return left_unused_bits;
}
uint8_t Unitig::get_right_unused(){
    return right_unused_bits;
}
int Unitig::unitig_length(){
    return 4*unitig_encoding.size()-(int)(left_unused_bits+right_unused_bits)/2;//number of elements in the bit vector * 4 bases - number of unused bits/2
}
std::vector<uint8_t> Unitig::get_encoding(){
    return unitig_encoding;
}
std::string Unitig::to_string(){
    /*
    recompute the string from the vector of bits
    useful to output the graph
    special care for the first and last encoding as they contain unused characters
    */
   std::string unitig="";
   if(left_unused_bits==0){
    unitig+=bits_to_seq_4(unitig_encoding[0],4);
   }
   else{
    unitig+=bits_to_seq_4(unitig_encoding[0],4).substr(left_unused_bits/2);
   }
   for(size_t i=1;i<unitig_encoding.size()-1;i++){
    unitig=unitig+bits_to_seq_4(unitig_encoding[i],4);
   }
   unitig=unitig+bits_to_seq_4(unitig_encoding[unitig_encoding.size()-1],4).substr(0,4-right_unused_bits/2);
   return unitig;
}
bool Unitig::operator==(const Unitig& other) const {
    // Check if member variables are equal
    // ...

    // Return true if equal, false otherwise
    //two unitigs are equal if they have the same left and right unused bits
    return left_unused_bits==other.left_unused_bits && right_unused_bits==other.right_unused_bits && unitig_encoding.size()==other.unitig_encoding.size() &&std::equal(unitig_encoding.begin(), unitig_encoding.end(),other.unitig_encoding.begin());
}
uint64_t Unitig::get_ith_mer(const int i,const int k){
    /*
        return the bit coding of the k-1 mer starting at i
    */
   int j = std::max(0, static_cast<int>(std::ceil((i-4 + (int)(left_unused_bits/2))/4)));//the position of the bit encoding of 4-mer in the vector which contains the letter i
   int pos_letter_in_bits=std::max(0,static_cast<int>(std::ceil((i-4 + (int)(left_unused_bits/2))%4)));//the position of ith letter in the bit encoding
   uint64_t answer=unitig_encoding[j]>>(2*pos_letter_in_bits);

   int letter_encoded=4-pos_letter_in_bits;
   while(letter_encoded<=k && j<unitig_encoding.size()){
    j=j+1;
    answer=(answer<<8)|unitig_encoding[j];
    letter_encoded+=4;
   }
}
uint64_t Unitig::get_next_mer(uint64_t mer,const int i,const int k){
    /*
    given the encoding of the (k-1)-mer at position i-1, and the position i we need to get the one at i
    */
   int j = std::max(0, static_cast<int>(std::ceil((i+k-5 + (int)(left_unused_bits/2))/4)));//the position of the bit encoding of 4-mer in the vector which contains the letter i
   int pos_letter_in_bits=std::max(0,static_cast<int>(std::ceil((i+k-5 + (int)(left_unused_bits/2))%4)));//the position of ith letter in the bit encoding
   uint8_t base=(unitig_encoding[j]<<2*pos_letter_in_bits)>>6;//shift by double the pos of the letter to the left then by 6, if we have 10110011 and pos=3 then we get 00000011
   uint64_t answer=((mer<<2)|base)&(0xffffffff>>(64-2*k));//set the rightmost two bits then mask the left most unused bits
   return answer;

}
void Unitig::insert_back(char base){
    /*
    add one base to the right of the unitig
    */
    uint8_t base_encoding=bit_encoding(std::string_view(std::string(1,base)));
    if(right_unused_bits!=0){//can we add directly to the last element
        base_encoding=base_encoding>>(8-right_unused_bits);
        unitig_encoding[unitig_encoding.size()-1]=(unitig_encoding[unitig_encoding.size()-1]<<right_unused_bits)|base_encoding;
        right_unused_bits-=2;
    }
    else{
        unitig_encoding.push_back(base_encoding);//otherwise insert it
        right_unused_bits=6;
    }
}
void Unitig::insert_front(char base){
    /*
    add one base at the beginning 
    */
    uint8_t base_encoding=bit_encoding(std::string_view(std::string(1,base)));
    if(left_unused_bits!=0){//can we add it directly
        base_encoding=base_encoding>>(left_unused_bits-2);
        unitig_encoding[0]=(unitig_encoding[0]&(0xff>>left_unused_bits))|base_encoding;//mask unnecessary bits
        left_unused_bits-=2;
    }
    else{
        unitig_encoding.emplace(unitig_encoding.begin(), base_encoding>>6);//otherwise add the encoding at the beginning
        left_unused_bits=6;
    }
}