#include "unitig.h"
#include "CommonUtils.h"
#include <string_view>
#include<algorithm> 
#include<iterator>
#include <cstdint>
#include <string>
#include <cassert>
#include <cmath>
Unitig::Unitig(std::string_view unitig){
    unused=0;//from left everything is used
    for(int i=0;i+4<=unitig.length();i+=4){
        unitig_encoding.push_back(bit_encoding(unitig.substr(i,4)));
    }
    //encoding the remaining bases, a substring of 0<length<4
    if(unitig.length()%4!=0){
        unitig_encoding.push_back(bit_encoding(unitig.substr(unitig.length()-unitig.length()%4)));
        unused=unused|(4-unitig.length()%4);
    }
}
Unitig::Unitig(uint8_t left_unused,uint8_t right_unused, std::vector<uint8_t> encoding){
    /*
    constructor to create a unitig from an already computed vector of bits
    useful when we need to split a unitig at a position i
    */
   unused=(left_unused<<2)|right_unused;
    std::copy(encoding.begin(),encoding.end(),std::back_inserter(unitig_encoding));//copy element from vector encoding to unitig_encoding
}
uint8_t Unitig::get_left_unused(){
    return unused>>2;
}
uint8_t Unitig::get_right_unused(){
    return unused&3;
}
int Unitig::unitig_length(){
    return 4*unitig_encoding.size()-(int)get_left_unused()-(int)get_right_unused();//number of elements in the bit vector * 4 bases - number of unused bits charaters from left and right
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
   char seq[unitig_length()+1];
   int i=0,k;
   std::string first_part=bits_to_seq_4(unitig_encoding[0],4).substr(get_left_unused());
   for(k=0;i<first_part.length();k++){
    seq[i]=first_part[k];
    i++;
   }
   //std::string unitig=bits_to_seq_4(unitig_encoding[0],4).substr(get_left_unused());//take the first characters
   for(size_t j=1;j<unitig_encoding.size()-1;j++){
    first_part=bits_to_seq_4(unitig_encoding[j],4);
    for(k=0;k<first_part.length();k++){
        seq[i]=first_part[k];
        i++;
    }
    //i++;
    //unitig=unitig+bits_to_seq_4(unitig_encoding[i],4);//take 4 characters at time and append
   }
   first_part=bits_to_seq_4(unitig_encoding[unitig_encoding.size()-1],4).substr(0,4-(int)get_right_unused());
   for(k=0;k<first_part.length();k++){
        seq[i]=first_part[k];
        i++;
    }
    seq[unitig_length()]='\0';
   //unitig=unitig+bits_to_seq_4(unitig_encoding[unitig_encoding.size()-1],4).substr(0,4-(int)get_right_unused());//add the last
   return seq;
}
bool Unitig::operator==(const Unitig& other) const {
    // Check if member variables are equal
    // ...

    // Return true if equal, false otherwise
    //two unitigs are equal if they have the same unused bits and same encoding of vectors
    //they can have different unused characters and different encoding, but the way we are encoding the unitigs, the split and join it cannot happen
    return unused==other.unused && unitig_encoding.size()==other.unitig_encoding.size() &&std::equal(unitig_encoding.begin(), unitig_encoding.end(),other.unitig_encoding.begin());
}
uint64_t Unitig::get_ith_mer(const int i,const int k){
    /*
        return the bit coding of the k-1 mer starting at i
    */
   int j = static_cast<int>(std::floor((i + (int)get_left_unused())/4));//the position of the bit encoding of 4-mer in the vector which contains the letter i
   int pos_letter_in_bits=(i+(int)get_left_unused())%4;//the position of ith letter in the bit encoding
   uint64_t answer=unitig_encoding[j]&(0xff>>(2*pos_letter_in_bits));//mask the left bits
   int letter_encoded=4-pos_letter_in_bits;
   
   while(letter_encoded+4<=k && j<unitig_encoding.size()){
    j=j+1;
    answer=(answer<<8)|unitig_encoding[j];
    letter_encoded+=4;
   }
   if(letter_encoded<k){
    j=j+1;
    answer=(answer<<(2*(k-letter_encoded)))|(unitig_encoding[j]>>(8-2*(k-letter_encoded)));
   }
   return answer;
}
uint64_t Unitig::get_next_mer(uint64_t mer,const int i,const int k){
    /*
    given the encoding of the (k-1)-mer at position i-1, and the position i we need to get the one at i
    */
   int j = static_cast<int>(std::floor((i+k-1 + (int)get_left_unused())/4));//the position of the bit encoding of 4-mer in the vector which contains the letter i+k-1 to be added to the  (i-1)'k-thmer
   int pos_letter_in_bits=(i+k-1+(int)get_left_unused())%4;//the position of ith letter in the bit encoding
   uint8_t base=(unitig_encoding[j]<<2*pos_letter_in_bits)>>6;//shift by double the pos of the letter to the left then by 6, if we have 10110011 and pos=3 then we get 00000011
   uint64_t answer=((mer<<2)|base)&(0xffffffffffffffff>>(64-2*k));//set the rightmost two bits then mask the left most unused bits
   return answer;

}
void Unitig::insert_back(char base){
    /*
    add one base to the right of the unitig
    */
    uint8_t base_encoding=baseToInt(base);
    uint8_t right_unused=get_right_unused();
    if(!(unused&3)){//can we add directly to the last element
        base_encoding=base_encoding<<((right_unused>>1)-2);
        unitig_encoding[unitig_encoding.size()-1]=(unitig_encoding[unitig_encoding.size()-1]&(0xff<<(get_right_unused()<<2)))|base_encoding;
        unused=unused-1;
    }
    else{
        unitig_encoding.push_back(base_encoding<<6);//otherwise insert it
        unused=unused|3;
    }
}
void Unitig::insert_front(char base){
    /*
    add one base at the beginning 
    */
    uint8_t base_encoding=baseToInt(base);
    if(!(unused&12&3)){//can we add it directly
        base_encoding=base_encoding<<(8-get_left_unused()<<2);
        unitig_encoding[0]=(unitig_encoding[0]&(0xff>>(get_left_unused()<<2)))|base_encoding;//mask unnecessary bits
        unused=unused-4;
    }
    else{
        unitig_encoding.emplace(unitig_encoding.begin(), base_encoding);//otherwise add the encoding at the beginning
        unused=unused|12;
    }
}