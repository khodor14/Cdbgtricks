#include "bit_vector.hpp"
#include <iostream>
bool BitVector::operator==(const BitVector& other) const{
    return bit_array==other.bit_array;
}
BitVector::BitVector(size_t num_of_bits){
    number_of_bits=num_of_bits;
    number_of_bits_set=0;
    number_left_bits=0;
    number_right_bits=0;
    int capacity=ceil(num_of_bits/64);
    bit_array.resize(capacity);
}
BitVector::BitVector(size_t number_bits,std::vector<uint64_t>& colors){
    number_of_bits=number_bits;
    bit_array.resize(colors.size());
    for(int i=0;i<colors.size();i++){
        bit_array[i]=colors[i];
    }
}
void BitVector::set(size_t pos){
    int pos_in_array=pos/64;
    number_of_bits_set+=1;
    //computing left and right bits
    if(pos*2<number_of_bits){
        number_left_bits=number_left_bits+1;
    }
    else{
        number_right_bits=number_right_bits+1;
    }
    if(pos_in_array>=bit_array.size()){
        for(int j=0;j<pos_in_array-bit_array.size()+1;j++){
            bit_array.push_back(0);
        }
    }
    uint64_t elem=bit_array[pos_in_array];
    elem=elem|((uint64_t)1<<(63-pos%64));
    bit_array[pos_in_array]=elem; 
}
void BitVector::unset(size_t pos){
    int pos_in_array=pos/64;
    uint64_t elem=bit_array[pos_in_array];
    uint64_t bitmask = 1ULL <<(63-pos%64);
    bitmask=~bitmask;
    bit_array[pos_in_array]=elem&bitmask;

}
void BitVector::reset(){
    for(size_t i=0;i<bit_array.size();i++){
        bit_array[i]=0;
    }
}
bool BitVector::is_set(size_t pos){
    int pos_in_array=pos/64;
    uint64_t elem=bit_array[pos_in_array];
    return elem&(1ULL<<(63-pos%64));
}
size_t BitVector::get_number_of_bits(){
    return number_of_bits;
}
std::vector<uint64_t> BitVector::getColors(){
    return bit_array;
}
uint64_t BitVector::get_identity(){
    uint64_t res=number_of_bits;
    res=res<<32;
    int min,max;
    if(number_left_bits<number_right_bits){
        min=number_left_bits;
        max=number_right_bits;
    }
    else{
        min=number_right_bits;
        max=number_left_bits;
    }
    int diff_Ratio=max-min+max/min;

    res=res|diff_Ratio;
    return res;
}