#include "bit_vector.hpp"
#include <iostream>
bool BitVector::operator==(const BitVector& other) const{
    return bit_array==other.bit_array;
}
BitVector::BitVector(size_t num_of_bits){
    number_of_bits=num_of_bits;
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
    std::cout<<elem<<"\n";
    return elem&(1ULL<<(63-pos%64));
}
size_t BitVector::get_number_of_bits(){
    return number_of_bits;
}
std::vector<uint64_t> BitVector::getColors(){
    return bit_array;
}