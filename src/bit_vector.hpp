#ifndef BVECTOR
#define BVECTOR
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <cstring>
class BitVector{
    private:
        size_t number_of_bits;
        int number_of_bits_set;
        int number_left_bits;
        int number_right_bits;
        std::vector<uint64_t> bit_array;
    public:
        BitVector()= default;
        BitVector(size_t number_bits);
        BitVector(size_t number_bits,std::vector<uint64_t>& colors);
        void set(size_t pos);
        void unset(size_t pos);
        void reset();
        void set_number_bits(int num_bits);
        bool is_set(size_t pos);
        std::vector<uint64_t> getColors();
        size_t get_number_of_bits();
        bool operator==(const BitVector& other) const;
        uint64_t get_identity();
        inline bool operator!=(const BitVector& other)const { 
            return !operator==(other);
        }
};
#endif // !BVECTOR
