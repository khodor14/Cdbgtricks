#ifndef COLORSET
#define COLORSET
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <string>
#include <ankerl/unordered_dense.h>
#include "bit_vector.hpp"
#include <unordered_map>
#include <zstd.h>
#include <iostream>
#include <fstream>
#include <vector>
class ColorSet{
    private:
        ankerl::unordered_dense::map<int,std::string> colors_names;
        std::unordered_map<int,BitVector> color_classes;
        std::unordered_map<uint64_t,std::vector<int>> hash_to_ids;
    public:
        ColorSet()= default;
        void insert_color_class_name(int id,std::string name);
        void insert_Color_class(int id,BitVector& color_class);
        BitVector get_color_class(int id);
        std::string get_color_class_name(int id);
        void read(std::string color_f_name);
        void write(std::string color_f_name);
        int get_highest_color_class_id();
        int get_number_color_names();
        inline void insert_id_hash(uint64_t hash,int id){
            std::vector<int> ids_vec=hash_to_ids[hash];
            ids_vec.push_back(id);
            hash_to_ids[hash]=ids_vec;
        }
        inline bool hash_exist(uint64_t hash){
            return hash_to_ids.count(hash)>0;
        }
        inline std::vector<int> get_similar_vector(uint64_t hash){
            return hash_to_ids[hash];
        }
        inline int same_vector(BitVector query,std::vector<int> ids){
            for(int id:ids){
                if(query==color_classes[id]){
                    return id;
                }
            }
            return -1;
        }
};
#endif // !COLORSET
