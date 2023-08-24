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
class ColorSet{
    private:
        ankerl::unordered_dense::map<int,std::string> colors_names;
        std::unordered_map<int,BitVector> color_classes;
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
};
#endif // !COLORSET
