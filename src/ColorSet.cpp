#include "Colorset.hpp"
template <typename T>
class CustomAllocator {
public:
    using value_type = T;

    T* allocate(std::size_t n) {
        if (n > std::numeric_limits<std::size_t>::max() / sizeof(T))
            throw std::bad_alloc();
        return static_cast<T*>(std::malloc(n * sizeof(T)));
    }

    void deallocate(T* p, std::size_t n) noexcept {
        std::free(p);
    }
};
using LargeVector = std::vector<char, CustomAllocator<char>>;
void ColorSet::insert_color_class_name(int id,std::string color_name){
    colors_names[id]=color_name;
}
void ColorSet::insert_Color_class(int id,BitVector& colors){
    color_classes[id]=colors;
}
BitVector ColorSet::get_color_class(int id){
    if(color_classes.count(id)>0){
        return color_classes[id];
    }
    return BitVector();
}
std::string ColorSet::get_color_class_name(int id){
    if(colors_names.count(id)>0){
        return colors_names[id];
    }
}
void ColorSet::read(std::string colors_f_name){
    //start reading color names
    std::ifstream inFile(colors_f_name, std::ios::binary);
    size_t vecSize;
    inFile.read(reinterpret_cast<char*>(&vecSize), sizeof(vecSize));
    size_t compressedSize;
    inFile.read(reinterpret_cast<char*>(&compressedSize), sizeof(compressedSize));
    LargeVector compressedData(static_cast<size_t>(compressedSize));
    inFile.read(compressedData.data(), compressedSize);
    // Decompress the compressed buffer using Zstd
    size_t decompressedSize = ZSTD_getFrameContentSize(compressedData.data(), compressedData.size());
    LargeVector decompressedData(decompressedSize);
    decompressedSize = ZSTD_decompress(decompressedData.data(), decompressedData.size(), compressedData.data(), compressedData.size());
        //kmer_occurences[i].reserve(vecSize);
    const std::pair<int, std::string>* mapData = reinterpret_cast<const std::pair<int, std::string>*>(decompressedData.data());
    size_t numPairs = vecSize/sizeof(std::pair<int, std::string>);
    for(size_t i=0;i<numPairs;i++){
        colors_names[mapData->first]=mapData->second;
        mapData++;
    }
    //end reading color names
    //start reading color classes
    size_t color_numbers;
    inFile.read(reinterpret_cast<char*>(&color_numbers), sizeof(color_numbers));
    for(size_t i=0;i<color_numbers;i++){
        int id;
        inFile.read(reinterpret_cast<char*>(&id), sizeof(id));
        size_t number_of_bits;
        inFile.read(reinterpret_cast<char*>(&number_of_bits), sizeof(number_of_bits));
        size_t vecSize;
        inFile.read(reinterpret_cast<char*>(&vecSize), sizeof(vecSize));
        size_t compressedSize_vec;
        inFile.read(reinterpret_cast<char*>(&compressedSize_vec), sizeof(compressedSize_vec));
        std::vector<char> compressedData_vec(compressedSize_vec);
        inFile.read(compressedData_vec.data(), compressedSize_vec);
         const size_t decompressedBufferSize = ZSTD_getFrameContentSize(compressedData_vec.data(), compressedSize_vec);
        std::vector<char> decompressedBuffer(decompressedBufferSize);
        decompressedSize = ZSTD_decompress(decompressedBuffer.data(), decompressedBuffer.size(), compressedData_vec.data(), compressedData_vec.size());
        std::vector<uint64_t> decompressedData(decompressedSize / sizeof(uint64_t));
        std::memcpy(decompressedData.data(), decompressedBuffer.data(), decompressedSize);
    }
}
void ColorSet::write(std::string colors_f_name){
    std::ofstream outFile(colors_f_name, std::ios::binary);
    //write color names
    std::vector<char> buffer_names;
    size_t size=colors_names.size();
    outFile.write(reinterpret_cast<char*>(&size), sizeof(size));
    const char* vecData_names = reinterpret_cast<const char*>(colors_names.values().data());
    size_t vecSize_names = size * sizeof(std::pair<int, std::string>);
    buffer_names.insert(buffer_names.end(), vecData_names, vecData_names + vecSize_names);
    std::vector<char> compressedData_names(ZSTD_compressBound(buffer_names.size()));
    size_t compressedSize = ZSTD_compress(compressedData_names.data(), compressedData_names.size(),
                                  buffer_names.data(), buffer_names.size(), 1);
    compressedData_names.resize(compressedSize);
    outFile.write(reinterpret_cast<char*>(&compressedSize), sizeof(compressedSize));
    // Write the compressed data to a file
    outFile.write(compressedData_names.data(), compressedData_names.size());
    //write color_classes
    size_t color_colass_numbers=color_classes.size();
    outFile.write(reinterpret_cast<char*>(&color_colass_numbers), sizeof(color_colass_numbers));
    for (auto& vec : color_classes) {
        std::vector<char> buffer;
        std::vector<uint64_t> colors=vec.second.getColors();
        size_t number_of_Bits=vec.second.get_number_of_bits();
        size_t vSize=colors.size();
        int id=vec.first;
        //write color class id
        outFile.write(reinterpret_cast<char*>(&id), sizeof(id));
        //write number of bits
        outFile.write(reinterpret_cast<char*>(&number_of_Bits), sizeof(number_of_Bits));
        //write size of color vector
        outFile.write(reinterpret_cast<char*>(&vSize), sizeof(vSize));
        //compress colors and write to disk
        const char* vecData = reinterpret_cast<const char*>(colors.data());
        size_t vecSize = colors.size() * sizeof(uint64_t);
        buffer.insert(buffer.end(), vecData, vecData + vecSize);
        std::vector<char> compressedData(ZSTD_compressBound(buffer.size()));
        compressedSize = ZSTD_compress(compressedData.data(), compressedData.size(),
                                  buffer.data(), buffer.size(), 1);
        compressedData.resize(compressedSize);
        outFile.write(reinterpret_cast<char*>(&compressedSize), sizeof(compressedSize));
        // Write the compressed data to a file
        outFile.write(compressedData.data(), compressedData.size());
    }
    
    outFile.close(); 
}
int ColorSet::get_highest_color_class_id(){
    return color_classes.size();
}
int ColorSet::get_number_color_names(){
    return colors_names.size();
}