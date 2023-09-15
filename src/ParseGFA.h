#include <vector>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include "unitig.h"
#include "FileSerializer.hpp"
#include "Colorset.hpp"
#include "bit_vector.hpp"
#ifndef ParseGFA_H
#define ParseGFA_H

class Node
{
private:
    int id;
    std::string unitig;
public:
    Node(int id,std::string unitig);
    int get_id();
    std::string get_unitig();
    ~Node() = default;
    bool operator!=(const Node& other);
};

class edge
{
private:
    int sink;
    int source;
    bool from=true;
    bool to=true;
public:
    edge(int sink,int source,bool from,bool to);
    int get_source();
    int get_sink();
    ~edge() = default;
    bool operator!=(edge &other);

};

class GfaGraph
{
private:
    std::unordered_map<int,Unitig> unitigs;
    std::vector<edge> edges;
    int max_node_id=0;
    ColorSet color_classes;
public:
    GfaGraph() = default;
    static GfaGraph LoadFromFile(std::string filename);
    static GfaGraph LoadFromStream(std::string file,bool gfa);
    void convertToFasta(std::string filename);
    std::vector<int> find_in_neighbors(int node_id);
    std::vector<int> find_out_neighbors(int node_id);
    std::unordered_map<int,Unitig> get_nodes();
    void fixe_edges(int node_id,int new_node, bool from, bool to);
    int get_max_node_id();
    void set_max_node_id(int id);
    void delete_unitig(int id);
    void insert_unitig(int id,Unitig u);
    Unitig get_unitig(int id);
    void insert_color_class(int id,BitVector& colors);
    void serialize(const std::string filename);
    void deserialize(const std::string filename);
    int get_highest_number_colors();
    void insert_color_name(int id,std::string color_name);
    int get_number_classes();
    void update_color_class_id(int id,int color_Class_id);
    BitVector get_color_class(int id);
    int get_color_class_of_unitig(int id_unitig);
    void write_colors(std::string filename);
    void read_colors(std::string filename);
    ~GfaGraph() = default;
    inline void insert_id_hash(uint64_t hash,int id){
        color_classes.insert_id_hash(hash,id);
    }
    inline bool hash_exist(uint64_t hash){
        return color_classes.hash_exist(hash);
    }
    inline std::vector<int> get_similar_vector(uint64_t hash){
        return color_classes.get_similar_vector(hash);
    }
    inline int same_vector(BitVector query,std::vector<int> ids){
        return color_classes.same_vector(query,ids);
    }
};

#endif // !ParseGFA_H

