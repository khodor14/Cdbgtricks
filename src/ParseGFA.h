#ifndef ParseGFA_H
#define ParseGFA_H
#include <vector>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
#include "unitig.h"
#include "CommonUtils.h"
#include "kseq.h"
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cstring>
#include <string_view>
#include <sys/stat.h>
#include <stdint.h>
#include <stdio.h>
#include <memory>
#include "zstr.hpp"

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
    int test_kmer_presence(uint64_t kmer,uint64_t suffix,uint64_t kmer_pos,int k);
    bool test_kmer(uint64_t kmer,uint64_t kmer_pos,int k);
    uint64_t get_kmer(uint64_t position,int k);
    void serialize(const std::string filename);
    void deserialize(const std::string filename);
    ~GfaGraph() = default;
};

#endif // !ParseGFA_H

