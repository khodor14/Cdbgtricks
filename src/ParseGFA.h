#include <vector>
#include <fstream>
#include <string>
#include <vector>
#include <unordered_map>
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
    std::unordered_map<int,std::string> unitigs;
    std::vector<edge> edges;
    int max_node_id=0;
public:
    GfaGraph() = default;
    static GfaGraph LoadFromFile(std::string filename);
    static GfaGraph LoadFromStream(std::ifstream & file,bool gfa);
    void convertToFasta(std::string filename);
    std::vector<int> find_in_neighbors(int node_id);
    std::vector<int> find_out_neighbors(int node_id);
    std::unordered_map<int,std::string> get_nodes();
    void fixe_edges(int node_id,int new_node, bool from, bool to);
    int get_max_node_id();
    void set_max_node_id(int id);
    ~GfaGraph() = default;
};

#endif // !ParseGFA_H

