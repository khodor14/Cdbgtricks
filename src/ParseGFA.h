#include <vector>
#include <fstream>
#include <string>
#include <vector>
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
    std::vector<Node> nodes;
    std::vector<edge> edges;
public:
    GfaGraph() = default;
    static GfaGraph LoadFromFile(std::string filename);
    static GfaGraph LoadFromStream(std::ifstream & file);
    void convertToFasta();
    std::vector<int> find_in_neighbors(int node_id);
    std::vector<int> find_out_neighbors(int node_id);
    std::vector<Node> get_nodes();
    void fixe_edges(int node_id,int new_node, bool from, bool to);
    void fixe_node(int node_id,int pos_modification);
    int create_nodes(std::string sequence);//return the node id
    void create_edge(int node_id,int new_node_id,bool from,bool to);
    ~GfaGraph() = default;
};

#endif // !ParseGFA_H

