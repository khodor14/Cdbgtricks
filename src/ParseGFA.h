class Node
{
private:
    int id;
    std::string unitig;
public:
    Node(int id,std::string unitig);
    ~Node();
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
    edge(/* args */);
    ~edge();
};

class GfaGraph
{
private:
    std::vector<Node> nodes;
    std::vector<edge> edges;
public:
    GfaGraph(/* args */);
    static GfaGraph LoadFromFile(std::string filename);
    void convertToFasta();
    ~GfaGraph();
};

