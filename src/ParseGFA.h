#ifndef GfaGraph
#define GfaGraph
#include <iostream>
#include <vector>
#include <string>

#ifndef structNode
#define structNode
struct structNode
{
    uint id;
    std::string sequence; 
};


#endif // !structNode

#ifndef edge
#define edge
struct edge
{
    u_int sink;
    u_int source;
    bool from(true);
    bool to(true);
};

#endif // !edge
class GfaGraph{
public:
    GfaGraph LoadFromFile(std::string filename);
    void convertToFasta();
    std::vector<structNode> nodes;
    std::vector<edge> edges;
}
#endif // !GfaGraph
