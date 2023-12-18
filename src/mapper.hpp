#ifndef mapper_H
#define mapper_H
#include <string>
#include "CommonUtils.h"
#include "ParseGFA.h"
#include <vector>
#include <string_view>
#include "index.hpp"
class Mapper
{
private:
    float ratio;
public:
    Mapper(float r);    
    std::vector<int> map_back(const std::string_view seq,Index_mphf& graph_ind,GfaGraph& graph);
    ~Mapper()=default;
};
#endif // !index_mphf_H