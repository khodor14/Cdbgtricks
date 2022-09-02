#include <CommonUtils.h>
#include <ParseGFA.h>
#include <iostream>
void search_modifications(std::string k_1_mer,std::vector<uint64_t> hash_kmers,std::vector<std::tuple<int,char,std::tuple<int,bool,bool>>> possible_modifications);
void modify_unitig(std::string unitig,std::vector<uint64_t> hash_kmers);
static GfaGraph modify_graph(GfaGraph g,std::vector<uint64_t> hash_kmers);
