#include <index_kmer.h>
#include <CommonUtils.h>
#include <algorithm>
#include <cstring>

Index::Index(int buckets){
    number_of_buckets=buckets;
    index_table.rehash(number_of_buckets);
}

void Index::create(GfaGraph& graph,int k){
    std::vector<Node> nodes=graph.get_nodes();
    auto iter_nodes=nodes.begin();
    while(iter_nodes!=nodes.end()){
        int id=iter_nodes->get_id();
        std::string unitig=iter_nodes->get_unitig();

        for(int i=0;i<unitig.length()-k+2;i++){
            std::string sub=unitig.substr(i,i+k-1);
            std::string rc=getCanonical(sub);
            std::string toStore=std::min(sub,rc);
            bool orientation=true;
            if(std::strcmp(sub.c_str(),rc.c_str())>0){
                orientation=false;
            }
            if(index_table.find(toStore)==index_table.end()){
                    std::vector<std::tuple<int,int,bool>> data;
                    data.push_back(std::tuple<int,int,bool>(id,i,orientation));
                    index_table[toStore]=data;
            }
            else{
                std::vector<std::tuple<int,int,bool>> data=index_table.at(toStore);
                data.push_back(std::tuple<int,int,bool>(id,i,orientation));
                index_table[toStore]=data;

            }
        }
	}
}
std::vector<std::tuple<int,int,bool>> Index::find(std::string kmer){
    return index_table.at(getCanonical(kmer));
}
void Index::update_k_1_mer(std::string k_1_mer,int prev_id,int current_id,int position){
    std::vector<std::tuple<int,int,bool>> data=index_table.at(getCanonical(k_1_mer));
    for (std::tuple<int, int,bool>& tup : data){
        if(std::get<0>(tup)==prev_id){
            std::get<0>(tup)=current_id;
            std::get<1>(tup)=position;
            break;
        }
    }
}
void Index::update_unitig(std::string seq,int id,int previous_id,int k){

    for(int i=0;i<seq.length()-k+2;i++){
        update_k_1_mer(seq.substr(i,i+k-1),previous_id,id,i);
    }
}
void Index::add_unitig(std::string unitig,int id,int k){

        for(int i=0;i<unitig.length()-k+2;i++){
            std::string sub=unitig.substr(i,i+k-1);
            std::string rc=getCanonical(sub);
            std::string toStore=std::min(sub,rc);
            bool orientation=true;
            if(std::strcmp(sub.c_str(),rc.c_str())>0){
                orientation=false;
            }
            if(index_table.find(toStore)==index_table.end()){
                    std::vector<std::tuple<int,int,bool>> data;
                    data.push_back(std::tuple<int,int,bool>(id,i,orientation));
                    index_table[toStore]=data;
            }
            else{
                std::vector<std::tuple<int,int,bool>> data=index_table.at(toStore);
                data.push_back(std::tuple<int,int,bool>(id,i,orientation));
                index_table[toStore]=data;

            }
        }
}