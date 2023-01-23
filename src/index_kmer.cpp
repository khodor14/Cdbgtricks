#include "index_kmer.h"
#include "CommonUtils.h"
#include <algorithm>
#include <cstring>
#include <vector>
#include <string>
#include <tuple>

Index::Index(int buckets,int k_size){
    k=k_size;
    number_of_buckets=buckets;
    index_table.rehash(number_of_buckets);
}
int Index::get_k(){
    return k;
}
void Index::create(GfaGraph& graph){
    std::unordered_map<int,std::string> nodes=graph.get_nodes();
   for (std::pair<int, std::string> node : nodes){
        int id=node.first;
        std::string unitig=node.second;
        int i=0;
        for(int i=0;i<=unitig.length()-k+1;++i){
            std::string sub=unitig.substr(i,k-1);
            std::string toStore=getCanonical(sub);
            bool orientation=isCanonical(sub);
            if(index_table.count(toStore)==0){
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
    if(index_table.count(getCanonical(kmer))==0){
        return std::vector<std::tuple<int,int,bool>>();
    }
    return index_table.at(getCanonical(kmer));
}
void Index::update_k_1_mer(std::string k_1_mer,int prev_id,int current_id,int position){
    /*
        the input are: a string of size (k-1)
        the previous id where this string occurs
        the current id where it occurs after the updte performed
        the new position where it occurs
    */
    std::vector<std::tuple<int,int,bool>> data=find(k_1_mer);//get the occurences of the string from the index

    //go over the occurences
    for (int i=0;i<data.size();i++){
        //which occurens has the previous id, the old occurence
        if(std::get<0>(data[i])==prev_id){
            data[i]=std::tuple<int,int,bool>(current_id,position,std::get<2>(data[i]));//update the values 
            break;//break since the (k-1)-mer occurs once in a unitig except for palindromic k-mers
        }
    }
    index_table[getCanonical(k_1_mer)]=data;
}
void Index::insert(std::string k_1_mer,int id,int position){
    /*
        the input is a string of size (k-1)
        the id of the unitig where this string occurs
        the position

        The k-1-mer gets inserted into the index
    */
   std::vector<std::tuple<int,int,bool>> occurences=find(k_1_mer);//get previous occurences

   bool orientation=isCanonical(k_1_mer);//find the orientation
   occurences.push_back(std::tuple<int,int,bool>(id,position,orientation));//insert the new occurence into vector
   index_table[getCanonical(k_1_mer)]=occurences;

}
void Index::insertSubUnitig(std::string unitig,int id,int starting_position,int ending_position){
    /*
        the inputs are:
                     a string unitig
                     its id
                     the starting position (int)
                     the ending position (int)
        Actions performed: all the (k-1)-mers occuring between starting_position and ending position gets inserted
    */
   for(int position=starting_position;position<=ending_position;position++){
        //call the function insert to insert the (k-1)-mer at position i
        insert(unitig.substr(position,k-1),id,position);
   }
}
void Index::update_unitig(std::string seq,int current_id,int previous_id,int starting_position,int ending_position){
    /*
        the inputs are:
                     a string unitig
                     its current id
                     its previous id
                     the starting position (int)
                     the ending position (int)
        Actions performed: the info all the (k-1)-mers occuring between starting_position and ending position get updated
    */
    for(int position=starting_position;position<=ending_position;position++){
        update_k_1_mer(seq.substr(position,k-1),previous_id,current_id,position);
    }
}