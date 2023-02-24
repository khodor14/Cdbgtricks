#include "index_kmer.h"
#include "CommonUtils.h"
#include <algorithm>
#include <cstring>
#include <vector>
#include <string>
#include <tuple>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <sparsehash/sparse_hash_map>
#include <sparsehash/dense_hash_map>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/string.hpp>
#include <unordered_map>
namespace boost {
namespace serialization {
template<class Archive, typename... Args>
void serialize(Archive & ar, std::tuple<Args...> & t, const unsigned int version)
{
    ar & std::get<0>(t);
    ar & std::get<1>(t);
    ar & std::get<2>(t);
}
template<class Archive, typename Key, typename T>
void save(Archive& ar, const google::dense_hash_map<Key, T>& map,
          const unsigned int version)
{
    std::vector<std::pair<Key, T>> vec(map.begin(), map.end());
    ar & vec;
}

template<class Archive, typename Key, typename T>
void load(Archive& ar, google::dense_hash_map<Key, T>& map,
          const unsigned int version)
{
    std::vector<std::pair<Key, T>> vec;
    ar & vec;
    map.set_empty_key(Key()); // set the empty key
    //map.set_deleted_key(Key()); // set the deleted key
    for (const auto& pair : vec) {
        map[pair.first] = pair.second;
    }
}

template<class Archive, typename Key, typename T>
void serialize(Archive& ar, google::dense_hash_map<Key, T>& map,
               const unsigned int version)
{
    split_free(ar, map, version);
}
}
}
Index::Index(int k_size){
    k=k_size;
}
int Index::get_k(){
    return k;
}
void Index::create(GfaGraph& graph){
    std::unordered_map<int,std::string> nodes=graph.get_nodes();
    index_table.set_empty_key("");
   for (std::pair<int, std::string> node : nodes){
        int id=node.first;
        std::string unitig=node.second;
        int i=0;
        for(i=0;i<=unitig.length()-k+1;++i){
            std::string sub=unitig.substr(i,k-1);
            std::string toStore=getCanonical(sub);
            bool orientation=isCanonical(sub);
            if(index_table.count(toStore)==0){
                    std::vector<std::tuple<int,int,bool>> data;
                    data.push_back(std::tuple<int,int,bool>(id,i,orientation));
                    index_table[toStore]=data;
            }
            else{
                std::vector<std::tuple<int,int,bool>> data=index_table[toStore];
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
    return index_table[getCanonical(kmer)];
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
void Index::serialize(const std::string filename){
  std::ofstream ofs(filename, std::ios::binary);
  boost::archive::binary_oarchive oa(ofs);
  oa << index_table;
  ofs.close();
}
void Index::deserialize(const std::string filename){
  std::ifstream ifs(filename, std::ios::binary);
  boost::archive::binary_iarchive ia(ifs);
  ia >> index_table;
  ifs.close();
}
