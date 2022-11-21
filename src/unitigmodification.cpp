#include "unitigmodification.h"
#include "CommonUtils.h"
#include <strings.h>
#include<algorithm>
#include <vector>
#include <unordered_map>

std::vector<char> possible_right_extension(std::string mer,std::unordered_map<std::string,bool> &kmers_to_add_to_graph){
    /*
    the input is :a string (k-1)-mer to be extended from the right
                 :an unordered map containing the canonical forms of k-mers as keys and boolean value as values
                                                                                        the boolean is used to say that the k-mer had been used to construct another unitig
    the output is a vector containing the possible characters that can extend the (k-1)-mer to a k-mer 

    example: if the input is mer=ACTG and the map is {ACTGA:false,ACTGT:false,TTGCT:false}
    then the ouput should be <A,T>
    */

   std::vector<char> possible_character_for_extension;
   std::string alphabet="ACTG";
   for(int i=0;i<alphabet.length();i++){
        //the canonical of mer+character from the alphabet is present
        if(kmers_to_add_to_graph.find(getCanonical(mer+alphabet[i]))!=kmers_to_add_to_graph.end()){
            //add the corresponding character to the vector
            possible_character_for_extension.push_back(alphabet[i]);
        }
   }
   //return the resultant vector
   return possible_character_for_extension;

}
std::string extend_right(std::string kmer, Index& unitigs_index,std::unordered_map<std::string,bool> &kmers_to_add_to_graph){
    /*
    the input is: a string k-mer representing the k-mer to be extended from right
                :the index of the unitigs, i.e. (k-1)-mer -> (unitig id,position,orientation)
                :an unordered map containing the canonical forms of k-mers as keys and boolean value as values
                                                                                        the boolean is used to say that the k-mer had been used to construct another unitig
    the outpus is a extended string
    */
   int k_mer_length=kmer.length();
   std::string extended_kmer=kmer;
   std::string suffix=kmer.substr(1,kmer.length());//taking the suffix of the kmer to extend it
  while (true)//as long as we can extended the (k-1)-mer suffix of the k-mer
  {
    /*
    if the (k-1)-mer suffix of the k-mer is not in the graph
    we can try to extends
    */
   
   if(unitigs_index.find(suffix).size()==0){
    
    
    //we cannot extend if the suffix of the k-mer has left extension

    /*
    if the k-mer is ACTGA and we already have TCTGA in the k-mer set(to be added to the graph)
    then we cannot extend ACTGA from the right
    N.B:to test this we send the reverse complement of the suffix and we try to extend it from right
    */
    if(possible_right_extension(reverseComplement(suffix),kmers_to_add_to_graph).size()==1){
        //now we find the extensions of the actual suffix
        std::vector<char> possible_character_to_extend=possible_right_extension(suffix,kmers_to_add_to_graph);
        //to extend we should have only one character in the vector
        if(possible_character_to_extend.size()==1){
            char character_to_add=possible_character_to_extend.front();//get the character to augment it to our string result from the vector
            std::string key=suffix+character_to_add;
            kmers_to_add_to_graph[getCanonical(key)]=true;
            extended_kmer=extended_kmer+character_to_add;//append it to the string only once, the append function takes as argument the number of times we need to append to a string
            suffix=suffix.substr(1,suffix.size())+character_to_add;
            //suffix=extended_kmer.substr(extended_kmer.length()-k_mer_length-1,extended_kmer.length());//takes the (k-1)-mer suffix of the extended k-mer
        }
        else{
            //we have zero or more than 1 possible right extension of the (k-1)-mer
            break;
        }
    }
    else{
        //the (k-1)-mer suffix of the string has more than one left parent extensions
        break;
    }
   }
   else{
    //the suffix already exists in the graph
    break;
   }
  }
   //return the resultant k-mer
    return extended_kmer;
}
std::string unitig_from_kmers(std::string kmer,Index &unitig_index, std::unordered_map<std::string,bool> &kmers){
    std::string right=extend_right(kmer,unitig_index,kmers);//extend the right of the kmer
    std::string left=reverseComplement(extend_right(reverseComplement(kmer),unitig_index,kmers));//extend the left, extend right of reverse then reverse the result
    return left+right.substr(kmer.length(),right.length());//return the concatination of left and right. Omit the first k characters of right because they are the suffix of left
}
std::unordered_map<int,std::string> construct_unitigs_from_kmer(Index &unitig_index, std::unordered_map<std::string,bool> &kmers){
    std::unordered_map<int,std::string> unitigs_map;//storing them in a map for fast access
    int i=1;
    for (std::pair<std::string, bool> element : kmers)
    {
        if(!element.second){//if this k-mer was not used in any unitig
            kmers[element.first]=true;
            std::string res=unitig_from_kmers(element.first,unitig_index,kmers);
            unitigs_map[i]=res;//create unitig from this k-mer
            i++;
        }
    }
    return unitigs_map;
}
void split(std::unordered_map<int,std::string> &graph_unitigs,Index &ind,int id,int position,int k){
    /*
    The input are:
                the unitigs of the graph
                the index
                the id of unitig to split
                the position
                the size of (k-1)-mer
    
    Actions performed: split the unitig whose id is the parameter id at the given position into left and right
                        create a new unitig whose id is the number of unitigs+1 and its sequence is the right portion
                        insert the first (k-1)-mer of right in the index
                        update the id and positions of all (k-1)-mers of right starting at 1
                        insert the new unitig into the unitig sets
    */
    int new_id=graph_unitigs.size();//new id
    std::string original_unitig=graph_unitigs.at(id);//get the unitig to split
    std::string left=original_unitig.substr(0,position+k-1);//take the left from position 0 into position=position+k
    std::string right=original_unitig.substr(position,original_unitig.length()-position);//take the right from position till the end
    graph_unitigs[id]=left;//keep the id but change the unitig
    ind.insert(right.substr(0,k),new_id,0);//insert first (k-1)-mer into index
    ind.update_unitig(right,new_id,id,1,right.length()-k-1);//update the id and position of all other (k-1)-mers
    graph_unitigs[new_id]=right;//insert the right portion into the data struction
    /*
        Note:this implementation needs to be oprimised
    */
}
void checkAndMerge(std::tuple<int,int,bool> occurence_graph,std::tuple<int,int,bool> occurence_unitig,Index& index_graph,std::unordered_map<int,std::string> &graph_unitigs,std::unordered_map<int,std::string> &constructed_unitigs){
    /*
        the inputs of the function:
                        the occurence of k-1 mer in the graph
                        the occurence of the same k-1 mer in the constructed unitigs
                        the unitigs of the graph
                        the unitigs constructed
        actions performed in this function: merge two unitig if we can merge them

    */
   int position_graph=std::get<1>(occurence_graph);
   int position_unitig=std::get<1>(occurence_unitig);
   bool orientation_graph=std::get<2>(occurence_graph);
   bool orientation_unitig=std::get<2>(occurence_unitig);
   int id_graph=std::get<0>(occurence_graph);
   int id_constructed=std::get<0>(occurence_unitig);
   std::string concatination="";
   if(orientation_graph==orientation_unitig){
        if(position_graph==0 && position_unitig!=0){//result=unitig constructed+unitig of the graph
            concatination=constructed_unitigs.at(id_constructed)+graph_unitigs.at(id_graph).substr(index_graph.get_k()-1,graph_unitigs.at(id_graph).length()-index_graph.get_k()+1);//concatination
            index_graph.insertSubUnitig(concatination,id_graph,0,concatination.length()-graph_unitigs.at(id_graph).length());//insert first (k-1)-mers from the constructed unitig
            index_graph.update_unitig(concatination,id_graph,id_graph,concatination.length()-graph_unitigs.at(id_graph).length(),concatination.length());//update the positions of the (k-1)-mers from the unitig of the graph
        }
        else if(position_graph!=0 && position_unitig==0){//result=unitig of the graph+constructed unitig 
            concatination=graph_unitigs.at(id_graph)+constructed_unitigs.at(id_constructed).substr(index_graph.get_k()-1,constructed_unitigs.at(id_constructed).length()-index_graph.get_k()+1);
            index_graph.insertSubUnitig(concatination,id_graph,concatination.length()-graph_unitigs.at(id_graph).length(),concatination.length());//insert (k-1)-mers from the constructed unitig
        }
        /*
        otherwise don't merge the unitigs
        */
   }
   else{
        /*
        either pos u=pos graph =0
        or pos u#0 and pos graph#0

        in these cases we merge
        */
   }
}