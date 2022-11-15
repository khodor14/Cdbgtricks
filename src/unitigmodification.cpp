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
std::vector<std::string> construct_unitigs_from_kmer(Index &unitig_index, std::unordered_map<std::string,bool> &kmers){
    std::vector<std::string> unitigs_set;
    for (std::pair<std::string, bool> element : kmers)
    {
        if(!element.second){//if this k-mer was not used in any unitig
            kmers[element.first]=true;
            std::string res=unitig_from_kmers(element.first,unitig_index,kmers);
            unitigs_set.push_back(res);//create unitig from this k-mer
        }
    }
    return unitigs_set;
}
