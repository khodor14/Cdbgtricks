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
   std::cout<<mer<<" "<<possible_character_for_extension.size()<<std::endl;
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
std::tuple<std::string,int,std::string,int> split_unitig(std::string unitig,int id,int new_id,int position,int k){
    std::string prefix=unitig.substr(0,unitig.length()-position+k-1);
    std::string suffix=unitig.substr(0,position);

    return std::tuple<std::string,int,std::string,int>(prefix,id,suffix,new_id);
}
std::tuple<std::string,int> join(std::string unitig1,std::string unitig2,int id,int k){
    std::string joined=unitig1.substr(0,unitig1.length()-k+2)+unitig2;
    return std::tuple<std::string,int>(joined,id);
}
std::string compact(std::string unitig,std::string kmer,bool compact_position){
    //left compact
    if(compact_position==0){
        return kmer[0]+unitig;
    }
    //right compact
    else{
        return unitig+kmer[kmer.length()-1];
    }
}
void modify_unitigs(std::string kmer,Index & table,std::unordered_map<int,std::string> unitigs){
    std::string prefix=kmer.substr(0,kmer.length()-1); //getting the (k-1) prefix of the kmer
    std::string suffix=kmer.substr(1,kmer.length());//getting the (k-1) suffix

    //retrieving the index of the each (k-1)-mer

    std::vector<std::tuple<int,int,bool>> prefix_vec=table.find(prefix);
    std::vector<std::tuple<int,int,bool>> suffix_vec=table.find(suffix);

    if((prefix_vec.size()==0 && suffix_vec.size()==0) || (prefix_vec.size()>1 && suffix_vec.size()>1)){ 
        //1-case the suffix and prefix of the k-mer don't exist, create a new unitig and index it
        //2-case 
            /*
                in this case we connect two disconnected component via this k-mer
                we don't need to spplit nor to join nor to compac
                    u1  \      / u5
                    u2 -- kmer -- u6
                    u3 --      -- u7
                    u4  /      \ u8
                */
        int id=unitigs.size()+1;
        unitigs[id]=kmer;
        table.add_unitig(kmer,id,kmer.length());
    }
    /*
    the (k-1)-mers of the k-mer occurs once in two different unitigs
    here we have three possibilies
            1-joining the unitigs if (k-1)-mers occur @ the exremeties of the unitigs. In this case the unitigs share (k-2) character overlap
            2- we must split the two unitigs if the (k-1)-mers occur at positions i and j respectively suhc that 0<i<|u1|-k+1 and 0<j<|u2|-k+1
            2- we must compact the first unitig and split the second one. In this case the suffix/prefix occur at the extremity of one and prefix/suffix occur in the middle of the second one
    */
    else if(prefix_vec.size()==1 && suffix_vec.size()==1){
        int id1=std::get<0>(prefix_vec.front());
        int pos1=std::get<1>(prefix_vec.front());
        int id2=std::get<0>(suffix_vec.front());
        int pos2=std::get<1>(suffix_vec.front());
        if((pos1==0 && pos2==unitigs[id2].length()-kmer.length()+1) ||  (pos2==0 && pos1==unitigs[id1].length()-kmer.length()+1)){
            std::tuple<std::string,int> joined=join(unitigs[id1],unitigs[id2],id1,kmer.length());

            unitigs[std::get<1>(joined)]=std::get<0>(joined);
            table.update_unitig(std::get<0>(joined),std::get<1>(joined),id1,kmer.length());
            table.update_unitig(std::get<0>(joined),std::get<1>(joined),id2,kmer.length());
            unitigs.erase(id1);
            unitigs.erase(id2);
        }
        else if(pos1==0 || pos1==unitigs[id1].length()-kmer.length()+1){
            //compact
            std::string compacted=compact(unitigs[id1],kmer,pos1);
            unitigs[id1]=compacted;
            
            //split
            std::tuple<std::string,int,std::string,int> splitted=split_unitig(unitigs[id2],id2,unitigs.size()+1,pos2,kmer.length());

            table.update_unitig(std::get<0>(splitted),std::get<1>(splitted),id2,kmer.size());
            table.update_unitig(std::get<2>(splitted),std::get<3>(splitted),id2,kmer.size());

        }
        else if (pos2==0 || pos2==unitigs[id2].length()-kmer.length()+1)
        {
            //compact
            std::string compacted=compact(unitigs[id2],kmer,pos2);
            unitigs[id2]=compacted;
            //split
            std::tuple<std::string,int,std::string,int> splitted=split_unitig(unitigs[id1],id1,unitigs.size()+1,pos1,kmer.length());

            table.update_unitig(std::get<0>(splitted),std::get<1>(splitted),id1,kmer.size());
            table.update_unitig(std::get<2>(splitted),std::get<3>(splitted),id1,kmer.size());
        }
        else{
            std::tuple<std::string,int,std::string,int> splitted1=split_unitig(unitigs[id1],id1,unitigs.size()+1,pos1,kmer.length());
            std::tuple<std::string,int,std::string,int> splitted2=split_unitig(unitigs[id1],id1,unitigs.size()+1,pos1,kmer.length());
            table.update_unitig(std::get<0>(splitted1),std::get<1>(splitted1),id1,kmer.size());
            table.update_unitig(std::get<2>(splitted1),std::get<3>(splitted1),id1,kmer.size());
            table.update_unitig(std::get<0>(splitted2),std::get<1>(splitted2),id2,kmer.size());
            table.update_unitig(std::get<2>(splitted2),std::get<3>(splitted2),id2,kmer.size());
        }
    }
    else if(prefix_vec.size()==1 && suffix_vec.size()>1){
        int id=std::get<0>(prefix_vec.front());
        int pos=std::get<1>(prefix_vec.front());
        std::tuple<std::string,int,std::string,int> splitted=split_unitig(unitigs[id],id,unitigs.size()+1,pos,kmer.length());
        table.update_unitig(std::get<0>(splitted),std::get<1>(splitted),id,kmer.size());
        table.update_unitig(std::get<2>(splitted),std::get<3>(splitted),id,kmer.size());
    }
    else if(prefix_vec.size()>1 && suffix_vec.size()==1){
        int id=std::get<0>(suffix_vec.front());
        int pos=std::get<1>(suffix_vec.front());
        std::tuple<std::string,int,std::string,int> splitted=split_unitig(unitigs[id],id,unitigs.size()+1,pos,kmer.length());
        table.update_unitig(std::get<0>(splitted),std::get<1>(splitted),id,kmer.size());
        table.update_unitig(std::get<2>(splitted),std::get<3>(splitted),id,kmer.size());
    }
    else if(prefix_vec.size()==0 && suffix_vec.size()>1){
        int id=unitigs.size()+1;
        unitigs[id]=kmer;
        table.add_unitig(kmer,id,kmer.length());
    }
    else if(prefix_vec.size()>1 && suffix_vec.size()==0){
        int id=unitigs.size()+1;
        unitigs[id]=kmer;
        table.add_unitig(kmer,id,kmer.length());
    }
    else if (prefix_vec.size()==1 && suffix_vec.size()==0)
    {
        int id=std::get<0>(suffix_vec.front());
        int pos=std::get<1>(suffix_vec.front());

        if(pos==0){
            int new_id=unitigs.size()+1;
            unitigs[new_id]=kmer;
            table.add_unitig(kmer,id,kmer.length());
        }
        else if(pos==unitigs[id].length()-kmer.length()+1){
            //compact
            std::string compacted=compact(unitigs[id],kmer,pos);
            table.add_unitig(compacted,id,kmer.length());

        }
        else{
            std::tuple<std::string,int,std::string,int> splitted=split_unitig(unitigs[id],id,unitigs.size()+1,pos,kmer.length());
            table.update_unitig(std::get<0>(splitted),std::get<1>(splitted),id,kmer.size());
            table.update_unitig(std::get<2>(splitted),std::get<3>(splitted),id,kmer.size());  
        }
    }
    else{
        int id=std::get<0>(prefix_vec.front());
        int pos=std::get<1>(prefix_vec.front());

        if(pos==0){
            int new_id=unitigs.size()+1;
            unitigs[new_id]=kmer;
            table.add_unitig(kmer,id,kmer.length());
        }
        else if(pos==unitigs[id].length()-kmer.length()+1){
            //compact
            std::string compacted=compact(unitigs[id],kmer,pos);
            table.add_unitig(compacted,id,kmer.length());

        }
        else{
            std::tuple<std::string,int,std::string,int> splitted=split_unitig(unitigs[id],id,unitigs.size()+1,pos,kmer.length());
            table.update_unitig(std::get<0>(splitted),std::get<1>(splitted),id,kmer.size());
            table.update_unitig(std::get<2>(splitted),std::get<3>(splitted),id,kmer.size());  
        }
    }
}