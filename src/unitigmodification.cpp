#include "unitigmodification.h"
#include "CommonUtils.h"
#include <strings.h>
#include<algorithm>
#include <vector>
#include <unordered_map>
#include <chrono>
#include "unitig.h"
std::unordered_map<uint64_t,std::vector<std::tuple<int,int,bool>>> index_constructed_unitigs(const std::unordered_map<int,Unitig>& constructed_unitigs,int k){
    /*
    This function takes the constructed unitigs and build an index for them
    It returns an index:
        for each unitig u:
                    (k-1) suffix of u ->(id of u,0,orientation)
                    (k-1) prefix of u ->(id of u,|u|-k+1,orientation)
        each suffix or prefix is associated with a vector of occurences in the constructed unitigs
    */
  std::unordered_map<uint64_t,std::vector<std::tuple<int,int,bool>>> index;//empty index
  for(std::pair<int,Unitig> elem:constructed_unitigs){//go over the constructed unitigs
    uint64_t prefix=elem.second.get_ith_mer(0,k-1);//the prefix of unitig u=elem.second
    std::tuple<uint64_t,bool> rev_canon_pref=reverseComplementCanonical(prefix,k-1);
    if(index.count(std::get<0>(rev_canon_pref))==0){//if it's not stored yet
        std::vector<std::tuple<int,int,bool>> v;//create an empty vector
        index[std::get<0>(rev_canon_pref)]=v;//associate v to this (k-1)-mer
    }
    uint64_t suffix=elem.second.get_ith_mer(elem.second.unitig_length()-k+1,k-1);//suffix of u
    std::tuple<uint64_t,bool> rev_canon_suff=reverseComplementCanonical(suffix,k-1);
    if(index.count(std::get<0>(rev_canon_suff))==0){//not stored yet
        std::vector<std::tuple<int,int,bool>> v;//empty vector
        index[std::get<0>(rev_canon_suff)]=v;//associate v to (k-1)-mer
    }
    std::vector<std::tuple<int,int,bool>> pref=index[std::get<0>(rev_canon_pref)];//get the occurences of prefix
    pref.push_back(std::tuple<int,int,bool>(elem.first,0,std::get<1>(rev_canon_pref)));//add the new occurence to the vector
    index[std::get<0>(rev_canon_pref)]=pref;//store back the vector
    std::vector<std::tuple<int,int,bool>> suf=index[std::get<0>(rev_canon_suff)];//get the occurences of suffix
    suf.push_back(std::tuple<int,int,bool>(elem.first,elem.second.unitig_length()-k+1,std::get<1>(rev_canon_suff)));//add the new occurence
    index[std::get<0>(rev_canon_suff)]=suf;//store back the occurences of suffix
  }
  return index;//return the constructed index
}
std::vector<char> possible_right_extension(uint64_t mer,std::unordered_map<uint64_t,bool> &kmers_to_add_to_graph){
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
        if(kmers_to_add_to_graph.find(canonical_bits((mer<<2)|baseToInt(alphabet[i])))!=kmers_to_add_to_graph.end()){
            //add the corresponding character to the vector
            possible_character_for_extension.push_back(alphabet[i]);
        }
   }
   //return the resultant vector
   return possible_character_for_extension;

}
std::string extend_right(std::string kmer,uint64_t mer, Index& unitigs_index,std::unordered_map<uint64,bool> &kmers_to_add_to_graph){
    /*
    the input is: a string k-mer representing the k-mer to be extended from right
                :the index of the unitigs, i.e. (k-1)-mer -> (unitig id,position,orientation)
                :an unordered map containing the canonical forms of k-mers as keys and boolean value as values
                                                                                        the boolean is used to say that the k-mer had been used to construct another unitig
    the outpus is a extended string
    */
   std::string extended_kmer=kmer;
   uint64_t suffix=mer>>2;//taking the suffix of the kmer to extend it
  while (true)//as long as we can extended the (k-1)-mer suffix of the k-mer
  {
    /*
    if the (k-1)-mer suffix of the k-mer is not in the graph
    we can try to extends
    */
   
   if(unitigs_index.find(hash(suffix)).size()==0){
    
    
    //we cannot extend if the suffix of the k-mer has left extension

    /*
    if the k-mer is ACTGA and we already have TCTGA in the k-mer set(to be added to the graph)
    then we cannot extend ACTGA from the right
    N.B:to test this we send the reverse complement of the suffix and we try to extend it from right
    */
    if(possible_right_extension(reverse_complement(suffix,unitigs_index.get_k()-1),kmers_to_add_to_graph).size()==1){
        //now we find the extensions of the actual suffix
        std::vector<char> possible_character_to_extend=possible_right_extension(suffix,kmers_to_add_to_graph);
        //to extend we should have only one character in the vector
        if(possible_character_to_extend.size()==1){
            char character_to_add=possible_character_to_extend.front();//get the character to augment it to our string result from the vector
            std::string key=suffix+character_to_add;
            kmers_to_add_to_graph[canonical_bits((suffix<<2)|baseToInt(character_to_add))]=true;
            extended_kmer=extended_kmer+character_to_add;//append it to the string only once, the append function takes as argument the number of times we need to append to a string
            suffix=(suffix<<2)|baseToInt(character_to_add);
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
std::string unitig_from_kmers_string(std::string kmer,uint64_t mer,Index& unitig_index, std::unordered_map<uint64_t,bool>& kmers){
    std::string right=extend_right(kmer,mer,unitig_index,kmers);//extend the right of the kmer
    std::string left=reverseComplement(extend_right(reverseComplement(kmer),reverse_complement(mer,unitig_index.get_k()),unitig_index,kmers));//extend the left, extend right of reverse then reverse the result
    return left+right.substr(kmer.length(),right.length());//return the concatenation of left and right. Omit the first k characters of right because they are the suffix of left
}
Unitig unitig_from_kmers(uint64_t kmer,Index &unitig_index, std::unordered_map<uint64_t,bool> &kmers){
    /*
    create a unitig from the bit encoding of the k-mer
    Input: uint64 representing the kmer, the index of graph and all kmers
    Return a unitig(8 bits per 4 bases,how many left unused bits(0 in this case) and how many right unused bits)
    */
    return Unitig(unitig_from_kmers_string(to_string(kmer,unitig_index.get_k()),kmer,unitig_index,kmers));
}
std::unordered_map<int,Unitig> construct_unitigs_from_kmer(Index &unitig_index, std::unordered_map<uint64_t,bool> &kmers,int k){
    /*
    simple modification is needed
    the input should be a map(int,bool) meaning that the hash values are already computed
    */
    std::unordered_map<int,Unitig> unitigs_map;//storing them in a map for fast access
    int i=1;
    for (std::pair<uint64_t, bool> element : kmers)
    {
        if(!element.second){//if this k-mer was not used in any unitig
            kmers[element.first]=true;
            Unitig res=unitig_from_kmers(element.first,unitig_index,kmers);//to be modified according to the map, the map shoul gives int->bool
            unitigs_map[i]=res;//create unitig from this k-mer
            i++;
        }
    }
    return unitigs_map;
}
int decision(const std::vector<std::tuple<int,int,bool>> &occ_graph,const std::vector<std::tuple<int,int,bool>> &occ_unitigs,const std::unordered_map<int,Unitig> &unitigs,int k){
    int decision=0;//0 for nothing,1 for split, 2 for merge
    if(occ_graph.size()==0){
        return decision;
    }
    std::tuple<int,int,bool> first_occ=occ_graph.front();
    Unitig unitig;
    if(unitigs.count(std::get<0>(first_occ))>0){
        unitig=unitigs.at(std::get<0>(first_occ));
    }
    else{
        std::cout<<"unitig not found "<<std::get<0>(first_occ)<<std::endl;
        exit(0);
    }
    int pos=std::get<1>(first_occ);
    if((pos==0 || pos==unitig.unitig_length()-k+1) && (occ_graph.size()==1 && occ_unitigs.size()==1)){
        decision=2;//merge
    }
    else if(pos>0 && pos<unitig.unitig_length()-k+1){
        decision=1;//split
    }
    return decision;
}
void split(std::unordered_map<int,Unitig> &graph_unitigs,Index &ind,int id,int position,int k,int *max_node_id,int *num_split,float *time_split,bool verbose){
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
    auto start=std::chrono::steady_clock::now();
    int new_id=*max_node_id+1;//new id to modify, not a good solution, should modify the merge function before
    Unitig original_unitig=graph_unitigs[id];//.at(id);
    *num_split=*num_split+1;
    Unitig left=Unitig(0,(8-2*(position+ind.get_k()-1)%4),std::vector<uint8_t>(original_unitig.get_encoding().begin(),original_unitig.get_encoding().begin()+(position+ind.get_k()-1)/2));
    //original_unitig.substr(0,position+k-1);//take the left from position 0 into position=position+k
    Unitig right=Unitig(2*(position%4),original_unitig.get_right_unused(),std::vector<uint8_t>(original_unitig.get_encoding().begin()+position/4,original_unitig.get_encoding().end()));//take the right from position till the end
    
    graph_unitigs[id]=left;//keep the id but change the unitig
    /*
    To update the index properly
    */
    //ind.insert(hash(right.substr(0,k-1)),new_id,0);//insert first (k-1)-mer into index
    //ind.update_unitig(right,new_id,id,1,right.length()-k+1);//update the id and position of all other (k-1)-mers
    ind.insert(left.get_ith_mer(0,ind.get_k()-1),new_id,0);//insert first (k-1)-mer into index
    ind.update_unitig(right,new_id,id,1,right.length()-k+1,true);//update the id and position of all other (k-1)-mers
    graph_unitigs[new_id]=right;//insert the right portion into the data struction
    auto end=std::chrono::steady_clock::now();
    *time_split=*time_split+std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()*1e-9;
    if(verbose){
        std::cout<<"Unitig "<<id<<" is splitted at "<<position<<"\n";
    }
    *max_node_id=*max_node_id+1;//incrementing the max nod id
    /*
        Note:this implementation needs to be oprimised
    */
}
void checkAndMerge(std::tuple<int,int,bool> occurence_graph,std::tuple<int,int,bool> occurence_unitig,std::unordered_map<uint64_t,std::vector<std::tuple<int,int,bool>>> &unitig_index,Index& index_graph,std::unordered_map<int,Unitig> &graph_unitigs,std::unordered_map<int,Unitig> &constructed_unitigs,int *max_node_id,int *num_split,int *num_join,float *time_split,float * time_join,bool verbose){
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
   Unitig graph_u=graph_unitigs[id_graph];//make sure this does not throw exception
   Unitig constructed_unitig=constructed_unitigs[id_constructed];
   std::tuple<std::string,bool> concat=can_we_merge(position_unitig,position_graph,orientation_unitig,orientation_graph,constructed_unitig,graph_u,index_graph);//concatinate the unitigs if possible
   std::string merged=std::get<0>(concat);
   if(merged.length()>0){//the strings get concatenated
        *num_join=*num_join+1;
        if(verbose){
            std::cout<<"Unitig "<<id_graph<<" is merged with a constructed unitig\n";
        }
        /*
            try to check if the second (k-1)-mer extremity of the constructed unitig is extremity of one unitig in the graph
            if yes check and merge this unitig with the constructed unitig then merge all three
            the result would be: unitig graph 1+ constructed unitig+unitig graph 2
                                ------------------- unitig graph 1
                                               ----------------------------- constructed unitig
                                                                        --------------------------------- unitig graph 2
        */
       int position2_u=abs(int(position_unitig-constructed_unitig.unitig_length()+index_graph.get_k()-1));//if its prefix then position1=0 and position2=|u|-k+1, if pos1=|u|-k+1 then pos2=0 so pos2=|pos1-|u|+k-1|
       std::vector<std::tuple<int,int,bool>> occ_g=index_graph.find(constructed_unitig.get_ith_mer(position2_u,index_graph.get_k()-1));
       std::vector<std::tuple<int,int,bool>> occ_u;
       std::tuple<uint64_t,bool> data=reverseComplementCanonical(constructed_unitig.get_ith_mer(position2_u,index_graph.get_k()-1));
       if(unitig_index.count(std::get<0>(data))>0){
            occ_u=unitig_index[std::get<0>(data)];
       }
       int dec=decision(occ_g,occ_u,graph_unitigs,index_graph.get_k());
       if(dec==1){
        /*
        decision one means we need to split the unitig
        */

        split(graph_unitigs,index_graph,std::get<0>(occ_g.front()),std::get<1>(occ_g.front()),index_graph.get_k(),max_node_id,num_split,time_split,verbose);
       }
       else if (dec==2){
            /*
            we may merge the unitig from the second extremity
            */
           int pos_u=std::get<1>(occ_u.front());
           int pos_g=std::get<1>(occ_g.front());
           bool orient_u=std::get<2>(occ_u.front());
           bool orient_g=std::get<2>(occ_g.front());
           std::tuple<std::string,bool> concat_other_extremity=can_we_merge(pos_u,pos_g,orient_u,orient_g,constructed_unitigs[std::get<0>(occ_u.front())],graph_unitigs[std::get<0>(occ_g.front())],index_graph);//check if we can merge
           if(std::get<0>(concat_other_extremity).length()>0){
            /*
            we merged the unitig from the second end-point
            */
           *num_join=*num_join+1;
           if(verbose){
            std::cout<<"Unitig "<<id_graph<<" is merged with a constructed unitig from the second extremity\n";
           }
            auto start=std::chrono::steady_clock::now();
            if(std::get<1>(concat_other_extremity)){//the unitig of the graph is first
                /*
                graph unitig 2+constructed unitig+graph unitig 1
                */
               merged=std::get<0>(concat_other_extremity)+merged.substr(constructed_unitig.length());
               //update the index of the first unitig, as the orientation of (k-1)-mers may change due to the reverse complement call in the merge function
               Unitig merged_u=Unitig(merged);
               index_graph.update_unitig(merged_u,std::get<0>(occ_g.front()),std::get<0>(occ_g.front()),0,graph_unitigs[std::get<0>(occ_g.front())].unitig_length()-index_graph.get_k()+1,false);
               //insert the (k-1)-mers from the funitig (future unitig)
               index_graph.insertSubUnitig(merged_u,std::get<0>(occ_g.front()),graph_unitigs[std::get<0>(occ_g.front())].unitig_length()-index_graph.get_k()+2,merged.length()-graph_unitigs[id_graph].unitig_length()-1);//insert k-1 mers from constructed unitig to the index as 
               //update the id as well as the positions of the (k-1)-mer of the second unitig
               index_graph.update_unitig(merged_u,std::get<0>(occ_g.front()),id_graph,merged.length()-graph_unitigs[id_graph].unitig_length(),merged.length()-index_graph.get_k()+1,false);//update the id and position of the k-1 mer of the unitig used for merge
               //associate the merged unitig (munitig) to the id of the left unitig
               graph_unitigs[std::get<0>(occ_g.front())]=merged_u;
               graph_unitigs.erase(id_graph);//remove the unitig used for merge
            }
            else{
                /*
                graph unitig 1+constructed unitig+graph unitig 2
                */
               merged=merged+std::get<0>(concat_other_extremity).substr(constructed_unitig.length());
               //update the index of the first unitig, as the orientation of (k-1)-mers may change due to the reverse complement call in the merge function
               Unitig merged_u=Unitig(merged);
               index_graph.update_unitig(merged_u,id_graph,id_graph,0,graph_unitigs[id_graph].unitig_length()-index_graph.get_k()+1);
                //insert the (k-1)-mers from the funitig (future unitig)
               index_graph.insertSubUnitig(merged_u,id_graph,graph_unitigs[id_graph].unitig_length()-index_graph.get_k()+2,merged.length()-graph_unitigs[std::get<0>(occ_g.front())].unitig_length()-1);//insert k-1 mers from constructed unitig to the index as 
               //update the id as well as the positions of the (k-1)-mer of the second unitig
               index_graph.update_unitig(merged,id_graph,std::get<0>(occ_g.front()),merged.length()-graph_unitigs[std::get<0>(occ_g.front())].unitig_length(),merged.length()-index_graph.get_k()+1);//update the id and position of the k-1 mer of the unitig used for merge
               //associate the merged unitig (munitig) to the id of the left unitig
               graph_unitigs[id_graph]=merged_u;
               graph_unitigs.erase(std::get<0>(occ_g.front()));//remove the unitig used for merge
            }
            auto end=std::chrono::steady_clock::now();
            *time_join=*time_join+std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()*1e-9;
            constructed_unitigs.erase(id_constructed);
        }
       }
       if(dec==1 || dec==0){
        auto start=std::chrono::steady_clock::now();
        if(std::get<1>(concat)){
            //update the orientation
            Unitig merged_u=Unitig(merged);
            index_graph.update_unitig(merged_u,id_graph,id_graph,0,graph_unitigs[id_graph].unitig_length()-index_graph.get_k()+1);
            //insert the (k-1)-mers from funitig
            index_graph.insertSubUnitig(merged,id_graph,graph_unitigs[id_graph].unitig_length()-index_graph.get_k()+2,merged.length()-index_graph.get_k()+1);//insert (k-1)-mers from constructed unitigs 
        }
        else{
            //insert the (k-1)-mers from funitig
            Unitig merged_u=Unitig(merged);
            index_graph.insertSubUnitig(merged_u,id_graph,0,constructed_unitig.unitig_length()-index_graph.get_k());//insert k-1 mers from constructed unitigs
            //update the position,the orientation of (k-1)-mers from the unitig used for merge
            index_graph.update_unitig(merged,id_graph,id_graph,constructed_unitig.unitig_length()-index_graph.get_k()+1,merged.length()-index_graph.get_k()+1);//shift positions
        }
        graph_unitigs[id_graph]=Unitig(merged);//to be changed by adding default and assignment operator =
        constructed_unitigs.erase(id_constructed);//remove this funitig
        auto end=std::chrono::steady_clock::now();
        *time_join=*time_join+std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()*1e-9;
       }
       unitig_index.erase(std::get<0>(reverseComplementCanonical(constructed_unitig.get_ith_mer(position2_u,index_graph.get_k()-1))));//remove the second extremity from the index
       
   }
}
std::tuple<std::string,bool> can_we_merge(int position_u,int position_g,bool orient_u,bool orient_g,Unitig unitig_constrct,Unitig unitig_graph,Index& index_table){
    /*
    the input are:
                position in the constructed unitig position_u
                position in the unitig of the graph position_g
                orientation in the constructed unitig
                orientation in the unitig of the graph
                the unitig of the graph
                the unitig constructed
                the size of the k-mer
    
    return:
            the concatenation if possible, otherwise empty string
            the order of the concatenation: true if the unitig of the graph is first, false otherwise
    */
    std::string concatenation="";
    bool order=true;//the first part is the unitig from the graph
    if(orient_g==orient_u){
        if(position_g==0 && position_u!=0){//result=unitig constructed+unitig of the graph
            //concatenate the unitig constructed to the one of the graph
            concatenation=unitig_constrct.to_string()+unitig_graph.to_string().substr(index_table.get_k()-1);
            order=false;
        }
        else if(position_g!=0 && position_u==0){//result=unitig of the graph+constructed unitig 
            concatenation=unitig_graph.to_string()+unitig_constrct.to_string().substr(index_table.get_k()-1);
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
       if(position_g==0 && position_u==0){
        //reverse the constructed unitig then add to it the unitig of the graph
        concatenation=reverseComplement(unitig_graph.to_string())+unitig_constrct.to_string().substr(index_table.get_k()-1);
       }
       else if(position_g!=0 && position_u!=0)
       {
            concatenation=unitig_constrct.to_string()+reverseComplement(unitig_graph.to_string()).substr(index_table.get_k()-1);
            order=false;
           //add to the unitig of the graph the reverse of the constructed unitig
       }
   }
    return std::tuple<std::string,bool>(concatenation,order);
}
void merge_unitigs(std::unordered_map<uint64_t,std::vector<std::tuple<int,int,bool>>>& unitig_index,Index& graph_index,std::unordered_map<int,Unitig>& graph_unitigs,std::unordered_map<int,Unitig>& constructed_unitigs,int *max_node_id,int *num_split,int *num_join,float *time_split,float * time_join,float * time_update,bool verbose,bool update_index){
    /*
        this function takes the unitigs of the graph, the constructed unitigs and their index

        it merge these set of unitigs by calling appropriate functions    
    */
    for(std::pair<uint64_t,std::vector<std::tuple<int,int,bool>>> elem:unitig_index){
        uint64_t mer_unitig=elem.first;
        std::vector<std::tuple<int,int,bool>> graph_occ=graph_index.find(mer_unitig);//all occurrences in the graph
        int dec=decision(graph_occ,elem.second,graph_unitigs,graph_index.get_k());//find the decision of split or merge or nothing
        if(dec==1){//split
            split(graph_unitigs,graph_index,std::get<0>(graph_occ.front()),std::get<1>(graph_occ.front()),graph_index.get_k(),max_node_id,num_split,time_split,verbose);
        }
        else if(dec==2){//merge
            checkAndMerge(graph_occ.front(),elem.second.front(),unitig_index,graph_index,graph_unitigs,constructed_unitigs,max_node_id,num_split,num_join,time_split,time_join,verbose);
        }    
   }
   //add the remaining constructed unitigs to the unitigs of the graph
    auto start=std::chrono::steady_clock::now();
   for(std::pair<int,Unitig> const_unitig:constructed_unitigs){
    if(const_unitig.second.unitig_length()>0){
        *max_node_id=*max_node_id+1;
        graph_unitigs[*max_node_id]=const_unitig.second;
        if(update_index){
            graph_index.insertSubUnitig(const_unitig.second,*max_node_id,0,const_unitig.second.length()-graph_index.get_k()+1);
        }
    }
   }
    auto end=std::chrono::steady_clock::now();
    *time_update=std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()*1e-9;
}