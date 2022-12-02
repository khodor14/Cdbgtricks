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
std::string unitig_from_kmers(std::string kmer,Index& unitig_index, std::unordered_map<std::string,bool>& kmers){
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
    int new_id=graph_unitigs.size()*2;//new id
    std::string original_unitig=graph_unitigs.at(id);//get the unitig to split
    std::string left=original_unitig.substr(0,position+k-1);//take the left from position 0 into position=position+k
    std::string right=original_unitig.substr(position,original_unitig.length()-position);//take the right from position till the end
    graph_unitigs[id]=left;//keep the id but change the unitig
    ind.insert(right.substr(0,k),new_id,0);//insert first (k-1)-mer into index
    ind.update_unitig(right,new_id,id,1,right.length()-k-1);//update the id and position of all other (k-1)-mers
    graph_unitigs[new_id]=right;//insert the right portion into the data struction
    std::cout<<"Unitig "<<id<<" is splitted at "<<position<<"\n";
    /*
        Note:this implementation needs to be oprimised
    */
}
int decision(const std::vector<std::tuple<int,int,bool>> &occ_graph,const std::vector<std::tuple<int,int,bool>> &occ_unitigs,const std::unordered_map<int,std::string> &unitigs,int k){
    int decision=0;//0 for nothing,1 for split, 2 for merge
    if(occ_graph.size()==0){
        return decision;
    }
    std::tuple<int,int,bool> first_occ=occ_graph.front();
    std::string unitig=unitigs.at(std::get<0>(first_occ));
    int pos=std::get<1>(first_occ);
    if((pos==0 || pos==unitig.length()-k+1) && (occ_graph.size()==1 && occ_unitigs.size()==1)){
        decision=2;//merge
    }
    else if(pos>0 && pos<unitig.length()-k+1){
        decision=1;//split
    }
    return decision;
}
void checkAndMerge(std::tuple<int,int,bool> occurence_graph,std::tuple<int,int,bool> occurence_unitig,std::unordered_map<std::string,std::vector<std::tuple<int,int,bool>>> &unitig_index,Index& index_graph,std::unordered_map<int,std::string> &graph_unitigs,std::unordered_map<int,std::string> &constructed_unitigs){
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
   std::string graph_u=graph_unitigs.at(id_graph);
   std::string constructed_unitig=constructed_unitigs.at(id_constructed);

   std::tuple<std::string,bool> concat=can_we_merge(position_unitig,position_graph,orientation_unitig,orientation_graph,constructed_unitig,graph_u,index_graph);//concatinate the unitigs if possible
   std::string merged=std::get<0>(concat);
   if(merged.length()>0){//the strings get concatinated
        std::cout<<"Unitig "<<id_graph<<" is merged with a constructed unitig\n";

        /*
            try to check if the second (k-1)-mer extremity of the constructed unitig is extremity of one unitig in the graph
            if yes check and merge this unitig with the constructed unitig then merge all three
            the result would be: unitig graph 1+ constructed unitig+unitig graph 2
                                ------------------- unitig graph 1
                                               ----------------------------- constructed unitig
                                                                        --------------------------------- unitig graph 2
        */
       int position2_u=abs(position_unitig-constructed_unitig.length()+index_graph.get_k()-1);//if its prefix then position1=0 and position2=|u|-k+1, if pos1=|u|-k+1 then pos2=0 so pos2=|pos1-|u|+k-1|
       std::vector<std::tuple<int,int,bool>> occ_g=index_graph.find(constructed_unitig.substr(position2_u,index_graph.get_k()-1));
       std::vector<std::tuple<int,int,bool>> occ_u;
       if(unitig_index.count(getCanonical(constructed_unitig.substr(position2_u,index_graph.get_k()-1)))>0){
            occ_u=unitig_index.at(getCanonical(constructed_unitig.substr(position2_u,index_graph.get_k()-1)));
       }
       int dec=decision(occ_g,occ_u,graph_unitigs,index_graph.get_k());
       if(dec==1){
        /*
        decision one means we need to split the unitig
        */

        split(graph_unitigs,index_graph,std::get<0>(occ_g.front()),std::get<1>(occ_g.front()),index_graph.get_k());
       }
       else if (dec==2){
            /*
            we may merge the unitig from the second extremity
            */
           int pos_u=std::get<1>(occ_u.front());
           int pos_g=std::get<1>(occ_g.front());
           bool orient_u=std::get<2>(occ_u.front());
           bool orient_g=std::get<2>(occ_g.front());
           std::tuple<std::string,bool> concat_other_extremity=can_we_merge(pos_u,pos_g,orient_u,orient_g,constructed_unitigs.at(std::get<0>(occ_u.front())),graph_unitigs.at(std::get<0>(occ_g.front())),index_graph);//check if we can merge
           if(std::get<0>(concat_other_extremity).length()>0){
            std::cout<<"Unitig "<<id_graph<<" is merged with a constructed unitig from the second extremity\n";
            if(std::get<1>(concat_other_extremity)){//the unitig of the graph is first
                /*
                graph unitig 2+constructed unitig+graph unitig 1
                */
               merged=graph_unitigs.at(std::get<0>(occ_g.front()))+merged.substr(index_graph.get_k()-1);
               index_graph.insertSubUnitig(merged,std::get<0>(occ_g.front()),graph_unitigs.at(std::get<0>(occ_g.front())).length()-index_graph.get_k()+2,merged.length()-graph_unitigs.at(id_graph).length()-1);//insert k-1 mers from constructed unitig to the index as 
               index_graph.update_unitig(merged,std::get<0>(occ_g.front()),id_graph,merged.length()-graph_unitigs.at(id_graph).length(),merged.length()-index_graph.get_k()-1);//update the id and position of the k-1 mer of the unitig used for merge
               graph_unitigs[std::get<0>(occ_g.front())]=merged;
               graph_unitigs.erase(id_graph);//remove the unitig used for merge
            }
            else{
                /*
                graph unitig 1+constructed unitig+graph unitig 2
                */
               merged=merged+graph_unitigs.at(std::get<0>(occ_g.front())).substr(index_graph.get_k()-1);
               index_graph.insertSubUnitig(merged,id_graph,graph_unitigs.at(id_graph).length()-index_graph.get_k()+2,merged.length()-graph_unitigs.at(std::get<0>(occ_g.front())).length()-1);//insert k-1 mers from constructed unitig to the index as 
               index_graph.update_unitig(merged,id_graph,std::get<0>(occ_g.front()),merged.length()-graph_unitigs.at(std::get<0>(occ_g.front())).length(),merged.length()-index_graph.get_k()-1);//update the id and position of the k-1 mer of the unitig used for merge
               graph_unitigs[std::get<0>(occ_g.front())]=merged;
               graph_unitigs.erase(id_graph);//remove the unitig used for merge
            }
            constructed_unitigs.erase(std::get<0>(occ_u.front()));
        }
       }
       if(dec==1 || dec==0){
        if(std::get<1>(concat)){
            index_graph.insertSubUnitig(merged,id_graph,graph_unitigs.at(id_graph).length()-index_graph.get_k()+2,merged.length()-index_graph.get_k()+1);//insert (k-1)-mers from constructed unitigs 
        }
        else{
            index_graph.insertSubUnitig(merged,id_graph,0,constructed_unitig.length()-index_graph.get_k());//insert k-1 mers from constructed unitigs
            index_graph.update_unitig(merged,id_graph,id_graph,constructed_unitig.length()-index_graph.get_k()+1,merged.length()-index_graph.get_k()+1);//shift positions
        }
        graph_unitigs[id_graph]=merged;

       }
       constructed_unitigs.erase(id_constructed);//remove this constructed unitig
       unitig_index.erase(getCanonical(constructed_unitig.substr(position2_u,index_graph.get_k()-1)));//remove the second extremity from the index
       
   }
}
std::tuple<std::string,bool> can_we_merge(int position_u,int position_g,bool orient_u,bool orient_g,std::string unitig_constrct,std::string unitig_graph,Index& index_table){
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
            the concatination if possible, otherwise empty string
            the order of the concatination: true if the unitig of the graph is first, false otherwise
    */
    std::string concatination="";
    bool order=true;//the first part is the unitig from the graph
    if(orient_g==orient_u){
        if(position_g==0 && position_u!=0){//result=unitig constructed+unitig of the graph
            //concatinate the unitig constructed to the one of the graph
            concatination=unitig_constrct+unitig_graph.substr(index_table.get_k()-1);
            order=false;
        }
        else if(position_g!=0 && position_u==0){//result=unitig of the graph+constructed unitig 
            concatination=unitig_graph.substr(0,unitig_graph.length()-index_table.get_k()+1)+unitig_constrct;
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
        concatination=reverseComplement(unitig_constrct)+unitig_graph.substr(index_table.get_k()-1);
        order=false;
       }
       else if(position_g!=0 && position_u!=0)
       {
            concatination=unitig_graph.substr(0,unitig_graph.length()-index_table.get_k()+1)+reverseComplement(unitig_constrct);
           //add to the unitig of the graph the reverse of the constructed unitig
       }
   }
    return std::tuple<std::string,bool>(concatination,order);
}
void merge_unitigs(std::unordered_map<std::string,std::vector<std::tuple<int,int,bool>>>& unitig_index,Index& graph_index,std::unordered_map<int,std::string>& graph_unitigs,std::unordered_map<int,std::string>& constructed_unitigs){
    /*
        this function takes the unitigs of the graph, the constructed unitigs and their index

        it merge these set of unitigs by calling appropriate functions    
    */
    for(std::pair<std::string,std::vector<std::tuple<int,int,bool>>> elem:unitig_index){
        std::string mer_unitig=elem.first;
        std::vector<std::tuple<int,int,bool>> graph_occ=graph_index.find(mer_unitig);//all occurrences in the graph
        int dec=decision(graph_occ,elem.second,graph_unitigs,graph_index.get_k());//find the decision of split or merge or nothing
        if(dec==1){//split
            split(graph_unitigs,graph_index,std::get<0>(graph_occ.front()),std::get<1>(graph_occ.front()),graph_index.get_k());
        }
        else if(dec==2){//merge
            checkAndMerge(graph_occ.front(),elem.second.front(),unitig_index,graph_index,graph_unitigs,constructed_unitigs);
        }    
   }

   //add the remaining constructed unitigs to the unitigs of the graph
   int id=100000;//graph_unitigs.size()+1;//id of the new added unitig
   for(std::pair<int,std::string> const_unitig:constructed_unitigs){
    graph_unitigs[id]=const_unitig.second;
    
    id++;
   }

}
std::unordered_map<std::string,std::vector<std::tuple<int,int,bool>>> index_constructed_unitigs(const std::unordered_map<int,std::string>& constructed_unitigs,int k){
    /*
    This function takes the constructed unitigs and build an index for them
    It returns an index:
        for each unitig u:
                    (k-1) suffix of u ->(id of u,0,orientation)
                    (k-1) prefix of u ->(id of u,|u|-k+1,orientation)
        each suffix or prefix is associated with a vector of occurences in the constructed unitigs
    */
  std::unordered_map<std::string,std::vector<std::tuple<int,int,bool>>> index;//empty index
  for(std::pair<int,std::string> elem:constructed_unitigs){//go over the constructed unitigs
    std::string prefix=elem.second.substr(0,k-1);//the prefix of unitig u=elem.second
    if(index.count(getCanonical(prefix))==0){//if it's not stored yet
        std::vector<std::tuple<int,int,bool>> v;//create an empty vector
        index[getCanonical(prefix)]=v;//associate v to this (k-1)-mer
    }
    std::string suffix=elem.second.substr(elem.second.length()-k+1);//suffix of u
    if(index.count(getCanonical(suffix))==0){//not stored yet
        std::vector<std::tuple<int,int,bool>> v;//empty vector
        index[getCanonical(suffix)]=v;//associate v to (k-1)-mer
    }
    std::vector<std::tuple<int,int,bool>> pref=index.at(getCanonical(prefix));//get the occurences of prefix
    pref.push_back(std::tuple<int,int,bool>(elem.first,0,isCanonical(prefix)));//add the new occurence to the vector
    index[getCanonical(prefix)]=pref;//store back the vector
    std::vector<std::tuple<int,int,bool>> suf=index.at(getCanonical(suffix));//get the occurences of suffix
    suf.push_back(std::tuple<int,int,bool>(elem.first,elem.second.length()-k+1,isCanonical(suffix)));//add the new occurence
    index[getCanonical(suffix)]=suf;//store back the occurences of suffix
  }
  return index;//return the constructed index
}