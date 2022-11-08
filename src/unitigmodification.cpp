#include <unitigmodification.h>
#include <CommonUtils.h>
#include <strings.h>
#include<algorithm>
#include <vector>
std::vector<std::tuple<std::string,bool>> extend_right(std::string kmer,bool reverse, std::unordered_map<std::string,bool> kmers){
    std::vector<std::tuple<std::string,bool>> results;
    std::string alphabet="ACTG";
    for(int i=0;i<4;i++){
        if(kmers.find(kmer+alphabet[i])!=kmers.end()){//check if we add a character from right to (k-1)-mer is in the kmers
            results.push_back(std::tuple<std::string,bool>(kmer+alphabet[i],reverse));
        }
    } 
    return results;
}
std::vector<std::tuple<std::string,bool>> extend_left(std::string kmer,bool reverse, std::unordered_map<std::string,bool> kmers){
    std::vector<std::tuple<std::string,bool>> results;
    std::string alphabet="ACTG";
    for(int i=0;i<4;i++){
        if(kmers.find(alphabet[i]+kmer)!=kmers.end()){//check if we add a character from left to (k-1)-mer is in the kmers
            results.push_back(std::tuple<std::string,bool>(kmer+alphabet[i],reverse));
        }
    }
    return results;
}
std::vector<std::tuple<std::string,bool>> possible_extension_left(std::string kmer, std::unordered_map<std::string,bool> kmers){
    std::vector<std::tuple<std::string,bool>> results_left=extend_left(kmer.substr(1,kmer.length()),false,kmers);
    std::vector<std::tuple<std::string,bool>> results_right=extend_left(reverseComplement(kmer.substr(1,kmer.length())),true,kmers);

    for(auto current_element:results_right){ //add right to left
        results_left.push_back(current_element);
    }

    return results_left;
}
std::vector<std::tuple<std::string,bool>> possible_extension_right(std::string kmer, std::unordered_map<std::string,bool> kmers){
    std::vector<std::tuple<std::string,bool>> results_left=extend_left(reverseComplement(kmer.substr(1,kmer.length())),true,kmers);
    std::vector<std::tuple<std::string,bool>> results_right=extend_left(kmer.substr(1,kmer.length()),false,kmers);
    for(auto current_element:results_left){ //add left to right
        results_right.push_back(current_element);
    }

    return results_right;
}
std::string extend(std::string kmer,Index &unitig_index, std::unordered_map<std::string,bool> kmers){
    //we extend from right
    std::string result=kmer;
    while(true){//as long as we can extend
        std::vector<std::tuple<std::string,bool>> results_left=possible_extension_left(kmer,kmers);
        std::vector<std::tuple<std::string,bool>> result_right=possible_extension_right(kmer,kmers);

        //to extend left should be empty and right should have one element

        if(results_left.size()==0 && result_right.size()==1){
            std::tuple<std::string,bool> extension=result_right.front();

            std::string to_extend=std::get<0>(extension);//the kmer to be extended with our string input kmer
            bool reverse=std::get<1>(extension);
            kmers[to_extend]=true;// set it to true because we used it
            if(reverse){
                to_extend=reverseComplement(to_extend);//if we find it in reverse complement reverse it
            }
            result=result+to_extend[to_extend.length()-1];// add the last character to the result
            if(unitig_index.find(to_extend.substr(1,to_extend.length())).size()>0){//the suffix of the k-mer found is in the graph stop extension
                break;
            }
        }
        else{
            //either we have more than one extension from right or we have left extensions
            break;
        }
    }
    return result;
}
std::string unitig_from_kmers(std::string kmer,Index &unitig_index, std::unordered_map<std::string,bool> kmers){
    std::string right=extend(kmer,unitig_index,kmers);//extend the right of the kmer
    std::string left=reverseComplement(extend(reverseComplement(kmer),unitig_index,kmers));//extend the left, extend right of reverse then reverse the result

    return left+right.substr(kmer.length(),right.length());//return the concatination of left and right. Omit the first k characters of right because they are the suffix of left
}
std::vector<std::string> construct_unitigs_from_kmer(std::string kmer,Index &unitig_index, std::unordered_map<std::string,bool> kmers){
    std::vector<std::string> unitigs_set;
    for (std::pair<std::string, bool> element : kmers)
    {
        if(!element.second){//if this k-mer was not used in any unitig
            kmers[element.first]=true;
            unitigs_set.push_back(unitig_from_kmers(element.first,unitig_index,kmers));//create unitig from this k-mer
        }
    }
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