#include <unitigmodification.h>
#include <CommonUtils.h>
#include <strings.h>
#include<algorithm>
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