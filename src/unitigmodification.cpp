#include <unitigmodification.h>
#include <CommonUtils.h>
#include <strings.h>
#include<algorithm>

void search_modifications(std::string k_1_mer,int position,std::vector<uint64_t> hash_kmers,std::vector<std::tuple<int,char,std::tuple<int,bool,bool>>> possible_modifications){
    std::string alphabet("ACTG");
    std::string temp1=k_1_mer;
    std::string temp2=k_1_mer;
    std::string rev1;
    std::string rev2;
    int count(0);
    for(int i=0;i<alphabet.length();i++){
        temp1=k_1_mer;
        temp2=k_1_mer;
        temp2.append(0,alphabet[i]);
        temp1.push_back(alphabet[i]);
        rev1=reverseComplement(temp1);
        rev2=reverseComplement(temp2);
        uint64_t h1=hash(temp1);
        uint64_t h2=hash(temp2);
        uint64_t rev_h1=hash(rev1);
        uint64_t rev_h2=hash(rev2);
        auto it_hash=hash_kmers.begin();
        while (it_hash!=hash_kmers.end())
        {
            if(*it_hash==h1){
                possible_modifications.push_back(std::make_tuple(count,alphabet[i],std::make_tuple(position,true,true)));//forward-forward
                count++;
                it_hash=hash_kmers.erase(it_hash);
            }
            else if(*it_hash==h2){
                possible_modifications.push_back(std::make_tuple(count,alphabet[i],std::make_tuple(position,true,false)));//forward-backward
                count++;
                it_hash=hash_kmers.erase(it_hash);

            }
           else if(*it_hash==rev_h1){
                possible_modifications.push_back(std::make_tuple(count,alphabet[i],std::make_tuple(position,false,true)));//backward-forward
                count++;
                it_hash=hash_kmers.erase(it_hash);

            }
           else if(*it_hash==rev_h2){
                possible_modifications.push_back(std::make_tuple(count,alphabet[i],std::make_tuple(position,false,false)));//backward-backward
                count++;
                it_hash=hash_kmers.erase(it_hash);

            }
            else{
                it_hash++;
            }
        }

    }

}
void modify_unitig(std::string unitig,std::vector<uint64_t> hash_kmers,int k){
    std::vector<std::tuple<int,char,std::tuple<int,bool,bool>>> possible_modifications;
    std::string k_1_mer;
    for(int i=0;i<=unitig.length()-k;i++){
        k_1_mer=unitig.substr(i,i+k-1);
        search_modifications(k_1_mer,i,hash_kmers,possible_modifications);
        if(possible_modifications.size()!=0){
            //elongation from left is found
            if(i==0 && possible_modifications.size()==1 && std::get<1>(std::get<2>(*possible_modifications.begin()))&& !std::get<2>(std::get<2>(*possible_modifications.begin()))){
                unitig=unitig.append(0,std::get<1>(*possible_modifications.begin()));//need to correct also all the in-neighbor nodes
                i=0;
                possible_modifications.clear();
            }
            //elongation from left is found
            else if(i==unitig.length()-k && possible_modifications.size()==1 && std::get<1>(std::get<2>(*possible_modifications.begin()))&& std::get<2>(std::get<2>(*possible_modifications.begin()))){
                unitig.push_back(std::get<1>(*possible_modifications.begin()));//need to correct also all the out-neighbor nodes
                i=unitig.length()-k;
                possible_modifications.clear();
            }
            //need to create new node without touching the unitig
            else if(i==0 && possible_modifications.size()==1 && std::get<1>(std::get<2>(*possible_modifications.begin()))&& std::get<2>(std::get<2>(*possible_modifications.begin()))){
              //create a node and set its incoming edges the same as the unitig  
            }
            else if (i==unitig.length()-k && possible_modifications.size()==1 && std::get<1>(std::get<2>(*possible_modifications.begin()))&& !std::get<2>(std::get<2>(*possible_modifications.begin())))
            {
                //create a node and set its outgoing edges the same as the unitig
            }
            //need to split
            else{

            }
        }
    }
}
static GfaGraph modify_graph(GfaGraph g,std::vector<uint64_t> hash_kmers){
    std::cout<<"to be implemented";
}

