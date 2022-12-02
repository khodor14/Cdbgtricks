#include "ParseGFA.h"
#include "index_kmer.h"
#include "CommonUtils.h"
#include "unitigmodification.h"
#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <algorithm>


void show_usage(){
    std::cerr<<"Usage:"<<"\n"
            <<"./main --input_graph <value> --input_genome <value> --k_mer_size <value>\n"
            <<"OR\n"
            <<"./main --input_graph <value> --k_mer_file <value> --k_mer_size <value>\n"
            <<"\n\n\n"
            <<"Options:"<<"\n"
            <<"\t"<<"-h[--help]"<<"\t prints this help message\n"
            <<"\t"<<"--input_graph"<<"\t the path to the pangenome graph in gfa or fasta format\n"
            <<"\t"<<"--input_genome"<<"\t the path to the genome used to augment the input graph\n"
            <<"\t"<<"--k_mer_size[-k]"<<" the size of the k-mer.\n\t\t\t It must be the same value used when constructing the input graph\n" 
            <<"\t"<<"--k_mer_file"<<"\t the file of absent k-mers from the graph if already computed\n"
            <<"\t"<<"--test\t\t "<<"used in testing mode to compare\n\t\t\t the augmented graph generated by this algo\n\t\t\t to an already augmented graph\n\t\t\t if this is set to true then the augmented graph must be given as input\n"
            <<"\t"<<"--augmented_graph "<<"the path to an already augmented graph in fasta or gfa\n\t\t\t it must be given if the argument --test is set to true\n"
            <<"\t"<<"--output_file_name[-o]"<<" the name of the output file\n" 
            ;
}
std::unordered_map<std::string,std::tuple<bool,std::string>> parseArgs(int argc,char **argv){
    std::unordered_map<std::string,std::tuple<bool,std::string>> arguments;
    arguments["graphfile"]=std::tuple<bool,std::string>(false,"");
    arguments["genomefile"]=std::tuple<bool,std::string>(false,"");
    arguments["kvalue"]=std::tuple<bool,std::string>(false,"");
    arguments["kmerfile"]=std::tuple<bool,std::string>(false,"");
    arguments["augmentedgraph"]=std::tuple<bool,std::string>(false,"");
    arguments["test"]=std::tuple<bool,std::string>(false,"");
    arguments["outputfilename"]=std::tuple<bool,std::string>(false,"");
    if(argc<2){
        show_usage();
        exit(0);
    }
    for(int i=1;i<argc;i++){
        if(!strcmp(argv[i],"--help") || !strcmp(argv[i],"-h")){
            show_usage();
        }
        else if(!strcmp(argv[i],"--test"))
        {
            arguments["test"]=std::tuple<bool,std::string>(true,"");
        }
        else if (!strcmp(argv[i],"-o")||!strcmp(argv[i],"--output_file_name"))
        {
            if(i+1<argc){
                arguments["outputfilename"]=std::tuple<bool,std::string>(true,"");
                i=i+1;
            }
            else{
                std::cerr<<"the output file name is missing, use --help or -h for details on how to use the program\n";
                show_usage();
                exit(0);
            }
        }
        else if (!strcmp(argv[i],"--augmented_graph"))
        {
            if(i+1<argc){
                arguments["augmentedgraph"]=std::tuple<bool,std::string>(true,argv[i+1]);
                i=i+1;
            }
            else{
                std::cerr<<"the path to an augmented graph is missing, use --help or -h for details on how to use the program\n";
                show_usage();
                exit(0);
            }
        }
        else if (!strcmp(argv[i],"--input_graph"))
        {
            if(i+1<argc){
                arguments["graphfile"]=std::tuple<bool,std::string>(true,argv[i+1]);
                i=i+1;
            }
            else{
                std::cerr<<"the path to a graph is missing, use --help or -h for details on how to use the program\n";
                show_usage();
                exit(0);
            }
        }
        else if(!strcmp(argv[i],"--input_genome"))
        {
            if(i+1<argc){
                arguments["genomefile"]=std::tuple<bool,std::string>(true,argv[i+1]);
                i=i+1;
            }
            else{
                std::cerr<<"the path to a genome is missing, use --help or -h for details on how to use the program\n";
                show_usage();
                exit(0);
            } 
        }
        else if(!strcmp(argv[i],"-k")||!strcmp(argv[i],"--k_mer_size"))
        {
            if(i+1<argc){
                arguments["kvalue"]=std::tuple<bool,std::string>(true,argv[i+1]);
                i=i+1;
            }
            else{
                std::cerr<<"the size of k-mer is missing, use --help or -h for details on how to use the program\n";
                show_usage();
                exit(0);
            } 
        }
        else if (!strcmp(argv[i],"--k_mer_file"))
        {
            if(i+1<argc){
                arguments["kmerfile"]=std::tuple<bool,std::string>(true,argv[i+1]);
                i=i+1;
            }
            else{
                std::cerr<<"the k-mers file is missing, use --help or -h for details on how to use the program\n";
                show_usage();
                exit(0);
            }
        }
        else{
            std::cerr<<"Unsupported argument "<<"\n";
            show_usage();
            exit(0);
        }
    }
    if(!(std::get<0>(arguments.at("kvalue"))&&std::get<0>(arguments.at("graphfile")))||!(std::get<0>(arguments.at("kmerfile"))||std::get<0>(arguments.at("genomefile")))){
        std::cerr<<"Some required arguments are missing\n";
        show_usage();
        exit(0);
    }
    if(std::get<0>(arguments.at("test"))&&!std::get<0>(arguments.at("augmentedgraph"))){
        std::cerr<<"When using testing mode, another graph is required\n";
        show_usage();
        exit(0);
    }
    return arguments;
}
std::vector<std::string> canonicalUnitigs(std::unordered_map<int,std::string> unitigs){
    std::vector<std::string> unitigs_vec;
    for(std::pair<int,std::string> elem:unitigs){
        unitigs_vec.push_back(getCanonical(elem.second));
    }
    std::sort(unitigs_vec.begin(),unitigs_vec.end());
    
    return unitigs_vec;
}

void validate_merging(const std::vector<std::string>& merged,const std::vector<std::string>& original_graph){
    /*
    two sorted vectors of unitigs
    the first one is the output of the algorithm we developed
    the second one is a graph that already constains the k-mers we added to the input graph of our algo
    */ 
    if(merged.size()!=original_graph.size()){
        std::cout<<"The augmented graph with my algo is not the same as the original graph (different number of unitigs)\n";
        std::cout<<merged.size()<<" "<<original_graph.size()<<"\n";
        return;
    }
    for(int i=0;i<merged.size();i++){
        if(strcmp(merged[i].c_str(),original_graph[i].c_str())){
            std::cout<<"The augmented graph with my algo is not the same as the original graph\n";
            return;
        }
    }
    std::cout<<"The augmented graph with my algo is the same as the original graph\n";

}
int main(int argc,char **argv){
    
    std::unordered_map<std::string,std::tuple<bool,std::string>> arguments=parseArgs(argc,argv);
    std::unordered_map<std::string,bool> k_mer;
    if(!std::get<0>(arguments.at("kmerfile"))){
        std::system("chmod +x utils.sh");
        std::system(("bash utils.sh "+std::get<1>(arguments.at("kvalue"))+" "+std::get<1>(arguments.at("graphfile"))+" "+std::get<1>(arguments.at("genomefile"))).c_str());
        k_mer=createHashTable("./kmers.txt");
    }
    else{
        k_mer=createHashTable(std::get<1>(arguments.at("kmerfile")));
    }

    GfaGraph g;
    GfaGraph g2=g.LoadFromFile(std::get<1>(arguments.at("graphfile")));
    Index ind=Index(1000,31);//stoi(std::get<1>(arguments.at("kvalue"))));
    ind.create(g2);
    
    std::unordered_map<int,std::string> constrct_unitigs=construct_unitigs_from_kmer(ind,k_mer);
    std::unordered_map<std::string,std::vector<std::tuple<int,int,bool>>> constrtc_index=index_constructed_unitigs(constrct_unitigs,ind.get_k());
    std::unordered_map<int,std::string> merged=g2.get_nodes();
    std::cout<<constrct_unitigs.size()<<" constructed unitigs\n";
    merge_unitigs(constrtc_index,ind,merged,constrct_unitigs);    
    if(std::get<0>(arguments.at("test"))){
        GfaGraph g3;
        GfaGraph to_compare=g3.LoadFromFile(std::get<1>(arguments.at("augmentedgraph")));
        validate_merging(canonicalUnitigs(merged),canonicalUnitigs(to_compare.get_nodes()));
    }

    if(std::get<0>(arguments.at("outputfilename"))){
        write_unitigs_to_fasta(merged,std::get<1>(arguments.at("outputfilename"))+".fa");
    }
    else{
        write_unitigs_to_fasta(merged,"augmented_by_our_algo_2.fa");
    }
}
