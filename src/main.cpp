#include "ParseGFA.h"
#include "index_kmer.h"
#include "CommonUtils.h"
#include "unitigmodification.h"
#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <chrono>
#include <sstream>
#include <sparsehash/sparse_hash_map>
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
            <<"\t"<<"--update_index[-u]"<<" index the constructed funitigs\n" 
            <<"\t"<<"-v"<<" verbosity\n" 
            <<"\t"<<"--output_kmers"<<" output kmers to temporary file\n"
            ;
}
google::sparse_hash_map<std::string,std::tuple<bool,std::string>> parseArgs(int argc,char **argv){
    google::sparse_hash_map<std::string,std::tuple<bool,std::string>> arguments;
    arguments["graphfile"]=std::tuple<bool,std::string>(false,"");
    arguments["genomefile"]=std::tuple<bool,std::string>(false,"");
    arguments["genomefile2"]=std::tuple<bool,std::string>(false,"");//remove me
    arguments["kvalue"]=std::tuple<bool,std::string>(false,"");
    arguments["kmerfile"]=std::tuple<bool,std::string>(false,"");
    arguments["augmentedgraph"]=std::tuple<bool,std::string>(false,"");
    arguments["test"]=std::tuple<bool,std::string>(false,"");
    arguments["outputfilename"]=std::tuple<bool,std::string>(false,"");
    arguments["verbosity"]=std::tuple<bool,std::string>(false,"");
    arguments["outputkmers"]=std::tuple<bool,std::string>(false,"");
    arguments["updateindex"]=std::tuple<bool,std::string>(false,"");
    if(argc<2){
        show_usage();
        exit(0);
    }
    for(int i=1;i<argc;i++){
        if(!strcmp(argv[i],"--help") || !strcmp(argv[i],"-h")){
            show_usage();
        }
        else if(!strcmp(argv[i],"-u")||!strcmp(argv[i],"--update_index"))
        {
            arguments["updateindex"]=std::tuple<bool,std::string>(true,"");
        }
        else if(!strcmp(argv[i],"--test"))
        {
            arguments["test"]=std::tuple<bool,std::string>(true,"");
        }
        else if(!strcmp(argv[i],"--output_kmers"))
        {
            if(i+1<argc){
                arguments["outputkmers"]=std::tuple<bool,std::string>(true,argv[i+1]);
                i=i+1;
            }
            else{
                std::cerr<<"the output kmers file name is missing, use --help or -h for details on how to use the program\n";
                show_usage();
                exit(0);
            }
        }
        else if(!strcmp(argv[i],"-v"))
        {
            arguments["verbosity"]=std::tuple<bool,std::string>(true,"");
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
            std::cerr<<"Unsupported argument "<<argv[i]<<"\n";
            show_usage();
            exit(0);
        }
    }
    if(!(std::get<0>(arguments["kvalue"])&&std::get<0>(arguments["graphfile"])&&std::get<0>(arguments["outputkmers"]))||!(std::get<0>(arguments["kmerfile"])||std::get<0>(arguments["genomefile"]))){
        std::cerr<<"Some required arguments are missing\n";
        show_usage();
        exit(0);
    }
    if(std::get<0>(arguments["test"])&&!std::get<0>(arguments["augmentedgraph"])){
        std::cerr<<"When using testing mode, another graph is required\n";
        show_usage();
        exit(0);
    }
    return arguments;
}
std::vector<std::string> canonicalUnitigs(google::sparse_hash_map<int,std::string> unitigs){
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
int total_kmers_in_unitigs(google::sparse_hash_map<int,std::string> unitigs,int k){
    int num_kmers_in_input=0;
    for(std::pair<int,std::string> unitig:unitigs){
        num_kmers_in_input=num_kmers_in_input+unitig.second.length()-k+1;
    }
    return num_kmers_in_input;
}
float average_unitig_length(google::sparse_hash_map<int,std::string> unitigs){
    //find average size of unitigs in the graph
    float average=0;
    for(std::pair<int,std::string> elem:unitigs){
        average=average+elem.second.length();
    }
    return average/unitigs.size();
}
int main(int argc,char **argv){
    //parse arguments
    google::sparse_hash_map<std::string,std::tuple<bool,std::string>> arguments=parseArgs(argc,argv);
    google::sparse_hash_map<std::string,bool> k_mer;
    //load the input graph
    GfaGraph g;
    GfaGraph g2=g.LoadFromFile(std::get<1>(arguments["graphfile"]));
    bool verbose=std::get<0>(arguments["verbosity"]);
    float time_kmtricks=0;
    //if we don't have a kmerfile, then a genome is passed
    if(!std::get<0>(arguments["kmerfile"])){
        std::string input_to_kmtricks=std::get<1>(arguments["graphfile"]);
        if(std::get<1>(arguments["graphfile"]).length()>4 && !std::strcmp(std::get<1>(arguments["graphfile"]).substr(std::get<1>(arguments["graphfile"]).length()-3).c_str(),"gfa")){
            //if the graph is not in fasta format, convert it to fasta as kmtricks accepts only fasta format (gzipped or not)
            g2.convertToFasta(std::get<1>(arguments["outputkmers"])+".fa");
            input_to_kmtricks=std::get<1>(arguments["outputkmers"])+".fa";
        }
        //call kmtricks
        std::system("chmod +x utils.sh");
        auto start=std::chrono::steady_clock::now();
        std::system(("bash utils.sh "+std::get<1>(arguments["kvalue"])+" "+input_to_kmtricks+" "+std::get<1>(arguments["genomefile"])+" "+std::get<1>(arguments["outputkmers"])+".txt").c_str());
        auto end =std::chrono::steady_clock::now();
        time_kmtricks=std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()*1e-9;
        //load the kmers
        k_mer=createHashTable(std::get<1>(arguments["outputkmers"])+".txt");
    }
    else{
        //load the kmers
        k_mer=createHashTable(std::get<1>(arguments["kmerfile"]));
    }
    //create the index of the graph
    Index ind=Index(1000,stoi(std::get<1>(arguments["kvalue"])));//
    auto start=std::chrono::steady_clock::now();
    ind.create(g2);
    auto end =std::chrono::steady_clock::now();
    float time_index = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()*1e-9;
    start=std::chrono::steady_clock::now();
    //construct the funitigs from the absent kmers
    google::sparse_hash_map<int,std::string> constrct_unitigs=construct_unitigs_from_kmer(ind,k_mer);   
    end=std::chrono::steady_clock::now();
    float time_construction=std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()*1e-9;
    /*
    debug mode to detect bugs
    */
    int num_kmer_absent=k_mer.size();
    int num_construct=constrct_unitigs.size();
    int original_unitigs=g2.get_nodes().size();
    std::cout<<"constructed "<<num_construct<<" from "<<num_kmer_absent<<" k-mers "<<std::endl;
    std::cout<<"# of unitigs in input is "<<original_unitigs<<std::endl;
    start=std::chrono::steady_clock::now();
    //index the funitigs (we store suffix->(id,position,orientation) and prefix->(id,position,orientation))
    google::sparse_hash_map<std::string,std::vector<std::tuple<int,int,bool>>> constrtc_index=index_constructed_unitigs(constrct_unitigs,ind.get_k());
    end=std::chrono::steady_clock::now();
    float time_index_constructed_unitigs=std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()*1e-9;
    google::sparse_hash_map<int,std::string> merged=g2.get_nodes();
    std::cout<<"kmers in unitigs="<<total_kmers_in_unitigs(merged,31)<<" kmers absent"<<k_mer.size()<<std::endl;
    std::cout<<"Originally we have: "<<merged.size()<<" unitigs\n";
    std::cout<<constrct_unitigs.size()<<" constructed unitigs\n";
    //statistics needed in testing
    int max_node_id=g2.get_max_node_id();
    int num_split=0;
    int num_join=0;
    float time_split=0;
    float time_join=0;
    float time_update=0;
    //merge the unitigs of the graph with the funitigs
    merge_unitigs(constrtc_index,ind,merged,constrct_unitigs,&max_node_id,&num_split,&num_join,&time_split,&time_join,&time_update,verbose,true);
    std::cout<<"After merging "<<total_kmers_in_unitigs(merged,31)<<std::endl;
    float average=average_unitig_length(g2.get_nodes());
    std::cout<<"Split\tJoin\t t index\t t kmtricks \t t construct\t t indexU \t t split \t t join\t number absent kmers \t average unitig length"<<std::endl;
    std::cout<<num_split<<"\t"<<num_join<<"\t"<<time_index<<"\t"<<time_kmtricks<<"\t"<<time_construction<<"\t"<<time_index_constructed_unitigs<<"\t"<<time_split<<"\t"<<time_join<<"\t"<<num_kmer_absent<<"\t"<<average<<std::endl;
    //if we are testing with an already augmented graph
    if(std::get<0>(arguments["test"])){
        GfaGraph g3;
        GfaGraph to_compare=g3.LoadFromFile(std::get<1>(arguments["augmentedgraph"]));
        validate_merging(canonicalUnitigs(merged),canonicalUnitigs(to_compare.get_nodes()));
    }
    //if the output file name prefix is passed as argument
    if(std::get<0>(arguments["outputfilename"])){
        write_unitigs_to_fasta(merged,std::get<1>(arguments["outputfilename"])+".fa");
    }
    //default autput file name
    else{
        write_unitigs_to_fasta(merged,"augmented_by_our_algo.fa");
    }
}
