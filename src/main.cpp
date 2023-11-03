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
#include <cmath>
#include "unitig.h"
#include "index.hpp"
#include "../external/pthash/external/essentials/include/essentials.hpp"
void show_usage(){
    std::cerr<<"Usage: ./ccdbgupdater [COMMAND] [PARAMETERS]"<<"\n\n"
            <<"./ccdbgupdater --input_graph <value> --input_genome <value> --k_mer_size <value>\n"
            <<"OR\n"
            <<"./ccdbgupdater --input_graph <value> --k_mer_file <value> --k_mer_size <value>\n"
            <<"\n\n\n"
            
            <<"[COMMAND]:\n"
            <<"\tupdate \t\t update a compacted de Bruijn graph by a new genome\n"
            <<"\tindex \t\t index a compacted de Bruijn graph by (k-1)-mers\n"
            <<"\tconvert \t\t convert a graph from GFA/FASTA/binary to binary/Fasta\n\n"
            
            <<"[PARAMETERS]: update\n\n"
            <<"\t"<<"-h[--help]"<<"\t prints this help message\n"
            <<"\t"<<"--input_graph[-ig]"<<"\t the path to the pangenome graph in gfa or fasta format\n"
            <<"\t"<<"--input_genome"<<"\t the path to the genome used to augment the input graph\n"
            <<"\t"<<"--k_mer_size[-k]"<<" the size of the k-mer.\n\t\t\t It must be the same value used when constructing the input graph\n"
            <<"\t"<<"--minimizer_size[-m]"<<" the size of the minimizer (m<k)\n"
            <<"\t"<<"--k_mer_file"<<"\t the file of absent k-mers from the graph if already computed\n"
            <<"\t"<<"--smallest_merge[-s]"<<"\t the threshold for merging buckets (all buckets below this threshold are merged)\n"
            <<"\t"<<"--output_file_name[-o]"<<" the name of the output file\n" 
            <<"\t"<<"--update_index[-u]"<<" index the constructed funitigs\n"
            <<"\t"<<"--load_index[-li]"<<" the path to the saved index\n"
            <<"\t"<<"--output_index[-oi]"<<" write the index to a binary file\n"
            <<"\t"<<"--output_graph_binary[-ogb]"<<" write the graph in binary format\n"
            <<"\t"<<"--load_graph_binary[-lgb]"<<" load the graph from binary file\n"
            <<"\t"<<"-v"<<" verbosity\n" 

            <<"[PARAMETERS]: index\n\n"
            <<"\t"<<"-h[--help]"<<"\t prints this help message\n"
            <<"\t"<<"--input_graph[-ig]"<<"\t the path to the pangenome graph in gfa/fasta/binary format\n"
            <<"\t"<<"--k_mer_size[-k]"<<" the size of the k-mer\n"
            <<"\t -it"<<"\t the input is in text format (fasta/GFA)\n"
            <<"\t -ib"<<"\t the input is in binary format\n"
            <<"\t"<<"--output_file_name[-o]"<<" the name of the output file\n" 
            <<"\t"<<"-v"<<" verbosity\n" 
            

            <<"[PARAMETERS]: convert\n\n"
            <<"\t"<<"-h[--help]"<<"\t prints this help message\n"
            <<"\t"<<"--input_graph[-ig]"<<"\t the path to the pangenome graph in gfa/fasta/binary format\n"
            <<"\t -it"<<"\t the input is in text format (fasta/GFA)\n"
            <<"\t -ib"<<"\t the input is in binary format\n"
            <<"\t -of"<<"\t output in fasta format\n"
            <<"\t -ob"<<"\t output in binary format\n"
            <<"\t"<<"--output_file_name[-o]"<<" the name of the output file\n" 
            <<"\t"<<"-v"<<" verbosity\n" 
            ;
}
std::unordered_map<std::string,std::tuple<bool,std::string>> parseArgs(int argc,char **argv){
    std::unordered_map<std::string,std::tuple<bool,std::string>> arguments;
    arguments["graphfile"]=std::tuple<bool,std::string>(false,"");
    arguments["genomefile"]=std::tuple<bool,std::string>(false,"");
    arguments["genomefile2"]=std::tuple<bool,std::string>(false,"");//remove me
    arguments["kvalue"]=std::tuple<bool,std::string>(false,"");
    arguments["kmerfile"]=std::tuple<bool,std::string>(false,"");
    arguments["test"]=std::tuple<bool,std::string>(false,"");
    arguments["outputfilename"]=std::tuple<bool,std::string>(false,"ccdbgupdater");
    arguments["verbosity"]=std::tuple<bool,std::string>(false,"");
    arguments["updateindex"]=std::tuple<bool,std::string>(false,"");
    arguments["loadindex"]=std::tuple<bool,std::string>(false,"");
    arguments["outputindex"]=std::tuple<bool,std::string>(false,"");
    arguments["loadgraphbinary"]=std::tuple<bool,std::string>(false,"");
    arguments["outputgraphbinary"]=std::tuple<bool,std::string>(false,"");
    arguments["minimizer"]=std::tuple<bool,std::string>(false,"");
    arguments["smallbucket"]=std::tuple<bool,std::string>(false,"");
    if(argc<2){
        show_usage();
        exit(0);
    }
    if (!strcmp(argv[1], "--help") || !strcmp(argv[1],"-h")){
        show_usage();
        exit(0);
    }
    else if(!strcmp(argv[1], "update")){
        arguments["option"]=std::tuple<bool,std::string>(false,"update");
    }
    else if(!strcmp(argv[1], "index"))
    {
        arguments["option"]=std::tuple<bool,std::string>(false,"index");
    }
    else if(!strcmp(argv[1], "convert"))
    {
        arguments["option"]=std::tuple<bool,std::string>(false,"convert");
    }
    else{
        std::cout<<"unsupported option"<<std::endl<<std::endl;
        show_usage();
        exit(0);
    }
    for(int i=2;i<argc;i++){
        if(!strcmp(argv[i],"--help") || !strcmp(argv[i],"-h")){
            show_usage();
        }
        else if(!strcmp(argv[i],"-u")||!strcmp(argv[i],"--update_index"))
        {
            arguments["updateindex"]=std::tuple<bool,std::string>(true,"");
        }
        else if(!strcmp(argv[i],"-it"))
        {   
            arguments["inputtext"]=std::tuple<bool,std::string>(true,"");
        }
        else if(!strcmp(argv[i],"-ib"))
        {   
            arguments["inputbinary"]=std::tuple<bool,std::string>(true,"");
        }
        else if(!strcmp(argv[i],"-of"))
        {   
            arguments["outputfasta"]=std::tuple<bool,std::string>(true,"");
        }
        else if(!strcmp(argv[i],"-ob"))
        {   
            arguments["outputbinary"]=std::tuple<bool,std::string>(true,"");
        }
        else if(!strcmp(argv[i],"--test"))
        {
            arguments["test"]=std::tuple<bool,std::string>(true,"");
        }
        else if(!strcmp(argv[i],"-v"))
        {
            arguments["verbosity"]=std::tuple<bool,std::string>(true,"");
        }
        else if(!strcmp(argv[i],"-oi")||!strcmp(argv[i],"--output_index"))
        {
            arguments["outputindex"]=std::tuple<bool,std::string>(true,"");
        }
        else if(!strcmp(argv[i],"-ogb")||!strcmp(argv[i],"--output_graph_binary"))
        {
            arguments["outputgraphbinary"]=std::tuple<bool,std::string>(true,"");
        }
        else if (!strcmp(argv[i],"-lgb")||!strcmp(argv[i],"--load_graph_binary"))
        {
            if(i+1<argc){
                arguments["loadgraphbinary"]=std::tuple<bool,std::string>(true,argv[i+1]);
                i=i+1;
            }
            else{
                std::cerr<<"the input file name of the graph in binary format is missing, use --help or -h for details on how to use the program\n";
                show_usage();
                exit(0);
            }
        }
        else if (!strcmp(argv[i],"-o")||!strcmp(argv[i],"--output_file_name"))
        {
            if(i+1<argc){
                arguments["outputfilename"]=std::tuple<bool,std::string>(true,argv[i+1]);
                i=i+1;
            }
            else{
                std::cerr<<"the output file name is missing, use --help or -h for details on how to use the program\n";
                show_usage();
                exit(0);
            }
        }
        else if (!strcmp(argv[i],"-li")||!strcmp(argv[i],"--load_index"))
        {
            if(i+1<argc){
                arguments["loadindex"]=std::tuple<bool,std::string>(true,argv[i+1]);
                i=i+1;
            }
            else{
                std::cerr<<"The path to the index to be loaded is missing, please use -h[--help] for details\n";
                show_usage();
                exit(0);
            }
        }
        else if (!strcmp(argv[i],"--input_graph")||!strcmp(argv[i],"-ig"))
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
        else if(!strcmp(argv[i],"-m")||!strcmp(argv[i],"--minimizer_size"))
        {
            if(i+1<argc){
                arguments["minimizer"]=std::tuple<bool,std::string>(true,argv[i+1]);
                i=i+1;
            }
            else{
                std::cerr<<"the size of minimizer is missing, use --help or -h for details on how to use the program\n";
                show_usage();
                exit(0);
            } 
        }
        else if(!strcmp(argv[i],"-s")||!strcmp(argv[i],"--smallest_merge"))
        {
            if(i+1<argc){
                arguments["smallbucket"]=std::tuple<bool,std::string>(true,argv[i+1]);
                i=i+1;
            }
            else{
                std::cerr<<"the threshold value is missing, use --help or -h for details on how to use the program\n";
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
    if(!strcmp(std::get<1>(arguments["option"]).c_str(),"update")){
        if(!(std::get<0>(arguments["kvalue"])&&std::get<0>(arguments["graphfile"]))||!(std::get<0>(arguments["kmerfile"])||std::get<0>(arguments["genomefile"]))){
            std::cerr<<"Some required arguments are missing\n";
            show_usage();
            exit(0);
        }
    }
    else if(!strcmp(std::get<1>(arguments["option"]).c_str(),"index"))
    {
        if(!(std::get<0>(arguments["kvalue"])&&std::get<0>(arguments["graphfile"]))||!(arguments.count("inputtext")>0 || arguments.count("inputbinary")>0)){
            std::cerr<<"Some required arguments are missing\n";
            show_usage();
            exit(0);
        }
    }
    else if(!strcmp(std::get<1>(arguments["option"]).c_str(),"convert"))
    {
        if(!std::get<0>(arguments["graphfile"])||!((arguments.count("inputtext")>0 && arguments.count("outputbinary")>0)||(arguments.count("inputbinary")>0 && arguments.count("outputfasta")>0))){
            std::cerr<<"Either the graph file is missing or you are trying to convert binary to binary or text to text\n";
            show_usage();
            exit(0);
        }
    }

    return arguments;
}
void save_index(std::string filename,Index_mphf& ind){
    essentials::save(ind,filename.c_str());
}
void load_index(std::string filename,Index_mphf& ind){
    essentials::load(ind,filename.c_str());
}
int main(int argc,char **argv){
    //parse arguments
    std::unordered_map<std::string,std::tuple<bool,std::string>> arguments=parseArgs(argc,argv);
    int minimizer_size=7;//default value
    int smallest_bucket_merge=500;//default value
    if(std::get<0>(arguments["minimizer"])){
        minimizer_size=stoi(std::get<1>(arguments["minimizer"]));
    }
    if(std::get<0>(arguments["smallbucket"])){
        smallest_bucket_merge=stoi(std::get<1>(arguments["smallbucket"]));
    }
    if(!strcmp(std::get<1>(arguments["option"]).c_str(),"index")){
        GfaGraph g;
        Index_mphf ind=Index_mphf(stoi(std::get<1>(arguments["kvalue"])),minimizer_size,smallest_bucket_merge);//
        if(std::get<0>(arguments["inputtext"])){
            GfaGraph g2=g.LoadFromFile(std::get<1>(arguments["graphfile"]));
            ind.build(g2);
            save_index((std::get<1>(arguments["outputfilename"])+"_index.bin"),ind);
        }
        else if(std::get<0>(arguments["inputbinary"])){
            g.deserialize(std::get<1>(arguments["graphfile"]));
            ind.build(g);
            save_index((std::get<1>(arguments["outputfilename"])+"_index.bin"),ind);
        }        
    }
    else if(!strcmp(std::get<1>(arguments["option"]).c_str(),"convert"))
    {
        GfaGraph g;
        if(std::get<0>(arguments["inputtext"]) && std::get<0>(arguments["outputbinary"])){
            GfaGraph g2=g.LoadFromFile(std::get<1>(arguments["graphfile"]));
            g2.serialize(std::get<1>(arguments["outputfilename"])+".bin");
        }
        else if(std::get<0>(arguments["inputbinary"]) && std::get<0>(arguments["outputfasta"])){
            g.deserialize(std::get<1>(arguments["graphfile"]));
            g.convertToFasta(std::get<1>(arguments["outputfilename"])+".fa");
        }
    }
    else{//update the graph
        std::unordered_map<uint64_t,bool> k_mer;
        //load the input graph
        auto start_g=std::chrono::steady_clock::now();
        GfaGraph g;
        GfaGraph g2;//=g.LoadFromFile("bifrost_graph_101.fasta");
        if(std::get<0>(arguments["loadgraphbinary"])){
            g.deserialize(std::get<1>(arguments["loadgraphbinary"]));
            g2=g;
        }
        else{
            g2=g.LoadFromFile(std::get<1>(arguments["graphfile"]));
        }
        auto end_G=std::chrono::steady_clock::now();
        std::cout<<"Time to load the graph is "<<std::chrono::duration_cast<std::chrono::nanoseconds>(end_G - start_g).count()*1e-9<<std::endl;
        bool verbose=std::get<0>(arguments["verbosity"]);
        float time_kmtricks=0;
        //if we don't have a kmerfile, then a genome is passed
        if(!std::get<0>(arguments["kmerfile"])){
            std::string input_to_kmtricks=std::get<1>(arguments["graphfile"]);
            if(std::get<1>(arguments["graphfile"]).length()>4 && !std::strcmp(std::get<1>(arguments["graphfile"]).substr(std::get<1>(arguments["graphfile"]).length()-3).c_str(),"gfa")){
                //if the graph is not in fasta format, convert it to fasta as kmtricks accepts only fasta format (gzipped or not)
                g2.convertToFasta("unitigs_fasta.fa");
                input_to_kmtricks="unitigs_fasta.fa";
            }
            //call kmtricks
            std::system("chmod +x ../src/utils.sh");
            auto start=std::chrono::steady_clock::now();
            std::system(("bash ../src/utils.sh "+std::get<1>(arguments["kvalue"])+" "+input_to_kmtricks+" "+std::get<1>(arguments["genomefile"])+" "+std::get<1>(arguments["outputfilename"])+".txt").c_str());
            auto end =std::chrono::steady_clock::now();
            time_kmtricks=std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()*1e-9;
            //load the absent k-mers found by kmtricks
            k_mer=createHashTable(std::get<1>(arguments["outputfilename"])+".txt");
        }
        else{
            //load the kmers
            k_mer=createHashTable(std::get<1>(arguments["kmerfile"]));
        }
        //create the index of the graph
        Index_mphf ind=Index_mphf(stoi(std::get<1>(arguments["kvalue"])),minimizer_size,smallest_bucket_merge);//
        auto start=std::chrono::steady_clock::now();
        if(!std::get<0>(arguments["loadindex"])){
            ind.build(g2);
        }
        else{
            load_index((std::get<1>(arguments["outputfilename"])+"_index.bin"),ind);
        }
        auto end =std::chrono::steady_clock::now();
        float time_index = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()*1e-9;
        start=std::chrono::steady_clock::now();
        //construct the funitigs from the absent kmers
        std::unordered_map<uint64_t,bool> candidate_splits;//save split info u id,position
        std::unordered_map<int,std::pair<int,int>> candidate_join;//save join id u->(funitig 1,funitig 2)
        std::unordered_map<int,std::pair<int,int>> funitigs_join;//save join id u->(funitig 1,funitig 2) 
        std::unordered_map<int,Unitig> constrct_unitigs=construct_unitigs_from_kmer(ind,k_mer,stoi(std::get<1>(arguments["kvalue"])),g2,candidate_splits,candidate_join,funitigs_join);   
        std::unordered_map<uint64_t,std::vector<std::tuple<int,int,bool>>> f_index=index_constructed_unitigs(constrct_unitigs,31);  
        end=std::chrono::steady_clock::now();
        float time_construction=std::chrono::duration_cast<std::chrono::nanoseconds>(end - start).count()*1e-9;
        update_CdBG(ind,g2,constrct_unitigs,candidate_splits,candidate_join,std::get<0>(arguments["verbosity"]),std::get<0>(arguments["updateindex"]));
        insert_funitigs(g2,constrct_unitigs);
        if(std::get<0>(arguments["updateindex"])){
            ind.extract_kmers_from_funitigs(constrct_unitigs,g2);
            ind.update(g2);
        }
        if(std::get<0>(arguments["outputindex"])){
            save_index(std::get<1>(arguments["outputfilename"])+"_index.bin",ind);
        }
        if(std::get<0>(arguments["outputgraphbinary"])){
            g2.serialize(std::get<1>(arguments["outputfilename"])+"_graph.bin");
        }
        else{
            g2.convertToFasta(std::get<1>(arguments["outputfilename"])+".fa");
        }
    }
}
