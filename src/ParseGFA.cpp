#include "ParseGFA.h"
#include "unitig.h"
#include <iostream>
#include <string>
#include <sstream>
#include <cassert>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <unordered_map>
#include <string_view>
#include <sys/stat.h>
#include <stdint.h>
#include <stdio.h>
#include <memory>
#include "zstr.hpp"
#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);
#endif
//Node::Node:
Node::Node(int id,std::string seq){
	this->id=id;
	this->unitig=seq;
}
int Node::get_id(){
	return id;
}
std::string Node::get_unitig(){
	return unitig;
}
//GfaGraph::GfaGraph:
edge::edge(int sink,int source, bool from,bool to){
	this->sink=sink;
	this->source=source;
	this->from=from;
	this->to=to;
}
int edge::get_sink(){
	return sink;
}
int edge::get_source(){
	return source;
}
void GfaGraph::set_max_node_id(int id){
	max_node_id=id;
}
GfaGraph GfaGraph::LoadFromFile(std::string filename){
	if(filename.length()>4 && !std::strcmp(filename.substr(filename.length()-3).c_str(),"gfa")){
		return LoadFromStream(filename,true);
	}
    else{
		return LoadFromStream(filename,false);
	}
}
GfaGraph GfaGraph::LoadFromStream(std::string filename,bool gfa){
	/*
	the input is a input stream file representing the graph in gfa or fasta
	if gfa is true then the graph is int gfa else it is in fasta
	*/
	GfaGraph graph;
	if(!gfa){
		int id=0;
		FILE* fp = fopen(filename.c_str(), "r");
		if(fp==0){
			std::cerr<<"Couldn't open the file "<<filename<<std::endl;
		}
		kseq_t* kseq;
    	kseq = kseq_init(gzopen(filename.c_str(),"r"));
		while(kseq_read(kseq)>=0){
			id++;
			graph.unitigs[id]=Unitig(kseq->seq.s);
		}
		graph.set_max_node_id(id);
		return graph;
	}
	std::unique_ptr<std::istream> graph_in=std::unique_ptr<std::istream>(new zstr::ifstream(filename, std::ios_base::in));
	std::istream& gin =*graph_in;
	std::string line;
	std::vector<std::string> line_fields;
	int max_id=0;
	std::ofstream graphfile;
	std::ostream graphf(0);
	graphfile.open(filename);
	graphf.rdbuf(graphfile.rdbuf());
	std::string seq;
	while (getline(gin, line).good())
	{
			//std::cout<<"Hello"<<std::endl;
		    if (line[0] == 'S'){ // Segment line

                const char* buffer = line.c_str() + 2;
                const char* end_buffer = line.c_str() + line.length();
                const char* prev_buffer = strchr(buffer, '\t');

                while (prev_buffer != NULL){

                    line_fields.push_back(std::string(buffer, prev_buffer - buffer));
                    buffer = prev_buffer + 1;
                    prev_buffer = strchr(buffer, '\t');
                }

                if (end_buffer - buffer != 0) line_fields.push_back(std::string(buffer, end_buffer - buffer));
                int id = stoi(move(line_fields[0]));
				if(id>max_id){
					max_id=id;
				}
				seq=move(line_fields[1]);
				graphf<<">sequence"+std::to_string(id)<<std::endl;
				graphf<<seq<<std::endl;
				graph.unitigs[id]=Unitig(seq);
				line_fields.clear();
		}
		/*if (line[0] == 'L')
		{
			std::stringstream sstr {line};
			std::string fromstr;
			std::string tostr;
			std::string fromstart;
			std::string toend;
			std::string dummy;
			std::string overlapstr;
			
			sstr >> dummy;
			assert(dummy == "L");
			sstr >> fromstr;
			size_t from = stoi(fromstr);
			sstr >> fromstart;
			sstr >> tostr;
			size_t to = stoi(tostr);
			sstr >> toend;
			assert(fromstart == "+" || fromstart == "-");
			assert(toend == "+" || toend == "-");
			sstr >> overlapstr;
			assert(overlapstr.size() >= 1);
			size_t charAfterIndex = 0;
			int overlap = std::stol(overlapstr, &charAfterIndex, 10);
			assert(overlap >= 0);
            bool fromP(true);
            bool toP(true);
            if(fromstart=="-"){
                fromP=false;
            }
            if(toend=="-"){
                toP=false;
            }
            edge e=edge(from,to,fromP,toP);

            graph.edges.push_back(e);
		}*/
	}
	graphfile.close();
	graph.set_max_node_id(max_id);
	return graph;
}
 void GfaGraph::convertToFasta(std::string filename){
	std::ofstream graphfile;
	std::ostream graph(0);

	graphfile.open(filename.c_str());
	graph.rdbuf(graphfile.rdbuf());
	for(std::pair<int,Unitig> unitig:unitigs){
		graph<<">sequence"+std::to_string(unitig.first)<<std::endl;
		graph<<unitig.second.to_string()<<std::endl;
	}
	graphfile.close();
}
int GfaGraph::get_max_node_id(){
	return max_node_id;
}
std::vector<int> GfaGraph::find_in_neighbors(int node_id){
	std::vector<int> in_neighbors;
	auto iter_edges=edges.begin();
	while(iter_edges!=edges.end()){
		if(iter_edges->get_source()==node_id){
			in_neighbors.push_back(iter_edges->get_sink());
		}
		iter_edges++;
	}
	return in_neighbors;
}
std::vector<int> GfaGraph::find_out_neighbors(int node_id){
	std::vector<int> out_neighbors;
	auto iter_edges=edges.begin();
	while(iter_edges!=edges.end()){
		if(iter_edges->get_sink()==node_id){
			out_neighbors.push_back(iter_edges->get_source());
		}
		iter_edges++;
	}
	return out_neighbors;	
}
void GfaGraph::delete_unitig(int id){
	unitigs.erase(id);
}
void GfaGraph::insert_unitig(int id,Unitig u){
	try{
		unitigs[id]=u;
	}
	catch(std::exception &e){
		std::cerr << "Exception caught : " << e.what() << std::endl;
	}
}
Unitig GfaGraph::get_unitig(int id){
	Unitig u=unitigs[id];
	if(unitigs.count(id)>0){
		u=unitigs[id];
	}
	else{
		std::cout<<"Unitig "<<id<< " does not exist"<<std::endl;
		exit(0);
	}
	return u;
}
void GfaGraph::serialize(const std::string filename) {
    std::ofstream ofs(filename, std::ios::binary);

    if (ofs.good()) {
        // Serialize the size of the map
        int size = static_cast<int>(unitigs.size());
		ofs.write(reinterpret_cast<const char*>(&max_node_id), sizeof(max_node_id));
        ofs.write(reinterpret_cast<const char*>(&size), sizeof(size));

        // Serialize each key-value pair in the map
        for (auto& pair : unitigs) {
            // Serialize the key
            int key = pair.first;

            ofs.write(reinterpret_cast<const char*>(&key), sizeof(key));

            // Serialize the value
            Unitig value = pair.second;
			uint8_t left_unused=value.get_left_unused();
			uint8_t right_unused=value.get_right_unused();
			std::vector<uint8_t> encode=value.get_encoding();
			ofs.write(reinterpret_cast<const char*>(&left_unused), sizeof(left_unused));
			ofs.write(reinterpret_cast<const char*>(&right_unused), sizeof(right_unused));
			int size_vec=encode.size();
			ofs.write(reinterpret_cast<const char*>(&size_vec), sizeof(size_vec));
			for(int i=0;i<size_vec;i++){
				uint8_t val=encode[i];
				ofs.write(reinterpret_cast<const char*>(&val), sizeof(val));
			}
        }

        ofs.close();
    }
    else {
        std::cout << "Unable to open binary file for serialization." << std::endl;
    }
}
void GfaGraph::deserialize(const std::string filename) {
    std::ifstream ifs(filename, std::ios::binary);

    if (ifs.good()) {
        // Deserialize the size of the map
        int size;
		ifs.read(reinterpret_cast<char*>(&max_node_id), sizeof(max_node_id));
        ifs.read(reinterpret_cast<char*>(&size), sizeof(size));

        // Deserialize each key-value pair in the map
        for (int i = 0; i < size; i++) {
            // Deserialize the key
            int key;
            ifs.read(reinterpret_cast<char*>(&key), sizeof(key));
			uint8_t left_unused;
			uint8_t right_unused;
			int size_vec;
			ifs.read(reinterpret_cast<char*>(&left_unused), sizeof(left_unused));
			ifs.read(reinterpret_cast<char*>(&right_unused), sizeof(right_unused));
			ifs.read(reinterpret_cast<char*>(&size_vec), sizeof(size_vec));
            // Deserialize the value
            Unitig value;
			std::vector<uint8_t> vec(size_vec);
			for(int j=0;j<size_vec;j++){
				uint8_t val;
				ifs.read(reinterpret_cast<char*>(&val), sizeof(val));
				vec[j]=val;
			}
            // Insert the key-value pair into the map
            unitigs[key] = Unitig(left_unused,right_unused,vec);
        }

        ifs.close();
    }
    else {
        std::cout << "Unable to open binary file for deserialization." << std::endl;
    } 
	 
}
void GfaGraph::fixe_edges(int node_id,int new_node, bool from, bool to){

}
std::unordered_map<int,Unitig> GfaGraph::get_nodes(){
	return unitigs;
}
int GfaGraph::test_kmer_presence(uint64_t kmer,uint64_t suffix,uint64_t kmer_pos,int k){
	/*
	Input:
		a k-mer
		its position info (64 bits):see below
		its size
	return if this k-mer is in the graph and if we encountered a split position along with the split position
	this function is needed as an mphf will return a valid identifier for a given k-mer, to validate that this k-mer
	is not an alien one, we compare it with the one at the provided position
	*/
	//assuming kmer is sent in its canonical form
	/*Representation of k-mer position
	|32 bits for unitig id|31 bits for position in this unitig|1 bit for orientation|
	*/
	int id_unitig_graph=(int)(kmer_pos>>32);//the id of the unitig having the first occurence
    int position_unitig=(int)((kmer_pos>>1)&0x7FFFFFFF);//the position of k-mer in this unitig
	bool orient=kmer_pos&1;//the orientation
	Unitig seq=unitigs[id_unitig_graph];
	uint64_t kmer_in_graph=seq.get_ith_mer(position_unitig,k);//take the k-mer from the unitig
	//either same orientation or different orientation
	//in the case of different orientation, we compare the kmer to the canonical form of the k-mer in graph
	if((orient && kmer_in_graph==kmer)||(!orient && kmer==canonical_bits(kmer_in_graph,k))){
		/*
		The k-mer (i.e. the suffix as well) is in the graph with 3-possible cases
					1- Case split easy to detect (we check that the suffix of kmer is neither the suffix nor the prefix of unitig)
					2- 
		*/
		uint64_t pref=kmer_in_graph>>2;
		uint64_t suff=kmer_in_graph&(0xffffffffffffffff>>(66-2*k));
		if(position_unitig>0&&position_unitig<seq.unitig_length()-k){
			//the kmer is not the k-suffix not the k-prefix => we encounter a split position=>let's compute the split position
			if(canonical_bits(suffix,k-1)==canonical_bits(suff,k-1)){
				return position_unitig+1;
			}
			else{
				return position_unitig;
			}
		}
		else{
			//making sure that we don't have a split position
			
			if(position_unitig==0){
				//compare suffix to suff
				if(canonical_bits(suffix,k-1)==canonical_bits(suff,k-1)){
					return 1;//we have to split at position=1
				}	
			}
			else{
				if(canonical_bits(pref,k-1)==canonical_bits(suffix,k-1)){
					return seq.unitig_length()-k;//we split at position |u|
				}
			}
			return 0;
		}

	}
	else{
		return -1;//the k-mer is not in the graph
	}
}