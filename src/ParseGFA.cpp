#include "ParseGFA.h"
#include <iostream>
#include <string>
#include <sstream>
#include <cassert>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cstring>

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
GfaGraph GfaGraph::LoadFromFile(std::string filename){
    std::ifstream file{filename, std::ios::in};
	if(filename.length()>4 && !std::strcmp(filename.substr(filename.length()-3).c_str(),"gfa")){
		std::cout<<filename<<" gfa\n";
		return LoadFromStream(file,true);
	}
    else{
		std::cout<<filename<<" fasta\n";
		return LoadFromStream(file,false);
	}
}
GfaGraph GfaGraph::LoadFromStream(std::ifstream &file,bool gfa){
	/*
	the input is a input stream file representing the graph in gfa or fasta
	if gfa is true then the graph is int gfa else it is in fasta
	*/
	GfaGraph graph;
	if(!gfa){
		int id=1;
		while (file.good())
		{
			std::string line = "";
			std::getline(file, line);
			if (line.size() == 0 && !file.good()) break;
			if (line.size() == 0) continue;
			if(line[0]=='>') continue;	
			graph.unitigs[id]=line;
			id++;
		}
		return graph;
	}
	while (file.good())
	{
		std::string line = "";
		std::getline(file, line);
		if (line.size() == 0 && !file.good()) break;
		if (line.size() == 0) continue;
		if (line[0] != 'S' && line[0] != 'L') continue;
		if (line[0] == 'S')
		{
			std::stringstream sstr;
			sstr.str(line);
			std::string idstr;
			std::string dummy;
			std::string seq;
			sstr >> dummy;
			assert(dummy == "S");
			sstr >> idstr;
			sstr >> seq;
			assert(seq.size() >= 1);
            int id=stoi(idstr);
			graph.unitigs[id]=seq;
		}
		if (line[0] == 'L')
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
		}
	}

	return graph;
}
 void GfaGraph::convertToFasta(){
    std::string filenameout("out_1.fa");
	std::ofstream out(filenameout);
	for(uint i(0);i<unitigs.size();++i){
		out<<">sequence"+std::to_string(i)<<std::endl;
		out<<unitigs[i]<<std::endl;
	}
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
void GfaGraph::fixe_edges(int node_id,int new_node, bool from, bool to){

}
std::unordered_map<int,std::string> GfaGraph::get_nodes(){
	return unitigs;
}