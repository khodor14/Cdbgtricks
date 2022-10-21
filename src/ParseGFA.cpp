#include <ParseGFA.hpp>
#include <iostream>

Node::Node :
int Node::get_id(){
	return id;
}
std::string Node::get_unitig(){
	return unitig;
}
GfaGraph::GfaGraph :

GfaGraph GfaGraph::LoadFromFile(std::string filename){
    std::ifstream file{filename};
    return LoadFromStream(file);
}

GfaGraph GfaGraph::LoadFromStream(std::istream& file){
	GfaGraph graph;
	while (file.good())
	{
		std::string line = "";
		std::getline(file, line);
		if (line.size() == 0 && !file.good()) break;
		if (line.size() == 0) continue;
		if (line[0] != 'S' && line[0] != 'L') continue;
		if (line[0] == 'S')
		{
			std::stringstream sstr {line};
			std::string idstr;
			std::string dummy;
			std::string seq;
			sstr >> dummy;
			assert(dummy == "S");
			sstr >> idstr;
			sstr >> seq;
            size_t id=stoi(idstr);
            structNode node;
            node.id=id;
            node.sequence=seq;
			assert(seq.size() >= 1);
			graph.nodes.push_back(node);
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
			int overlap = 0;
			sstr >> dummy;
			assert(dummy == "L");
			sstr >> fromstr;
			size_t from = stoi(fromstr);
			sstr >> fromstart;
			sstr >> tostr;
			size_t to = stoi(tostr)
			sstr >> toend;
			assert(fromstart == "+" || fromstart == "-");
			assert(toend == "+" || toend == "-");
			sstr >> overlapstr;
			assert(overlapstr.size() >= 1);
			size_t charAfterIndex = 0;
			overlap = std::stol(overlapstr, &charAfterIndex, 10);
			assert(overlap >= 0);
            bool fromP(true);
            bool toP(true);
            if(fromstart=="-"){
                fromP=false;
            }
            if(toend=="-"){
                toP=false;
            }
            edge e;
            e.sink=from;
            e.source=to;
            e.from=fromP;
            e.to=toP;
            graph.edges.push_back(e);
		}
	}

	return result;
}
 void GfaGraph::convertToFasta(){
    string filenameout("out_1.fa");
	ofstream out(filenameout);
	for(uint i(0),i<this.nodes.size(),++i){
		out<<">sequence"+to_string(i)<<endl;
		out<<this.nodes[i].sequence<<endl;
	}
}
std::vector<int> GfaGraph::find_in_neighbors(int node_id){
	std::vector<int> in_neighbors;
	auto iter_edges=edges.begin();
	while(iter_edges!=edges.end()){
		if(iter_edges->source==node_id){
			in_neighbors.push_back(iter_edges->sink);
		}
		iter_edges++;
	}
	return in_neighbors;
}
std::vector<int> GfaGraph::find_out_neighbors(int node_id){
	std::vector<int> out_neighbors;
	auto iter_edges=edges.begin();
	while(iter_edges!=edges.end()){
		if(iter_edges->sink==node_id){
			out_neighbors.push_back(iter_edges->source);
		}
		iter_edges++;
	}
	return out_neighbors;	
}
void GfaGraph::fixe_edges(int node_id,int new_node, bool from, bool to){

}
void GfaGraph::fixe_node(int node_id,int pos_modification){

}
int GfaGraph::create_nodes(std::string sequence){//return the node id
	int new_id=nodes.size()+1;
	Node new_node;
	new_node.ide=new_id;
	new_node.unitig=sequence;
	nodes.push_back(new_id);
	return new_id;
}
void GfaGraph::create_edge(int node_id,int new_node_id,bool from,bool to){

}
std::vector<Node> GfaGraph::get_nodes(){
	return nodes;
}