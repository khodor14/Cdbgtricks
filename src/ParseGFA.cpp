#include <ParseGFA.hpp>
#include <iostream>
GfaGraph::GfaGraph:

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
            size_t id=stoi(idstr)
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
void convertToFasta(){
    string filenameout("out_1.fa");
	ofstream out(filenameout);
	for(uint i(0),i<this.nodes.size(),++i){
		out<<">sequence"+to_string(i)<<endl;
		out<<this.nodes[i].sequence<<endl;
	}
}