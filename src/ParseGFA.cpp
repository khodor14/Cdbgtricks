#include "ParseGFA.h"


#ifndef KSEQ_INIT_READY
#define KSEQ_INIT_READY


KSEQ_INIT(gzFile, gzread);
#endif
//Node::Node:
Node::Node(int id, std::string seq) {
  this -> id = id;
  this -> unitig = seq;
}
int Node::get_id() {
  return id;
}
std::string Node::get_unitig() {
  return unitig;
}
//GfaGraph::GfaGraph:
edge::edge(int sink, int source, bool from, bool to) {
  this -> sink = sink;
  this -> source = source;
  this -> from = from;
  this -> to = to;
}
int edge::get_sink() {
  return sink;
}
int edge::get_source() {
  return source;
}
void GfaGraph::set_max_node_id(int id) {
  max_node_id = id;
}
GfaGraph GfaGraph::LoadFromFile(std::string filename) {
  if (filename.length() > 4 && !std::strcmp(filename.substr(filename.length() - 3).c_str(), "gfa")) {
    return LoadFromStream(filename, true);
  } else {
    return LoadFromStream(filename, false);
  }
}
GfaGraph GfaGraph::LoadFromStream(std::string filename, bool gfa) {
  /*
  the input is a input stream file representing the graph in gfa or fasta
  if gfa is true then the graph is int gfa else it is in fasta
  */
  GfaGraph graph;
  if (!gfa) {
    int id = 0;
    FILE * fp = fopen(filename.c_str(), "r");
    if (fp == 0) {
      std::cerr << "Couldn't open the file " << filename << std::endl;
    }
    kseq_t * kseq;
    kseq = kseq_init(gzopen(filename.c_str(), "r"));
    while (kseq_read(kseq) >= 0) {
      id++;
      graph.unitigs[id] = Unitig(kseq -> seq.s);
    }
    graph.set_max_node_id(id);
    return graph;
  }
  std::unique_ptr < std::istream > graph_in = std::unique_ptr < std::istream > (new zstr::ifstream(filename, std::ios_base::in));
  std::istream & gin = * graph_in;
  std::string line;
  std::vector < std::string > line_fields;
  int max_id = 0;
  std::string seq;
  while (getline(gin, line).good()) {
    if (line[0] == 'S') { // Segment line

      const char * buffer = line.c_str() + 2;
      const char * end_buffer = line.c_str() + line.length();
      const char * prev_buffer = strchr(buffer, '\t');

      while (prev_buffer != NULL) {

        line_fields.push_back(std::string(buffer, prev_buffer - buffer));
        buffer = prev_buffer + 1;
        prev_buffer = strchr(buffer, '\t');
      }

      if (end_buffer - buffer != 0) line_fields.push_back(std::string(buffer, end_buffer - buffer));
      move(line_fields[0]);
      max_id++;
      seq = move(line_fields[1]);
      graph.unitigs[max_id] = Unitig(seq);
      line_fields.clear();
    }
  }
  graph.set_max_node_id(max_id);
  return graph;
}
void GfaGraph::convertToFasta(std::string filename) {
  std::ofstream graphfile;
  std::ostream graph(0);

  graphfile.open(filename.c_str());
  graph.rdbuf(graphfile.rdbuf());
  for (int i=1;i<=unitigs.size();i++) {
    graph << "> "<<std::to_string(i) << std::endl;
    graph << unitigs[i].to_string() << std::endl;
  }
  graphfile.close();
}
int GfaGraph::get_max_node_id() {
  return max_node_id;
}
std::vector < int > GfaGraph::find_in_neighbors(int node_id) {
  std::vector < int > in_neighbors;
  auto iter_edges = edges.begin();
  while (iter_edges != edges.end()) {
    if (iter_edges -> get_source() == node_id) {
      in_neighbors.push_back(iter_edges -> get_sink());
    }
    iter_edges++;
  }
  return in_neighbors;
}
std::vector < int > GfaGraph::find_out_neighbors(int node_id) {
  std::vector < int > out_neighbors;
  auto iter_edges = edges.begin();
  while (iter_edges != edges.end()) {
    if (iter_edges -> get_sink() == node_id) {
      out_neighbors.push_back(iter_edges -> get_source());
    }
    iter_edges++;
  }
  return out_neighbors;
}
void GfaGraph::delete_unitig(int id) {
  unitigs.erase(id);
}
void GfaGraph::insert_unitig(int id, Unitig u) {
  try {
    unitigs[id] = u;
  } catch (std::exception & e) {
    std::cerr << "Exception caught : " << e.what() << std::endl;
  }
}
Unitig GfaGraph::get_unitig(int id) {
  Unitig u = unitigs[id];
  if (unitigs.count(id) > 0) {
    u = unitigs[id];
  } else {
    std::cout << "Unitig " << id << " does not exist" << std::endl;
    exit(0);
  }
  return u;
}
void GfaGraph::serialize(const std::string filename) {
  std::ofstream ofs(filename, std::ios::binary);

  if (ofs.good()) {
    // Serialize the size of the map
    int size = static_cast < int > (unitigs.size());
    ofs.write(reinterpret_cast <
      const char * > ( & max_node_id), sizeof(max_node_id));
    ofs.write(reinterpret_cast <
      const char * > ( & size), sizeof(size));

    // Serialize each key-value pair in the map
    for (auto & pair: unitigs) {
      // Serialize the key
      int key = pair.first;

      ofs.write(reinterpret_cast <
        const char * > ( & key), sizeof(key));

      // Serialize the value
      Unitig value = pair.second;
      uint8_t left_unused = value.get_left_unused();
      uint8_t right_unused = value.get_right_unused();
      std::vector < uint8_t > encode = value.get_encoding();
      ofs.write(reinterpret_cast <
        const char * > ( & left_unused), sizeof(left_unused));
      ofs.write(reinterpret_cast <
        const char * > ( & right_unused), sizeof(right_unused));
      int size_vec = encode.size();
      ofs.write(reinterpret_cast <
        const char * > ( & size_vec), sizeof(size_vec));
      for (int i = 0; i < size_vec; i++) {
        uint8_t val = encode[i];
        ofs.write(reinterpret_cast <
          const char * > ( & val), sizeof(val));
      }
    }

    ofs.close();
  } else {
    std::cout << "Unable to open binary file for serialization." << std::endl;
  }
}
void GfaGraph::deserialize(const std::string filename) {
  std::ifstream ifs(filename, std::ios::binary);

  if (ifs.good()) {
    // Deserialize the size of the map
    int size;
    ifs.read(reinterpret_cast < char * > ( & max_node_id), sizeof(max_node_id));
    ifs.read(reinterpret_cast < char * > ( & size), sizeof(size));

    // Deserialize each key-value pair in the map
    for (int i = 0; i < size; i++) {
      // Deserialize the key
      int key;
      ifs.read(reinterpret_cast < char * > ( & key), sizeof(key));
      uint8_t left_unused;
      uint8_t right_unused;
      int size_vec;
      ifs.read(reinterpret_cast < char * > ( & left_unused), sizeof(left_unused));
      ifs.read(reinterpret_cast < char * > ( & right_unused), sizeof(right_unused));
      ifs.read(reinterpret_cast < char * > ( & size_vec), sizeof(size_vec));
      // Deserialize the value
      Unitig value;
      std::vector < uint8_t > vec(size_vec);
      for (int j = 0; j < size_vec; j++) {
        uint8_t val;
        ifs.read(reinterpret_cast < char * > ( & val), sizeof(val));
        vec[j] = val;
      }
      // Insert the key-value pair into the map
      unitigs[key] = Unitig(left_unused, right_unused, vec);
    }

    ifs.close();
  } else {
    std::cout << "Unable to open binary file for deserialization." << std::endl;
  }

}
void GfaGraph::fixe_edges(int node_id, int new_node, bool from, bool to) {

}
std::unordered_map < int, Unitig > GfaGraph::get_nodes() {
  return unitigs;
}
bool GfaGraph::test_kmer(uint64_t kmer,uint64_t kmer_pos,int k){
  int id_unitig_graph = kmer_pos >> 32; //the id of the unitig having the first occurence
  int position_unitig = (kmer_pos >> 1) & 0x7FFFFFFF; //the position of k-mer in this unitig
  Unitig seq = unitigs[id_unitig_graph];
  uint64_t kmer_in_graph = seq.get_ith_mer(position_unitig, k); //take the k-mer from the unitig
  return canonical_bits(kmer_in_graph,k)==canonical_bits(kmer,k);
}
int GfaGraph::test_kmer_presence(uint64_t kmer, uint64_t suffix, uint64_t kmer_pos, int k) {
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
  int id_unitig_graph = kmer_pos >> 32; //the id of the unitig having the first occurence
  int position_unitig = (kmer_pos >> 1) & 0x7FFFFFFF; //the position of k-mer in this unitig
  Unitig seq = unitigs[id_unitig_graph];
  uint64_t kmer_in_graph = seq.get_ith_mer(position_unitig, k); //take the k-mer from the unitig
  std::tuple < uint64_t, bool > kmer_query_canonical = reverseComplementCanonical(kmer, k);
  std::tuple < uint64_t, bool > kmer_graph_canonical = reverseComplementCanonical(kmer_in_graph, k);
  if (std::get < 0 > (kmer_query_canonical) == std::get < 0 > (kmer_graph_canonical)) {
    //the k-mer is in the graph
    uint64_t pref_kmer = kmer_in_graph >> 2;
    if(position_unitig == 0){
      if(suffix == pref_kmer){
        return 0;//maybe a join with the unitig
      }
      else if (canonical_bits(suffix,k-1) != canonical_bits(pref_kmer,k-1))
      {
        return 1;
      }
      return -2;
    }
    else if (position_unitig == seq.unitig_length()-k)
    {
      uint64_t suff_kmer = kmer_in_graph & (0xFFFFFFFFFFFFFFFF >> (66-2*k));
      if (suff_kmer == suffix){
        return -2;
      }
      else if (canonical_bits(pref_kmer, k - 1) == canonical_bits(suffix, k - 1))
      {
        return position_unitig;
      }
      else{
        return -1;// we may join
      }
    }
    else{
      pref_kmer = canonical_bits(pref_kmer, k - 1);
      if (suffix == pref_kmer) {
        return position_unitig;
      }
      return position_unitig + 1;
    }
  } else {
    return -3; //the k-mer is not in the graph
  }

}
uint64_t GfaGraph::get_kmer(uint64_t position, int k) {
  //takes the position of the k-mer (unitig id,unitig position,orientation)
  //return the k-mer at that position
  int id_unitig_graph = (int)(position >> 32); //the id of the unitig having the first occurence
  int position_unitig = (int)((position >> 1) & 0x7FFFFFFF); //the position of k-mer in this unitig
  Unitig seq = unitigs[id_unitig_graph];
  if (position_unitig > seq.unitig_length() - k) {
    std::cout << "out of range\n";
  }
  uint64_t kmer_in_graph = seq.get_ith_mer(position_unitig, k); //take the k-mer from the unitig
  return kmer_in_graph;
}