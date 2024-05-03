#include "unitigmodification.h"

const char alphabet[] = {
  'A',
  'C',
  'G',
  'T'
};
const char reversealphabet[] = {
  'T',
  'G',
  'C',
  'A'
};
std::vector < char > possible_right_extension(uint64_t mer, int k, std::unordered_map < uint64_t, bool > & kmers_to_add_to_graph, Index_mphf & ind, GfaGraph & graph, std::unordered_map < uint64_t, bool > & candidate_splits, uint64_t * possible_join) {
  /*
  the input is :a string (k-1)-mer to be extended from the right
               :an unordered map containing the canonical forms of k-mers as keys and boolean value as values
                                                                                      the boolean is used to say that the k-mer had been used to construct another unitig
  the output is a vector containing the possible characters that can extend the (k-1)-mer to a k-mer 

  example: if the input is mer=ACTG and the map is {ACTGA:false,ACTGT:false,TTGCT:false}
  then the ouput should be <A,T>
  */

  std::vector < char > possible_character_for_extension;
  uint64_t query;
  size_t number_join = 0;
  uint64_t join_position = 0;
  uint64_t canonical_km = canonical_bits(mer, k - 1);
  for (int i = 0; i < 4; i++) {
    //the canonical of mer+character from the alphabet is present
    uint64_t next_kmer = ((mer << 2) & (0xffffffffffffffff >> (64 - 2 * k))) | baseToInt(alphabet[i]);
    query = canonical_bits(next_kmer, k);
    if (kmers_to_add_to_graph.count(query) > 0) { //&& !kmers_to_add_to_graph[query]
      //add the corresponding character to the vector
      if (query == next_kmer) {
        possible_character_for_extension.push_back(alphabet[i]);
      } else {
        possible_character_for_extension.push_back(reversealphabet[i]);
      }
    } else {
      uint64_t kmer_pos = ind.kmer_position(query); //get the k-mer position(id,position in unitig,orientation)
      if (kmer_pos != 0) { //we found an mphf corresponding to the minimizer of the k-mer
        int may_exist = graph.test_kmer_presence(query, canonical_km, kmer_pos, ind.get_k_length());
        if (may_exist > 0) { //we encountered a split position we add the position to a vector which will be sorted and processed
          kmer_pos = kmer_pos >> 32;
          kmer_pos = (kmer_pos << 32) | may_exist;
          candidate_splits[kmer_pos] = true;
          return std::vector < char > ();
        } else if (may_exist == 0 && number_join < 2) { //the k-mer is in the graph and we encountered a possible join
          /*
          add some conditions to check if it is a join
          how many possible joins are found=> if 1 then maybe join => let from left extention decide
          when the value of may_exist is 0 then we encounter a suffix prefix overlap
          which may result in a join case
          we stop at the second occurence of suffix-prefix overlap because it is no more a join
          the conditions outside the loop will correct the value of join_position
          */
          number_join++;
          join_position = kmer_pos;
        }
      }
    }
  }
  if (number_join > 1 || possible_character_for_extension.size() > 1) {
    join_position = 0;
  } else if (number_join == 0 && possible_character_for_extension.size() == 1) {
    join_position = 2;
  } else if (number_join == 1 && possible_character_for_extension.size() == 1) {
    join_position = 3;
  }
  * possible_join = join_position; //maybe a join to be tested in the calling function
  return possible_character_for_extension;
}
std::string extend_right(std::string kmer, uint64_t mer, Index_mphf & unitigs_index, std::unordered_map < uint64_t, bool > & kmers_to_add_to_graph, GfaGraph & graph, std::unordered_map < uint64_t, bool > & candidate_splits, std::unordered_map < int, std::pair < int, int >> & candidate_join, int id, bool extended_from_right, std::unordered_map < int, std::pair < int, int >> & funitigs_joins) {
  /*
  the input is: a string k-mer representing the k-mer to be extended from right
              :the index of the unitigs, i.e. (k-1)-mer -> (unitig id,position,orientation)
              :an unordered map containing the canonical forms of k-mers as keys and boolean value as values
                                                                                      the boolean is used to say that the k-mer had been used to construct another unitig
  the outpus is a extended string
  */
  std::string extended_kmer = kmer;
  uint64_t suffix = mer & (0xffffffffffffffff >> (66 - 2 * unitigs_index.get_k_length())); //taking the suffix of the kmer to extend it
  bool flag = true;
  while (flag) //as long as we can extended the (k-1)-mer suffix of the k-mer
  {
    /*
    if the (k-1)-mer suffix of the k-mer is not in the graph
    we can try to extends
    if possible_character_to_extend contains nothing then the suffix is in the graph=>here we tested 4 cases the other 4 are tested in the second if
    if it contains more than one character then we have more possibilies ==> we can't extend
    otherwise we try to extend if this suffix is not extendable from left
    */
    uint64_t possible_join = 1; //this is an undefined join position as a valid join position will start at 0x100000000
    uint64_t possible_join_2 = 1;
    /*
     if a join is encountered the function possible_right_extension will change the value to a valid join position
     if more than one possible join or more than one possible extension the function will change the value to zero which is undefined
    */
    std::vector < char > possible_character_to_extend = possible_right_extension(suffix, unitigs_index.get_k_length(), kmers_to_add_to_graph, unitigs_index, graph, candidate_splits, & possible_join);
    std::vector < char > possible_character_to_extend_left = possible_right_extension(reverse_complement(suffix, unitigs_index.get_k_length() - 1), unitigs_index.get_k_length(), kmers_to_add_to_graph, unitigs_index, graph, candidate_splits, & possible_join_2);
    if (possible_character_to_extend.size() == 1) {
      //we cannot extend if the suffix of the k-mer has left extension

      /*
      if the k-mer is ACTGA and we already have TCTGA in the k-mer set(to be added to the graph)
      then we cannot extend ACTGA from the right
      N.B:to test this we send the reverse complement of the suffix and we try to extend it from right
      */
      if (possible_character_to_extend_left.size() == 1) {
        //now we find the extensions of the actual suffix
        //to extend we should have only one character in the vector
        if (possible_join == 2 && possible_join_2 == 2) {
          char character_to_add = possible_character_to_extend.front(); //get the character to augment it to our string result from the vector
          //std::string key=suffix+character_to_add; 
          uint64_t k_mer_from_funitigs = canonical_bits((suffix << 2) | baseToInt(character_to_add), unitigs_index.get_k_length());
          if (!kmers_to_add_to_graph[k_mer_from_funitigs]) {
            kmers_to_add_to_graph[k_mer_from_funitigs] = true;
            extended_kmer = extended_kmer + character_to_add; //append it to the string only once, the append function takes as argument the number of times we need to append to a string
            suffix = ((suffix << 2) | baseToInt(character_to_add)) & (0xffffffffffffffff >> (66 - 2 * unitigs_index.get_k_length()));
            //suffix=extended_kmer.substr(extended_kmer.length()-k_mer_length-1,extended_kmer.length());//takes the (k-1)-mer suffix of the extended k-mer
          } else {
            flag = false;
          }
        } else {
          //we have zero or more than 1 possible right extension of the (k-1)-mer
          break;
        }
      } else {
        //the (k-1)-mer suffix of the string has more than one left parent extensions
        break;
      }
    } else {
      if (possible_character_to_extend.size() == 0 && possible_character_to_extend_left.size() == 1 && possible_join_2 == 2 && possible_join > 3) {
        //we encountered a join
        int id_u = possible_join >> 32;
        int position_unitig = (int)((possible_join >> 1) & 0x7FFFFFFF);
        if (funitigs_joins.count(id) == 0) {
          funitigs_joins[id] = std::pair < int, int > (0, 0);
        }
        if (candidate_join.count(id_u) == 0) {
          candidate_join[id_u] = std::pair < int, int > (0, 0);
        }
        if (position_unitig == 0) {
          candidate_join[id_u].first = id;
        } else {
          candidate_join[id_u].second = id;
        }
        if (extended_from_right) {
          funitigs_joins[id].second = id_u;
        } else {
          funitigs_joins[id].first = id_u;
        }
      }
      break;
    }
  }
  //return the resultant k-mer
  return extended_kmer;
}
std::string unitig_from_kmers_string(std::string kmer, uint64_t mer, Index_mphf & unitig_index, std::unordered_map < uint64_t, bool > & kmers, GfaGraph & graph, int id, std::unordered_map < uint64_t, bool > & candidate_splits, std::unordered_map < int, std::pair < int, int >> & candidate_join, std::unordered_map < int, std::pair < int, int >> & funitigs_joins) {

  std::string right = extend_right(kmer, mer, unitig_index, kmers, graph, candidate_splits, candidate_join, id, true, funitigs_joins); //extend the right of the kmer
  std::string left = reverseComplement(extend_right(reverseComplement(kmer), reverse_complement(mer, unitig_index.get_k_length()), unitig_index, kmers, graph, candidate_splits, candidate_join, id, false, funitigs_joins)); //extend the left, extend right of reverse then reverse the result
  return left + right.substr(kmer.length(), right.length()); //return the concatenation of left and right. Omit the first k characters of right because they are the suffix of left
}
Unitig unitig_from_kmers(uint64_t kmer, Index_mphf & unitig_index, std::unordered_map < uint64_t, bool > & kmers, GfaGraph & graph, int id, std::unordered_map < uint64_t, bool > & candidate_splits, std::unordered_map < int, std::pair < int, int >> & candidate_join, std::unordered_map < int, std::pair < int, int >> & funitigs_joins) {
  /*
  create a unitig from the bit encoding of the k-mer
  Input: uint64 representing the kmer, the index of graph and all kmers
  Return a unitig(8 bits per 4 bases,how many left unused bits(0 in this case) and how many right unused bits)
  */
  return Unitig(unitig_from_kmers_string(to_string(kmer, unitig_index.get_k_length()), kmer, unitig_index, kmers, graph, id, candidate_splits, candidate_join, funitigs_joins));
}
std::unordered_map < int, Unitig > construct_unitigs_from_kmer(Index_mphf & unitig_index, std::unordered_map < uint64_t, bool > & kmers, int k,
  GfaGraph & graph, std::unordered_map < uint64_t, bool > & candidate_splits,
  std::unordered_map < int, std::pair < int, int >> & candidate_join, std::unordered_map < int, std::pair < int, int >> funitigs_join) {
  /*
  from k-mers represented as 64 bits
  */
  std::unordered_map < int, Unitig > unitigs_map; //storing them in a map for fast access
  //std::unordered_map<uint64_t,bool> candidate_splits;//save split info u id,position
  //std::unordered_map<int,std::pair<int,int>> candidate_join;//save join id u->(funitig 1,funitig 2)
  //std::unordered_map<int,std::pair<int,int>> funitigs_join;//save join id u->(funitig 1,funitig 2) 
  int i = 1;
  for (std::pair < uint64_t, bool > element: kmers) {
    if (!element.second) { //if this k-mer was not used in any unitig
      kmers[element.first] = true;
      Unitig res = unitig_from_kmers(element.first, unitig_index, kmers, graph, i, candidate_splits, candidate_join, funitigs_join); //create unitig from this k-mer
      unitigs_map[i] = res; //put it back to map
      i++;
    }
  }
  return unitigs_map;
}
void update_CdBG(Index_mphf & graph_index, GfaGraph & graph, std::unordered_map < int, Unitig > & constructed_unitigs, std::unordered_map < uint64_t, bool > & candidate_splits, std::unordered_map < int, std::pair < int, int >> & candidate_join, bool update_index) {
  std::vector < uint64_t > sorted_splits(candidate_splits.size()); //to be used for sorting
  int i = 0;
  for (auto elem: candidate_splits) {
    sorted_splits[i] = elem.first;
    i++;
  }
  std::sort(sorted_splits.begin(), sorted_splits.begin() + sorted_splits.size()); //sort all splits so the positions of splits in the same unitig will be in consecutive order
  std::unordered_map < int, std::pair < int, int >> split_to_join; //if a unitig is a canididate join and split, then the join will occur from left or right (we need to save temporarily these info)
  Unitig seq; //unitig to be splitted at different positions
  std::vector < uint8_t > seq_bases; //its bases in binary
  uint8_t left_unused, right_unsed; //left and right unused of this unitig
  int id = 0; //its id
  bool first_split = true; //helpful to take the the prefix of the unitig
  int new_id = graph.get_max_node_id();
  std::cout << "splits " << candidate_splits.size() << "\n";
  for (int i = 0; i < sorted_splits.size(); i++) {
    uint64_t split = sorted_splits[i];
    int current_id = split >> 32;
    int pos = split & 0xFFFFFFFF;
    //update if we observed different unitig
    if (current_id != id) { //for the first iteration this is true since current_id!=0
      id = current_id;
      seq = graph.get_unitig(id);
      seq_bases = seq.get_encoding();
      left_unused = seq.get_left_unused();
      right_unsed = seq.get_right_unused();
      first_split = true; // we encountered a new unitig so we have a first split
      if (candidate_join.count(current_id) > 0) {
        split_to_join[current_id] = std::pair < int, int > (0, 0);
      }
    }
    int next_pos = 0;
    //try to get next split position. we need to be at most in the second last position and the next split position should belong to the same unitig(same id)
    if (i < sorted_splits.size() - 1) {
      uint64_t next_split = sorted_splits[i + 1];
      int next_id = next_split >> 32;
      if (next_id == current_id) {
        next_pos = next_split & 0xFFFFFFFF;
      }
    }
    if (first_split) {
      first_split = false;
      Unitig left = Unitig(left_unused, 3 - (pos + graph_index.get_k_length() - 2 + (int)(left_unused)) % 4, std::vector < uint8_t > (seq_bases.begin(), seq_bases.begin() + 1 + static_cast < int > (std::floor((pos + graph_index.get_k_length() - 2 + (int)(left_unused)) / 4))));
      graph.insert_unitig(current_id, left); //the left takes the same id to avoid updating the index
      split_to_join[current_id].first = current_id; //unitig having the current id is the left part of the unitig involved in split
    }
    if (next_pos != 0) { // take from pos to next_pos+(k-1)
      new_id++;
      Unitig middle = Unitig((pos + (int)(left_unused)) % 4, 3 - (next_pos + graph_index.get_k_length() - 2 + (int)(left_unused)) % 4, std::vector < uint8_t > (seq_bases.begin() + static_cast < int > (std::floor((pos + (int)(left_unused)) / 4)), seq_bases.begin() + 1 + static_cast < int > (std::floor((next_pos + graph_index.get_k_length() - 2 + (int)(left_unused)) / 4))));
      graph.insert_unitig(new_id, middle);
      if (update_index) {
        graph_index.update_unitig(seq, new_id, current_id, pos, next_pos - 1, true);
      }
    } else { //take from pos to the end
      new_id++;
      Unitig right = Unitig((pos + (int)(left_unused)) % 4, right_unsed, std::vector < uint8_t > (seq_bases.begin() + static_cast < int > (std::floor((pos + (int)(left_unused)) / 4)), seq_bases.end())); //take the right from position till the end
      graph.insert_unitig(new_id, right);
      split_to_join[current_id].second = new_id; //this is the last position of split, so unitig with new_id is the right part of the unitig involved in split
      if (update_index) {
        graph_index.update_unitig(seq, new_id, current_id, pos, seq.unitig_length() - graph_index.get_k_length(), true);
      }
    }
  }
  graph.set_max_node_id(new_id);
}
void insert_funitigs(GfaGraph & graph, std::unordered_map < int, Unitig > & constructed_unitigs) {
  int new_id = graph.get_max_node_id();
  for (auto elem: constructed_unitigs) {
    new_id++;
    graph.insert_unitig(new_id, elem.second); //insert funitig to graph
  }
  graph.set_max_node_id(new_id);
}