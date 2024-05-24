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
    uint64_t query = canonical_bits(((mer << 2) | baseToInt(alphabet[i])) & (0xffffffffffffffff >> (64 - 2 * k)),k);
    if (kmers_to_add_to_graph.count(query) > 0) { //&& !kmers_to_add_to_graph[query]
      possible_character_for_extension.push_back(alphabet[i]);
    } else {
      uint64_t kmer_pos = ind.kmer_position(query); //get the k-mer position(id,position in unitig,orientation)
      if (kmer_pos != 0) { //we found an mphf corresponding to the minimizer of the k-mer
        int may_exist = graph.test_kmer_presence(query, mer, kmer_pos, ind.get_k_length());
        if (may_exist > 0) { //we encountered a split position we add the position to a vector which will be sorted and processed
          kmer_pos = kmer_pos >> 32;
          kmer_pos = (kmer_pos << 32) | may_exist;
          candidate_splits[kmer_pos] = true;
          return std::vector < char > ();
        } else if (may_exist > -3 && number_join < 2) { //the k-mer is in the graph and we encountered a possible join
          /*
          add some conditions to check if it is a join
          how many possible joins are found=> if 1 then maybe join => let from left extention decide
          when the value of may_exist is 0 then we encounter a suffix prefix overlap
          which may result in a join case
          we stop at the second occurence of suffix-prefix overlap because it is no more a join
          the conditions outside the loop will correct the value of join_position
          */
         if(may_exist==0 || may_exist==-1){
          join_position = kmer_pos;
         }
          number_join++;
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
          uint64_t k_mer_from_funitigs = canonical_bits(((suffix << 2) | baseToInt(character_to_add) & (0xffffffffffffffff >> (64 - 2 * unitigs_index.get_k_length()))), unitigs_index.get_k_length());
          if (kmers_to_add_to_graph.count(k_mer_from_funitigs)>0 && !kmers_to_add_to_graph[k_mer_from_funitigs]) {
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
        std::cout<<id_u<<","<<position_unitig<<","<<id<<" funId\n";
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
  std::unordered_map < int, std::pair < int, int >> & candidate_join, std::unordered_map < int, std::pair < int, int >> & funitigs_join) {
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
void update_CdBG(Index_mphf & graph_index, GfaGraph & graph, std::unordered_map < uint64_t, bool > & candidate_splits, 
                std::unordered_map < int, std::pair < int, int >> & candidate_join, std::unordered_map<int,int> &map_on_split_if_join,bool update_index) {
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
      if(candidate_join.count(current_id) > 0){
        map_on_split_if_join[current_id] = new_id;
      }
      if (update_index) {
        graph_index.update_unitig(seq, new_id, current_id, pos, seq.unitig_length() - graph_index.get_k_length(), true);
      }
    }
  }
  graph.set_max_node_id(new_id);
}
void insert_funitigs(GfaGraph & graph, std::unordered_map < int, Unitig > & constructed_unitigs, std::vector< int > & unused_ids) {
  std::vector< int > skipped_ids(unused_ids.size());
  int i = 0;
  int new_id = 0;
  // use the ids of removed unitigs (those unitigs that are merged with other unitigs)
  for(auto it = constructed_unitigs.begin(); it != constructed_unitigs.end() && i < unused_ids.size(); it++){
    new_id = unused_ids[i];
    skipped_ids[i] = it->first;
    i++;
    graph.insert_unitig(new_id, it->second);
  }
  //delte the funitigs that took unused id (delete them from funitigs)
  for(int elem : skipped_ids){
    constructed_unitigs.erase(elem);
  }
  //insert the remaining funitigs
  new_id = graph.get_max_node_id();
  for (auto elem: constructed_unitigs) {
    new_id++;
    graph.insert_unitig(new_id, elem.second); //insert funitig to graph
  }
  graph.set_max_node_id(new_id);
}
void merge_all(GfaGraph& graph,std::unordered_map<int,Unitig>& funitigs,Index_mphf &unitig_index,std::unordered_map<int,std::pair<int,int>> &funitig_joins,
            std::unordered_map<int,std::pair<int,int>> &unitig_joins,std::unordered_map<int,int> &right_after_split,
            std::unordered_map<int,std::vector<std::pair<int,int>>> &track_funitig_offsets,std::vector<int> &unused_id){
    std::unordered_map<uint64_t,uint64_t> suffix_prefix_funitigs;
    std::unordered_map<uint64_t,uint64_t> suffix_prefix_unitigs;
    index_concerned_funitigs(suffix_prefix_funitigs,funitigs,funitig_joins,unitig_index.get_k_length());
    index_concerned_unitigs(suffix_prefix_unitigs,graph,unitig_joins,right_after_split,unitig_index.get_k_length());

    for(auto &elem:suffix_prefix_funitigs){
      merge_one_funitig(elem.first,elem.second,suffix_prefix_funitigs,suffix_prefix_unitigs,funitigs,graph,unitig_index,track_funitig_offsets,unused_id);
    }
}
void index_concerned_funitigs(std::unordered_map<uint64_t,uint64_t> &suffix_prefix_funitigs,std::unordered_map<int,Unitig>& funitigs,std::unordered_map<int,std::pair<int,int>> &funitig_joins,int k){
  for(auto &elem:funitig_joins){
    Unitig u=funitigs[elem.first];
    index_suffix_prefix(elem.first,elem.second.first,elem.second.second,u,suffix_prefix_funitigs,k);
  }
}
void index_concerned_unitigs(std::unordered_map<uint64_t,uint64_t> &suffix_prefix_unitigs,GfaGraph& graph,std::unordered_map<int,std::pair<int,int>> &unitig_joins,
                            std::unordered_map<int,int> right_after_split,int k){
    for(auto &elem:unitig_joins){
      Unitig u=graph.get_unitig(elem.first);
      if(elem.second.first != 0){
        insert_suff_pref(u.get_ith_mer(0,k-1),elem.first,0,k,suffix_prefix_unitigs);
      }
      if(elem.second.second != 0){
        if(right_after_split.count(elem.first) > 0){
          u = graph.get_unitig(right_after_split[elem.first]);
        }
        insert_suff_pref(u.get_ith_mer(u.unitig_length() - k + 1,k-1),elem.first,u.unitig_length() - k + 1,k,suffix_prefix_unitigs);
      }
  }
}
void index_suffix_prefix(const int id,const int left,const int right,Unitig &u,std::unordered_map<uint64_t,uint64_t> &suffix_prefix_funitigs,int k){
      if(left!=0){
        insert_suff_pref(u.get_ith_mer(0,k-1),id,0,k,suffix_prefix_funitigs);
      }
      if(right!=0){
        insert_suff_pref(u.get_ith_mer(u.unitig_length() - k + 1,k-1),id,u.unitig_length() - k + 1,k,suffix_prefix_funitigs);
      }
}
void merge_one_funitig(uint64_t kmmer,uint64_t position,std::unordered_map<uint64_t,uint64_t> &suffix_prefix_funitigs,
          std::unordered_map<uint64_t,uint64_t> &suffix_prefix_unitigs,std::unordered_map<int,Unitig>& funitigs,GfaGraph& graph,Index_mphf &unitig_index,
          std::unordered_map<int,std::vector<std::pair<int,int>>> &track_funitig_offsets,std::vector<int> &unused_id){
    uint64_t pos_unitig=suffix_prefix_unitigs[kmmer];
    int funitig_id = position >> 32;
    int kmm_offset_funitig= position & 0x7FFFFFFF;
    int unitig_id = pos_unitig >> 32;
    int kmm_offset_unitig= pos_unitig & 0x7FFFFFFF;
    Unitig funitig = funitigs[funitig_id];
    Unitig uni = graph.get_unitig(unitig_id);
    bool order1 = false;
    bool reversed1 = false;
    std::string merged = concat(funitig,uni,kmm_offset_funitig,kmm_offset_unitig,unitig_index.get_k_length(),order1,reversed1);
    int kmm_offset = abs(kmm_offset_funitig - funitig.unitig_length() + unitig_index.get_k_length() -1);
    uint64_t kmm_extremity = funitig.get_ith_mer(kmm_offset,unitig_index.get_k_length() - 1);
    kmm_extremity = canonical_bits(kmm_extremity,unitig_index.get_k_length() - 1);
    if(reversed1 && track_funitig_offsets.count(unitig_id) > 0){
      track_funitig_offsets[unitig_id]=reverseOffsets(track_funitig_offsets[unitig_id],uni.unitig_length()-1);
    }
    if(suffix_prefix_funitigs.count(kmm_extremity) > 0){
      uint64_t position_ext_uni = suffix_prefix_unitigs[kmm_extremity];
      int unitig_ext_id = position_ext_uni >> 32;
      int kmm_offset_ext_uni = position_ext_uni & 0x7FFFFFFF;
      bool order2=false;
      bool reversed2= false;
      Unitig uni2 = graph.get_unitig(unitig_ext_id);
      std::string merged2 = concat(funitig,uni2,kmm_offset_funitig,kmm_offset_ext_uni,unitig_index.get_k_length(),order2,reversed2);
      if(reversed2 && track_funitig_offsets.count(unitig_ext_id) > 0){
        track_funitig_offsets[unitig_ext_id]=reverseOffsets(track_funitig_offsets[unitig_ext_id],uni2.unitig_length()-1);
      }
      if(order2){
        merge_and_track(merged,merged2,unitig_ext_id,unitig_id,uni.unitig_length(),uni2.unitig_length(),funitig.unitig_length(),unitig_index,graph,track_funitig_offsets,unused_id,suffix_prefix_unitigs);
      }
      else{
        merge_and_track(merged2,merged,unitig_id,unitig_ext_id,uni2.unitig_length(),uni.unitig_length(),funitig.unitig_length(),unitig_index,graph,track_funitig_offsets,unused_id,suffix_prefix_unitigs);
      }
    }
    else{
      std::vector<std::pair<int,int>> offsets;
      offsets.push_back(std::pair<int,int>(0,funitig.unitig_length()-1));
      std::vector<std::pair<int,int>> offset_in_unitig;
      Unitig res = Unitig(merged);
      if(track_funitig_offsets.count(unitig_id) > 0){
        offset_in_unitig = track_funitig_offsets[unitig_id];
      }
      if(order1){
        shift_and_merge(offset_in_unitig,offsets,uni.unitig_length() - unitig_index.get_k_length() + 1);
        track_funitig_offsets[unitig_id] = offset_in_unitig;
        update_suff_pref_index(res.get_ith_mer(0,unitig_index.get_k_length()-1),unitig_id,0,unitig_index.get_k_length(),suffix_prefix_unitigs);
      }
      else{
        shift_and_merge(offsets,offset_in_unitig,funitig.unitig_length() - unitig_index.get_k_length() + 1);
        track_funitig_offsets[unitig_id] = offsets;
        update_suff_pref_index(res.get_ith_mer(res.unitig_length() - unitig_index.get_k_length() + 1,unitig_index.get_k_length()-1),unitig_id,res.unitig_length() - unitig_index.get_k_length() + 1,unitig_index.get_k_length(),suffix_prefix_unitigs);
      }
      graph.insert_unitig(unitig_id,res);
    }
    funitigs.erase(funitig_id);
}
std::string concat(Unitig funitig, Unitig uni,int offset_funitig,int offset_unitig,int k,bool &order, bool &reversed){
  std::string result("");
  if(offset_funitig == 0 && offset_unitig > 0){
    result = uni.to_string() + funitig.to_string().substr(k-1);
    order = true;
  }
  else if (offset_funitig > 0 && offset_unitig == 0)
  {
    result = funitig.to_string() + uni.to_string().substr(k-1);
    order = false;
  }
  else if (offset_funitig > 0 && offset_unitig > 0)
  {
    order = false;
    reversed =true;
    result = funitig.to_string() + reverseComplement(uni.to_string()).substr(k-1);
  }
  else{
    order = true;
    reversed = true;
    result = reverseComplement(uni.to_string()) + funitig.to_string().substr(k-1);
  }
  return result;
}
void shift_and_merge(std::vector<std::pair<int, int>> &destination, const std::vector<std::pair<int, int>> &source, int shift) {
    // shift elements in source
    destination.reserve(destination.size()+source.size());
    auto shift_pair = [shift](const std::pair<int, int>& p) {
        return std::make_pair(p.first + shift, p.second + shift);
    }; 
    // append elements in source to destination
    std::transform(source.begin(), source.end(), std::back_inserter(destination), shift_pair);
}
void merge_and_track(std::string merged,std::string merged2,int unitig1_id,int unitig2_id,
                      int unitig1_length,int unitig2_length,int funitig_length,
                      Index_mphf &unitig_index,GfaGraph& graph,std::unordered_map<int,std::vector<std::pair<int,int>>> &track_funitig_offsets,
                      std::vector<int> &unused_id,std::unordered_map<uint64_t,uint64_t> &suffix_prefix_unitigs){
  std::vector<std::pair<int,int>> offsets;
  if(track_funitig_offsets.count(unitig1_id)){
    offsets = track_funitig_offsets[unitig1_id];
  }
  offsets.push_back(std::pair<int,int>(unitig2_length - unitig_index.get_k_length() + 1,unitig2_length - unitig_index.get_k_length() + funitig_length));
  merged2 = merged2 + merged.substr(funitig_length);
  Unitig res = Unitig(merged2);
  graph.delete_unitig(unitig2_id);
  graph.insert_unitig(unitig1_id,res);
  shift_and_merge(offsets,track_funitig_offsets[unitig2_id],unitig1_length);
  track_funitig_offsets[unitig1_id]=offsets;
  unused_id.push_back(unitig2_id);
  update_suff_pref_index(res.get_ith_mer(0,unitig_index.get_k_length()-1),unitig1_id,0,unitig_index.get_k_length(),suffix_prefix_unitigs);
  update_suff_pref_index(res.get_ith_mer(res.unitig_length()-unitig_index.get_k_length()+1,unitig_index.get_k_length()-1),unitig1_id,res.unitig_length()-unitig_index.get_k_length()+1,unitig_index.get_k_length(),suffix_prefix_unitigs);
}
void update_suff_pref_index(uint64_t kmmer,int id,int offset,int k,std::unordered_map<uint64_t,uint64_t> &suffix_prefix_unitigs){
  std::tuple<uint64_t,bool> canon_kmmer = reverseComplementCanonical(kmmer,k-1);
  if(suffix_prefix_unitigs.count(std::get<0>(canon_kmmer)) > 0){
      uint64_t position=id;
      position=(position<<31)|offset;
      position=(position<<1)|std::get<1>(canon_kmmer);
      suffix_prefix_unitigs[std::get<0>(canon_kmmer)]=position;
  }
}
void insert_suff_pref(uint64_t kmmer,int id,int offset,int k,std::unordered_map<uint64_t,uint64_t> &suffix_prefix_index){
  std::tuple<uint64_t,bool> canon_kmmer = reverseComplementCanonical(kmmer,k-1);
  uint64_t position=id;
  position=(position<<31)|offset;
  position=(position<<1)|std::get<1>(canon_kmmer);
  suffix_prefix_index[std::get<0>(canon_kmmer)]=position;
}
void update_graph(GfaGraph &graph, Index_mphf &unitig_index,std::unordered_map<int,Unitig> &funitigs, std::unordered_map<int , std::pair<int, int>> &unitig_joins,
                  std::unordered_map<int , std::pair<int, int>> &funitig_joins,std::unordered_map < uint64_t, bool > & candidate_splits, bool update_index){
    std::unordered_map<int,int> map_on_split_if_join;// remember the id of the right part of a unitig split of it can join a funitig
    update_CdBG(unitig_index,graph,candidate_splits,unitig_joins,map_on_split_if_join,update_index);
    std::vector<int> unused_ids;
    std::unordered_map< int, std::vector<std::pair< int, int >>> track_funitig_offsets;
    //merge(graph,funitigs,unitig_index,funitig_joins,unitig_joins,,map_on_split_if_join);
    merge_all(graph,funitigs,unitig_index,funitig_joins,unitig_joins,map_on_split_if_join,track_funitig_offsets,unused_ids);
    if(update_index){
      update_all_with_skips(track_funitig_offsets, graph, unitig_index);
      unitig_index.update_index(funitigs, graph, unused_ids, track_funitig_offsets);
    }
    else{
      insert_funitigs(graph, funitigs, unused_ids);
    }
}
void update_all_with_skips(const std::unordered_map< int, std::vector<std::pair< int, int >>> &track_funitig_offsets, GfaGraph &graph, Index_mphf &unitig_index){
  for(auto &skips_offsets: track_funitig_offsets){
    unitig_index.update_unitig_with_skips(graph.get_unitig(skips_offsets.first),skips_offsets.first,skips_offsets.first,skips_offsets.second);
  }
}
std::vector<std::pair<int,int>> reverseOffsets(std::vector<std::pair<int,int>> funiOffsets,int shift){
    std::vector<std::pair<int,int>> reversed(funiOffsets.size());
    for(int i = funiOffsets.size()-1;i>=0;i--){
        reversed[funiOffsets.size()-i-1] = std::pair<int,int>(shift-funiOffsets[i].second,shift-funiOffsets[i].first);
    }
    return reversed;
}