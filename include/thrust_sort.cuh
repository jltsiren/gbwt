#pragma once
#include <gbwt/utils.h>
//#include <thrust/sort.h>

using gbwt::node_type;
using gbwt::range_type;
using gbwt::size_type;
using gbwt::text_type;

void 
radix_sort(const text_type &source,
          std::vector<size_type> &sequence_id,
          const std::unique_ptr<std::unordered_map<size_type, size_type>> &start_pos_map,
          std::vector<std::vector<std::pair<size_type, node_type>>> &sorted_seqs,
          const std::uint64_t total_nodes);
