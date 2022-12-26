#pragma once
#include <gbwt/utils.h>
//#include <thrust/sort.h>

using gbwt::node_type;
using gbwt::range_type;
using gbwt::size_type;
using gbwt::text_type;

std::vector<std::vector<std::pair<size_type, node_type>>>
radix_sort(const text_type &source,
           const std::vector<size_type> &start_position,
           const std::uint64_t total_nodes);
