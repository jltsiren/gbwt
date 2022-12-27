#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <thrust/sort.h>
#include <thrust_sort.cuh>

void print_vec(const thrust::host_vector<size_type> &vec) {
  for (auto &item : vec) {
    std::cout << item << " ";
  }
  std::cout << "\n";
}

std::vector<std::vector<std::pair<size_type, node_type>>>
radix_sort(const text_type &source,
           const std::vector<size_type> &start_position,
           const std::uint64_t total_nodes) {
  std::vector<std::vector<std::pair<size_type, node_type>>> sorted_seqs(
      total_nodes);
  size_type seqs_size = start_position.size();
  thrust::host_vector<size_type> seq_id(seqs_size);
  thrust::host_vector<node_type> keys(seqs_size);
  std::iota(seq_id.begin(), seq_id.end(), 0);

  // emplace_back the ENDMARKER's outgoing node
  size_type i = 0;
  //for (auto &first_node_pos : start_position) {
  //  sorted_seqs[0].emplace_back(std::make_pair(i++, source[first_node_pos]));
  //}

  for (size_type position = 0; seqs_size > 0; ++position) {
    i = 0;
    for (auto &idx : seq_id) {
      keys[i++] = source[start_position[idx] + position];
    }

    // sort sequence_id according to node_id
    thrust::device_vector<node_type> d_keys = keys;
    thrust::device_vector<size_type> d_seq_id = seq_id;
    thrust::stable_sort_by_key(d_keys.begin(), d_keys.end(), d_seq_id.begin());
    seq_id = d_seq_id;
    keys = d_keys;

    // remove paths that reaches the ENDMARKER
    size_type end_counter = 0;
    for (auto &key : keys) {
      if (key != gbwt::ENDMARKER) {
        break;
      } else {
        ++end_counter;
      }
    }
    if (end_counter > 0) {
      keys.assign(keys.begin() + end_counter, keys.end());
      seq_id.assign(seq_id.begin() + end_counter, seq_id.end());
    }

    i = 0;
    for (auto &idx : seq_id) {
      node_type next_node_id = source[start_position[idx] + position + 1];
      sorted_seqs[keys[i++] - 1].emplace_back(std::make_pair(idx, next_node_id));
    }

    seqs_size = seq_id.size();
  }
  return sorted_seqs;
}
