#include <chrono>
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
radix_sort(const text_type &source, std::vector<size_type> &sequence_id,
           const std::unique_ptr<std::unordered_map<size_type, size_type>>
               &start_pos_map,
           const std::uint64_t total_nodes) {

  std::chrono::steady_clock::time_point begin, end;
  double init_time = 0;
  begin = std::chrono::steady_clock::now();

  size_type seqs_size = (*start_pos_map).size();
  std::vector<std::vector<std::pair<size_type, node_type>>> sorted_seqs;
  sorted_seqs.reserve(total_nodes);

  thrust::host_vector<size_type> seq_id(sequence_id);
  // thrust::host_vector<size_type> seq_id(seqs_size);
  thrust::host_vector<node_type> keys(seqs_size);
  thrust::device_vector<node_type> d_keys;
  thrust::device_vector<size_type> d_seq_id;
  // std::iota(seq_id.begin(), seq_id.end(), 0);

  // emplace_back the ENDMARKER's outgoing node
  size_type i = 0;
  end = std::chrono::steady_clock::now();
  init_time = std::chrono::duration_cast<std::chrono::microseconds>(end - begin)
                  .count();

  double key_time = 0, h2d_copy_time = 0, sort_time, d2h_copy_time = 0,
         remove_time = 0, place_time = 0;

  for (size_type position = 0; seqs_size > 0; ++position) {
    begin = std::chrono::steady_clock::now();
    i = 0;
    for (auto &idx : seq_id) {
      keys[i++] = source[(*start_pos_map)[idx] + position];
    }
    end = std::chrono::steady_clock::now();
    key_time +=
        std::chrono::duration_cast<std::chrono::microseconds>(end - begin)
            .count();

    // sort sequence_id according to node_id
    begin = std::chrono::steady_clock::now();
    d_keys = keys;
    d_seq_id = seq_id;
    end = std::chrono::steady_clock::now();
    h2d_copy_time +=
        std::chrono::duration_cast<std::chrono::microseconds>(end - begin)
            .count();

    begin = std::chrono::steady_clock::now();
    thrust::stable_sort_by_key(d_keys.begin(), d_keys.end(), d_seq_id.begin());
    end = std::chrono::steady_clock::now();
    sort_time +=
        std::chrono::duration_cast<std::chrono::microseconds>(end - begin)
            .count();

    begin = std::chrono::steady_clock::now();
    keys = d_keys;
    seq_id = d_seq_id;
    end = std::chrono::steady_clock::now();
    d2h_copy_time +=
        std::chrono::duration_cast<std::chrono::microseconds>(end - begin)
            .count();

    // remove paths that reaches the ENDMARKER
    begin = std::chrono::steady_clock::now();
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
    end = std::chrono::steady_clock::now();
    remove_time +=
        std::chrono::duration_cast<std::chrono::microseconds>(end - begin)
            .count();

    begin = std::chrono::steady_clock::now();
    i = 0;
    for (auto &idx : seq_id) {
      node_type next_node_id = source[(*start_pos_map)[idx] + position + 1];
      // sorted_seqs[keys[i] - 1].emplace_back(std::make_pair(idx,
      // next_node_id));
      sorted_seqs[keys[i] - 1].push_back({idx, next_node_id});
      ++i;
    }

    seqs_size = seq_id.size();
    end = std::chrono::steady_clock::now();
    place_time +=
        std::chrono::duration_cast<std::chrono::microseconds>(end - begin)
            .count();
  }
  std::cout << "init_time: " << init_time << " [μs]\n";
  std::cout << "key_time: " << key_time << " [μs]\n";
  std::cout << "h2d_copy_time: " << h2d_copy_time << " [μs]\n";
  std::cout << "sort_time: " << sort_time << " [μs]\n";
  std::cout << "d2h_copy_time: " << d2h_copy_time << " [μs]\n";
  std::cout << "remove_time: " << remove_time << " [μs]\n";
  std::cout << "place_time: " << place_time << " [μs]\n";
  return sorted_seqs;
}

/*
begin = std::chrono::steady_clock::now();
end = std::chrono::steady_clock::now();
auto exe_time =
    std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count();
*/
