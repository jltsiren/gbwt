#include <chrono>
#include <iostream>
#include <thrust/device_free.h>
#include <thrust/device_malloc.h>
#include <thrust/device_ptr.h>
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

__device__ uint64_t get_int(const uint64_t *m_data, const size_type idx,
                            const uint8_t len) {
  const uint64_t *word = m_data + (idx >> 6);
  const uint8_t offset = idx & 0x3F;
  uint64_t w1 = (*word) >> offset;
  if ((offset + len) > 64) { // if offset+len > 64
                             // w1 or w2 adepted:
    return w1 | ((*(word + 1) &
                  (1 << ((offset + len) & 0x3F)) - 1) // set higher bits zero
                 << (64 - offset));                   // move bits to the left
  } else {
    return w1 & ((1 << len) - 1);
  }
}

__global__ void test_get_int(const uint64_t *d_source, node_type *d_test,
                             const size_type vec_size,
                             const size_type start_idx, const int t_width) {
  size_type idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < vec_size) {
    d_test[idx] =
        get_int(d_source, size_type((start_idx + idx)) * t_width, t_width);
  }
  return;
}

__global__ void assign_key(const uint64_t *d_source,
                           const size_type *d_start_pos, node_type *d_keys,
                           size_type *d_seq_id, const size_type position,
                           const size_type first_seq_id,
                           const size_type seqs_size, const int t_width) {
  size_type idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < seqs_size) {
    d_keys[idx] = get_int(
        d_source,
        (d_start_pos[d_seq_id[idx] - first_seq_id] + position) * t_width,
        t_width);
  }
  return;
}

struct is_zero {
  __host__ __device__ bool operator()(int x) { return (x == 0); }
};

void 
radix_sort(const text_type &source, std::vector<size_type> &sequence_id,
           const std::unique_ptr<std::unordered_map<size_type, size_type>>
               &start_pos_map,
           std::vector<std::vector<std::pair<size_type, node_type>>> &sorted_seqs,
           const std::uint64_t total_nodes) {
  // ---- Prepare values to be used ---- //
  // width of the integers which are accessed via the [] operator
  const uint8_t t_width = source.width();
  const size_type source_size_byte = source.bit_size() / 8 + 1;
  const size_type seqs_size = (*start_pos_map).size(); // total sequence size
  size_type first_seq_id = sequence_id[0];

  // copy source to device memory
  uint64_t *d_source;
  cudaMalloc(&d_source, source_size_byte);
  cudaMemcpyAsync(d_source, source.data(), source_size_byte,
                  cudaMemcpyHostToDevice);

  // copy start_position
  thrust::host_vector<size_type> start_pos;
  start_pos.reserve(seqs_size);
  for (auto &seq_id : sequence_id) {
    start_pos.push_back((*start_pos_map)[seq_id]);
  }
  thrust::device_vector<size_type> start_pos_vec = start_pos;
  size_type *d_start_pos = thrust::raw_pointer_cast(&start_pos_vec[0]);

  
  // copy sequence_id
  thrust::host_vector<size_type> h_seq_id(sequence_id);
  thrust::device_vector<size_type> seq_id_vec = h_seq_id;
  thrust::device_ptr<size_type> d_seq_id =
      thrust::device_pointer_cast(&seq_id_vec[0]);
  size_type *d_seq_id_raw = thrust::raw_pointer_cast(&seq_id_vec[0]);
  
  // use gpu to assign keys
  // there about 2000~ sequences passed in
  thrust::host_vector<node_type> h_keys_vec(seqs_size);
  thrust::device_ptr<node_type> d_keys =
      thrust::device_malloc(sizeof(node_type) * seqs_size);
  node_type *d_keys_raw = thrust::raw_pointer_cast(d_keys);

  size_type arr_start_idx = 0;     // first index which is not an ENDMARKER
  size_type seqs_left = seqs_size; // sequences that have not
                                   // reached the ENDMARKER
  const int thread_per_block = 512;
  cudaDeviceSynchronize();
  

  for (size_type position = 0; seqs_left > 0; ++position) {
    //---- Assign Keys ----//
    const int block_per_grid =
        (seqs_left + thread_per_block - 1) / thread_per_block;
    assign_key<<<block_per_grid, thread_per_block>>>(
        d_source, d_start_pos, d_keys_raw + arr_start_idx,
        d_seq_id_raw + arr_start_idx, position, first_seq_id, seqs_left,
        t_width);

    //---- Radix Sort ----//
    thrust::stable_sort_by_key(d_keys + arr_start_idx, d_keys + seqs_size,
                               d_seq_id + arr_start_idx);

    //---- Copy keys and sequence id back to host ----//
    thrust::copy(d_seq_id + arr_start_idx, d_seq_id + seqs_size,
                 h_seq_id.begin());
    thrust::copy(d_keys + arr_start_idx, d_keys + seqs_size,
                 h_keys_vec.begin());

    //---- Remove paths that reaches the ENDMARKER(zero):  version2 ----//
    // remove: find the first index in h_keys_vec that is not an ENDMARKER
    ///*
    size_type end_counter = 0;
    for (size_type i = 0; i < seqs_left; ++i) {
      if (h_keys_vec[i] != gbwt::ENDMARKER) {
        break;
      } else {
        ++end_counter;
      }
    }
    arr_start_idx += end_counter;
    seqs_left = seqs_size - arr_start_idx;
    if (seqs_left <= 0)
      break;
    //*/
    for (size_type i = end_counter; i < seqs_left+end_counter; ++i) {
      size_type seq_id = h_seq_id[i];
      node_type next_node_id = source[(*start_pos_map)[seq_id] + position + 1];
      sorted_seqs[h_keys_vec[i] - 1].push_back({seq_id, next_node_id});
    }
  }
  cudaFree(d_source);
  thrust::device_free(d_keys);
}