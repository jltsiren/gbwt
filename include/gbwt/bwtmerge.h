/*
  Copyright (c) 2018 Jouni Siren
  Copyright (c) 2015 Genome Research Ltd.

  Author: Jouni Siren <jouni.siren@iki.fi>

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.
*/

#ifndef GBWT_BWTMERGE_H
#define GBWT_BWTMERGE_H

#include <condition_variable>
#include <mutex>
#include <thread>

#include "internal.h"

namespace gbwt
{

/*
  bwtmerge.h: Internal support structures adapted from BWT-merge.
*/

//------------------------------------------------------------------------------

/*
  Two-level array allocated in 8-megabyte blocks using mmap().
*/

class BlockArray
{
public:
  typedef gbwt::size_type size_type;
  typedef gbwt::byte_type value_type;

  constexpr static size_type BLOCK_BITS = 23;
  constexpr static size_type BLOCK_SIZE = static_cast<size_type>(1) << BLOCK_BITS;
  constexpr static size_type BLOCK_MASK = BLOCK_SIZE - 1;

  BlockArray();
  BlockArray(const BlockArray& source);
  BlockArray(BlockArray&& source);
  ~BlockArray();

  void swap(BlockArray& source);
  BlockArray& operator=(const BlockArray& source);
  BlockArray& operator=(BlockArray&& source);

  size_type size() const { return this->bytes; }
  size_type blocks() const { return this->data.size(); }
  bool empty() const { return (this->size() == 0); }

  static size_type block(size_type i) { return i >> BLOCK_BITS; }
  static size_type offset(size_type i) { return i & BLOCK_MASK; }

  void clear();

  // Removes the block before block(i).
  void clearUntil(size_type i)
  {
    if(block(i) > 0) { this->clear(block(i) - 1); }
  }

  // Do not use beyond size().
  value_type operator[](size_type i) const
  {
    return this->data[block(i)][offset(i)];
  }

  // Do not use beyond size().
  value_type& operator[](size_type i)
  {
    return this->data[block(i)][offset(i)];
  }

  void push_back(value_type value)
  {
    if(offset(this->bytes) == 0) { this->allocateBlock(); }
    (*this)[this->bytes] = value;
    this->bytes++;
  }

  // The serialized format is compatible with sdsl::int_vector_buffer<8>.
  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  std::vector<value_type*> data;
  size_type                bytes;

private:
  void copy(const BlockArray& source);
  void allocateBlock();
  void clear(size_type block_index);
};  // class BlockArray

//------------------------------------------------------------------------------

/*
  A gap-encoded non-decreasing edge_type array, based on any byte array with
  operator[] and member function push_back(). Intended usage is GapArray<BlockArray>
  in memory and GapArray<sdsl::int_vector_buffer<8>> on disk. Note that the iterator
  is destructive if the array type is BlockArray.
*/

template<class ByteArray>
class GapIterator;

template<class ByteArray>
class GapArray
{
public:
  typedef gbwt::size_type size_type;

  typedef GapIterator<ByteArray> iterator;

  GapArray() { this->value_count = 0; }
  GapArray(const GapArray& source) { this->copy(source); }
  GapArray(GapArray&& source) { *this = std::move(source); }
  ~GapArray() { }

  // Builds a GapArray from the source vector. The vector is sorted during construction.
  // This constructor requires a specialization for initializing the ByteArray.
  explicit GapArray(std::vector<edge_type>& source)
  {
    std::cerr << "GapArray::GapArray(): Unsupported constructor" << std::endl;
  }

  // Merges the input arrays and clears them.
  // This constructor requires a specialization for initializing the ByteArray.
  GapArray(GapArray& a, GapArray& b)
  {
    std::cerr << "GapArray::GapArray(): Unsupported constructor" << std::endl;
  }

  void swap(GapArray& source)
  {
    if(this != &source)
    {
      this->data.swap(source.data);
      std::swap(this->value_count, source.value_count);
    }
  }

  GapArray& operator=(const GapArray& source)
  {
    if(this != &source) { this->copy(source); }
    return *this;
  }

  GapArray& operator=(GapArray&& source)
  {
    if(this != &source)
    {
      this->data = std::move(source.data);
      this->value_count = std::move(source.value_count);
    }
    return *this;
  }

  size_type size() const { return this->value_count; }
  size_type bytes() const { return this->data.size(); }
  bool empty() const { return (this->size() == 0); }

  // Note that there is a specialized clear() for the relevant ByteArray types.
  void clear()
  {
    this->value_count = 0;
  }

  // Store as sdsl::int_vector_buffer<8>. Specialized version exists for BlockArray.
  void write(const std::string& filename)
  {
    sdsl::int_vector_buffer<8> out(filename, std::ios::out);
    for(size_type i = 0; i < this->bytes(); i++) { out.push_back(this->data[i]); }
    out.close();
  }

  // Write the specified node ranges to separate files. The ranges must be contiguous.
  // Calls clear() after the write completes. Returns the total size of the files.
  // This function requires a specialized version for the particular ByteArray.
  size_type write(const std::vector<std::string>& filenames,
             const std::vector<range_type>& node_ranges,
             std::vector<size_type>& value_counts)
  {
    std::cerr << "GapArray::write(): Unsupported ByteArray type" << std::endl;
    this->clear();
    return 0;
  }

  ByteArray data;
  size_type value_count;

  static void writeValue(ByteArray& data, edge_type curr, edge_type& prev)
  {
    ByteCode::write(data, curr.first - prev.first);
    if(curr.first != prev.first) { prev.second = 0; } // Node changed, set previous offset to 0.
    ByteCode::write(data, curr.second - prev.second);
    prev = curr;
  }

private:
  void copy(const GapArray& source)
  {
    this->data = source.data;
    this->value_count = source.value_count;
  }

  void writeValue(edge_type curr, edge_type& prev)
  {
    writeValue(this->data, curr, prev);
  }
};  // class GapArray

template<> GapArray<BlockArray>::GapArray(std::vector<edge_type>& source);
template<> GapArray<BlockArray>::GapArray(GapArray& a, GapArray& b);
template<> void GapArray<BlockArray>::clear();
template<> void GapArray<BlockArray>::write(const std::string& filename);
template<> size_type GapArray<BlockArray>::write(const std::vector<std::string>& filenames,
                                            const std::vector<range_type>& node_ranges,
                                            std::vector<size_type>& value_counts);

void open(GapArray<sdsl::int_vector_buffer<8>>& array, const std::string filename, size_type values);
template<> void GapArray<sdsl::int_vector_buffer<8>>::clear();

//------------------------------------------------------------------------------

template<class ByteArray>
class GapIterator
{
public:
  typedef typename GapArray<ByteArray>::size_type size_type;

  GapIterator() :
    array(nullptr), pos(0), data_pointer(0), value(ENDMARKER, 0)
  {
  }

  GapIterator(GapArray<ByteArray>& data) :
    array(&data), pos(0), data_pointer(0), value(ENDMARKER, 0)
  {
    this->read();
  }

  GapIterator(const GapIterator& source)
  {
    this->copy(source);
  }

  GapIterator& operator=(const GapIterator& source)
  {
    if(this != &source) { this->copy(source); }
    return *this;
  }

  // Iterator operations.
  edge_type operator*() const { return this->value; }
  const edge_type* operator->() const { return &(this->value); }
  void operator++() { this->pos++; this->read(); }
  bool end() const { return (this->pos >= this->array->size()); }

  // Access to the GapArray.
  size_type size() const { return this->array->size(); }
  bool empty() const { return this->array->empty(); }

  // For ProducerBuffer.
  void open() {}
  void close() {}

  GapArray<ByteArray>* array;
  size_type pos, data_pointer;
  edge_type value;

private:
  // There is a specialized destructive read() for BlockArray.
  void read()
  {
    this->readCommon();
  }

  void readCommon()
  {
    if(this->end()) { this->value = invalid_edge(); return; }
    size_type node_diff = ByteCode::read(this->array->data, this->data_pointer);
    if(node_diff != 0) { this->value.second = 0; } // Node changed, set previous offset to 0.
    this->value.first += node_diff;
    this->value.second += ByteCode::read(this->array->data, this->data_pointer);
  }

  void copy(const GapIterator& source)
  {
    this->array = source.array;
    this->pos = source.pos;
    this->data_pointer = source.data_pointer;
    this->value = source.value;
  }
};  // class GapIterator

template<> void GapIterator<BlockArray>::read();

//------------------------------------------------------------------------------

/*
  A producer-consumer buffer. The producer writes to the producer buffer in a background
  thread, while the consumer iterates over the produced data using the consumer buffer.

  Iterator value invalid_edge() is equivalent to end().

  The producer must support the following operations:
  - open(), close()
  - empty()
  - ++, *, end()
*/

template<class Producer>
class ProducerBuffer;

template<class Producer>
void
produce(ProducerBuffer<Producer>* buffer)
{
  while(!(buffer->fill()));
}

// FIXME Maybe this should be a parameter?
constexpr static size_type PRODUCER_BUFFER_SIZE = 8 * MEGABYTE; // Positions.

template<class Producer>
class ProducerBuffer
{
public:
  constexpr static size_type BUFFER_SIZE = PRODUCER_BUFFER_SIZE;

  ProducerBuffer(Producer& source, size_type buffer_size = BUFFER_SIZE) :
    producer(source), buffer_capacity(buffer_size), finished(source.empty()),
    offset(0), value(invalid_edge())
  {
    this->producer_buffer.reserve(this->buffer_capacity);
    this->consumer_buffer.reserve(this->buffer_capacity);
    this->producer.open();
    this->producer_thread = std::thread(produce<Producer>, this);
    this->read();
  }

  ~ProducerBuffer()
  {
    // Stop the producer thread.
    this->mtx.lock();
    this->finished = true;
    this->mtx.unlock();
    if(this->producer_thread.joinable()) { this->producer_thread.join(); }
    // Close the producer.
    this->producer.close();
  }

  // Producer thread.
  Producer&               producer;
  std::vector<edge_type>  producer_buffer;
  size_type               buffer_capacity;
  bool                    finished;

  // Consumer thread.
  std::vector<edge_type>  consumer_buffer;
  size_type               offset;
  edge_type               value;

  // Multithreading.
  std::mutex              mtx;   // For producer data.
  std::condition_variable empty; // Is producer_buffer empty?
  std::thread             producer_thread;

  // Iterator operations.
  edge_type operator*() const { return this->value; }
  const edge_type* operator->() const { return &(this->value); }
  void operator++()
  {
    this->offset++;
    if(this->bufferEnd()) { this->read(); }
    else { this->value = this->consumer_buffer[this->offset]; }
  }
  bool end() { return (this->value == invalid_edge()); }

private:
  bool bufferEnd() const { return (this->offset >= this->consumer_buffer.size()); }

  /*
    Fill producer_buffer and return finished. The second version assumes that the current
    thread holds the mutex and that producer_buffer is empty.
  */
  bool fill()
  {
    // We need the mutex and producer_buffer must be empty.
    std::unique_lock<std::mutex> lock(this->mtx);
    this->empty.wait(lock, [this]() { return producer_buffer.empty(); });
    return this->forceFill();
  }

  bool forceFill()
  {
    if(this->finished) { return true; }
    while(!(this->producer.end()) && this->producer_buffer.size() < this->buffer_capacity)
    {
      this->producer_buffer.push_back(*(this->producer));
      ++(this->producer);
    }
    if(this->producer.end()) { this->finished = true; }
    return this->finished;
  }

  /*
    Swap the buffers and clear producer_buffer. The second version assumes that the current
    thread holds the mutex. Both assume that consumer_buffer is empty.
  */
  void read()
  {
    std::unique_lock<std::mutex> lock(this->mtx);
    this->forceRead();
    this->empty.notify_one();
  }

  void forceRead()
  {
    if(this->producer_buffer.empty())
    {
      this->forceFill();
    }
    this->consumer_buffer.swap(this->producer_buffer);
    this->producer_buffer.clear();
    this->offset = 0;
    this->value = (this->consumer_buffer.empty() ? invalid_edge() : this->consumer_buffer[this->offset]);
  }

  // Producer thread.
  friend void produce<Producer>(ProducerBuffer<Producer>* buffer);

  ProducerBuffer(const ProducerBuffer&) = delete;
  ProducerBuffer& operator=(const ProducerBuffer&) = delete;
};

//------------------------------------------------------------------------------

/*
  A RankArray is a number of temporary files containing GapArrays on disk. Each file is read
  using a separate thread with ProducerBuffer. When the object is deleted, the files are also
  deleted.

  Iterator value invalid_edge() is equivalent to end().
*/

class RankArray
{
public:
  typedef GapArray<sdsl::int_vector_buffer<8>> array_type;
  typedef array_type::size_type                size_type;
  typedef array_type::iterator                 iterator;

#ifdef GBWT_SAVE_MEMORY
  // Avoid comparisons by packing the pair of 32-bit values into a 64-bit integer.
  typedef std::pair<size_type, size_type>      tree_type; // (packed value, source)
#else
  typedef std::pair<edge_type, size_type>      tree_type; // (value, source)
#endif

  const static std::string TEMP_FILE_PREFIX;  // "ranks"

  RankArray();
  ~RankArray();

  /*
    Add a new file and return its name and number. The caller must ensure that only one thread
    calls addFile() at the same time. Do not call when the RankArray is open.
  */
  std::pair<std::string, size_type> addFile();

  /*
    Sets the value count for the last added file. Do not call when the RankArray is open.
  */
  void addValueCount(size_type file, size_type value_count);

  // For ProducerBuffer.
  void open();
  void close();

  // Iterator operations.
  edge_type operator*() const { return this->value; }
  const edge_type* operator->() const { return &(this->value); }
  void operator++();
  bool end() { return (this->value == invalid_edge()); }

  size_type size() const { return this->filenames.size(); }
  size_type empty() const { return (this->size() == 0); }

  std::vector<std::string> filenames;
  std::vector<size_type>   value_counts;

  std::vector<array_type>                inputs;
  std::vector<iterator>                  iterators;
  std::vector<ProducerBuffer<iterator>*> buffers;

  // Use a tournament tree instead of a priority queue.
  // The number of leaves is a power of two.
  std::vector<tree_type> tournament_tree;
  size_type              leaves;

  // Cached value.
  edge_type value;

private:
  // Tournament tree operations.
  void initTree();
  tree_type smaller(size_type tree_offset) const
  {
    if(this->tournament_tree[tree_offset].first <= this->tournament_tree[tree_offset ^ 0x1].first)
    {
      return this->tournament_tree[tree_offset];
    }
    else
    {
      return this->tournament_tree[tree_offset ^ 0x1];
    }
  }

#ifdef GBWT_SAVE_MEMORY
  static size_type pack(edge_type edge) { return (static_cast<size_type>(edge.first) << S_WORD_BITS) | edge.second; }
  static edge_type unpack(size_type packed) { return edge_type(packed >> S_WORD_BITS, packed & LOW_MASK); }
#endif

  void cacheValue()
  {
#ifdef GBWT_SAVE_MEMORY
    this->value = unpack(this->tournament_tree.back().first);
#else
    this->value = this->tournament_tree.back().first;
#endif
  }

  RankArray(const RankArray&) = delete;
  RankArray& operator=(const RankArray&) = delete;
};

//------------------------------------------------------------------------------

/*
  A structure for building RankArray objects. A search thread calls insert() when it wants
  to insert a new element into the buffers.

  1) The element is inserted into pos_buffer. If the buffer is small enough, insert() returns.
  2) We merge pos_buffer with thread_buffer and clear it.
  3) If thread_buffer is small enough, insert() returns.
  4) We merge thread_buffer with the global merge buffers until there is an empty slot
     or we run out of merge buffers.
  5) If there is an empty slot, we insert thread_buffer there. Otherwise we write it to files
     and add the files into the RankArrays. In either case, we clear the thread_buffer.

  Once all elements have been inserted, we need to call flush().
*/

class MergeBuffers
{
public:
  typedef GapArray<BlockArray> buffer_type;

  MergeBuffers(size_type expected_size, size_type num_threads, const MergeParameters& params, const std::vector<range_type>& node_ranges);
  ~MergeBuffers();

  MergeParameters parameters;

  // Thread-specific buffers.
  std::vector<std::vector<edge_type>> pos_buffers;
  std::vector<buffer_type>            thread_buffers;

  // Global merge buffers.
  std::mutex               buffer_lock;
  std::vector<buffer_type> merge_buffers;

  // Rank array.
  std::mutex              ra_lock;
  std::vector<range_type> job_ranges;
  std::vector<RankArray*> ra;
  size_type               ra_values, ra_bytes;
  size_type               final_size;

  std::mutex stderr_access;

  size_type threads() const { return this->pos_buffers.size(); }
  size_type jobs() const { return this->job_ranges.size(); }

  void insert(edge_type element, size_type thread)
  {
    this->pos_buffers[thread].push_back(element);
    if(this->pos_buffers[thread].size() >= this->parameters.posBufferPositions())
    {
      this->insert(this->pos_buffers[thread], this->thread_buffers[thread], false);
    }
  }

  void flush();

private:
  void insert(std::vector<edge_type>& pos_buffer, buffer_type& thread_buffer, bool force_merge);
  void write(buffer_type& buffer);

  MergeBuffers(const MergeBuffers&) = delete;
  MergeBuffers& operator=(const MergeBuffers&) = delete;
};

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_BWTMERGE_H
