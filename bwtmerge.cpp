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

#include <gbwt/bwtmerge.h>

#include <cstring>
#include <sys/mman.h>

namespace gbwt
{

//------------------------------------------------------------------------------

BlockArray::BlockArray() :
  bytes(0)
{
}

BlockArray::BlockArray(const BlockArray& source)
{
  this->copy(source);
}

BlockArray::BlockArray(BlockArray&& source)
{
  *this = std::move(source);
}

BlockArray::~BlockArray()
{
  this->clear();
}

void
BlockArray::copy(const BlockArray& source)
{
  this->clear();

  this->bytes = source.bytes;
  this->data.reserve(source.data.size());
  for(size_type i = 0; i < source.data.size(); i++)
  {
    if(source.data[i] == nullptr) { this->data.push_back(nullptr); }
    else
    {
      this->allocateBlock();
      std::memcpy(static_cast<void*>(this->data[i]), static_cast<void*>(source.data[i]), BLOCK_SIZE);
    }
  }
}

void
BlockArray::swap(BlockArray& source)
{
  if(this != &source)
  {
    this->data.swap(source.data);
    std::swap(this->bytes, source.bytes);
  }
}

BlockArray&
BlockArray::operator=(const BlockArray& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

BlockArray&
BlockArray::operator=(BlockArray&& source)
{
  if(this != &source)
  {
    this->clear();
    this->swap(source);
  }
  return *this;
}

BlockArray::size_type
BlockArray::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  // Header.
  size_type bits = this->size() * BYTE_BITS;
  written_bytes += sdsl::write_member(bits, out, child, "bits");

  // Data.
  size_type data_bytes = 0;
  sdsl::structure_tree_node* data_node = sdsl::structure_tree::add_child(child, "data", "gbwt::byte_type*");
  for(size_type i = 0; i < this->blocks(); i++)
  {
    size_type block_bytes = BLOCK_SIZE;
    if(i + 1 == this->blocks())
    {
      block_bytes = this->size() - data_bytes;
      block_bytes += sizeof(std::uint64_t) - this->size() % sizeof(std::uint64_t);
    }
    DiskIO::write(out, this->data[i], block_bytes);
    data_bytes += block_bytes;
  }
  sdsl::structure_tree::add_size(data_node, data_bytes);
  written_bytes += data_bytes;

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
BlockArray::load(std::istream& in)
{
  this->clear();

  // Header.
  size_type bits = 0;
  sdsl::read_member(bits, in);
  this->bytes = bits / BYTE_BITS;

  // Data.
  size_type data_bytes = 0;
  size_type total_blocks = (this->size() + BLOCK_SIZE - 1) / BLOCK_SIZE;
  this->data.reserve(total_blocks);
  for(size_type i = 0; i < total_blocks; i++)
  {
    this->allocateBlock();
    size_type block_bytes = BLOCK_SIZE;
    if(i + 1 == total_blocks)
    {
      block_bytes = this->size() - data_bytes;
      block_bytes += sizeof(std::uint64_t) - this->size() % sizeof(std::uint64_t);
    }
    DiskIO::read(in, this->data[i], block_bytes);
    data_bytes += block_bytes;
  }
}

void
BlockArray::clear()
{
  for(size_type i = 0; i < this->blocks(); i++)
  {
    this->clear(i);
  }
  this->data.clear();
  this->bytes = 0;
}

void
BlockArray::allocateBlock()
{
  value_type* ptr = static_cast<value_type*>(mmap(0, BLOCK_SIZE, PROT_READ | PROT_WRITE, MAP_ANON | MAP_PRIVATE, -1, 0));
  std::memset(ptr, 0, BLOCK_SIZE);
  this->data.push_back(ptr);
}

void
BlockArray::clear(size_type block_index)
{
  if(this->data[block_index] == nullptr) { return; }
  munmap(static_cast<void*>(this->data[block_index]), BLOCK_SIZE);
  this->data[block_index] = nullptr;
}

//------------------------------------------------------------------------------

template<>
GapArray<BlockArray>::GapArray(std::vector<edge_type>& source)
{
  this->value_count = source.size();
  if(source.empty()) { return; }

  sequentialSort(source.begin(), source.end());
  edge_type prev(ENDMARKER, 0);
  for(edge_type value : source) { this->writeValue(value, prev); }
}

template<>
GapArray<BlockArray>::GapArray(GapArray& a, GapArray& b)
{
  this->value_count = 0;
  if(a.empty()) { this->swap(b); return; }
  if(b.empty()) { this->swap(a); return; }

  iterator a_iter(a), b_iter(b);
  edge_type prev(ENDMARKER, 0);
  while(!(a_iter.end()) || !(b_iter.end()))
  {
    edge_type curr;
    if(*a_iter <= *b_iter) { curr = *a_iter; ++a_iter; }
    else { curr = *b_iter; ++b_iter; }
    this->writeValue(curr, prev);
    this->value_count++;
  }

  a.clear(); b.clear();
}

template<>
void
GapArray<BlockArray>::clear()
{
  this->data.clear();
  this->value_count = 0;
}

template<>
void
GapArray<BlockArray>::write(const std::string& filename)
{
  std::ofstream out(filename, std::ios_base::binary);
  if(!out)
  {
    std::cerr << "GapArray::write(): Cannot open output file " << filename << std::endl;
    std::exit(EXIT_FAILURE);
  }
  this->data.serialize(out);
  out.close();
}


void
open(GapArray<sdsl::int_vector_buffer<8>>& array, const std::string filename, size_type values)
{
  array.data = sdsl::int_vector_buffer<8>(filename);
  array.value_count = values;
}

template<>
void
GapArray<sdsl::int_vector_buffer<8>>::clear()
{
  this->data.close();
  this->value_count = 0;
}

//------------------------------------------------------------------------------

template<>
void
GapIterator<BlockArray>::read()
{
  this->readCommon();
  this->array->data.clearUntil(this->data_pointer);
}

//------------------------------------------------------------------------------

const std::string RankArray::TEMP_FILE_PREFIX = "ranks";

RankArray::RankArray()
{
}

RankArray::~RankArray()
{
  this->close();
  for(size_type i = 0; i < this->size(); i++) { TempFile::remove(this->filenames[i]); }
}

void
RankArray::addFile(GapArray<BlockArray>& array)
{
  std::string filename = TempFile::getName(TEMP_FILE_PREFIX);
  this->filenames.push_back(filename);
  this->value_counts.push_back(array.size());
  array.write(filename);
  array.clear();
}

void
RankArray::open()
{
  this->close();
  this->inputs = std::vector<array_type>(this->size());
  this->iterators = std::vector<iterator>(this->size());

  for(size_type i = 0; i < this->size(); i++)
  {
    gbwt::open(this->inputs[i], this->filenames[i], this->value_counts[i]);
    this->iterators[i] = iterator(this->inputs[i]);
  }

  this->heapify();
}

void
RankArray::close()
{
  this->iterators.clear();
  this->inputs.clear();
}

void
RankArray::heapify()
{
  if(this->size() <= 1) { return; }

  size_type i = parent(this->size() - 1);
  while(true)
  {
    this->down(i);
    if(i == 0) { break; }
    i--;
  }
}

//------------------------------------------------------------------------------

void
produce(RankArrayBuffer* buffer)
{
  while(!(buffer->fill()));
}

RankArrayBuffer::RankArrayBuffer(RankArray& data) :
  array(data), finished(data.empty()),
  offset(0)
{
  this->producer_buffer.reserve(BUFFER_SIZE);
  this->consumer_buffer.reserve(BUFFER_SIZE);

  this->array.open();
  this->producer_thread = std::thread(produce, this);
  this->read();
}

RankArrayBuffer::~RankArrayBuffer()
{
  // Stop the producer thread.
  this->mtx.lock();
  this->finished = true;
  this->mtx.unlock();
  if(this->producer_thread.joinable()) { this->producer_thread.join(); }

  this->array.close();
}

bool
RankArrayBuffer::producerEnd()
{
  std::unique_lock<std::mutex> lock(this->mtx);
  return (this->finished && this->producer_buffer.empty());
}

bool
RankArrayBuffer::fill()
{
  // We need the mutex and producer_buffer must be empty.
  std::unique_lock<std::mutex> lock(this->mtx);
  this->empty.wait(lock, [this]() { return producer_buffer.empty(); });

  return this->forceFill();
}

bool
RankArrayBuffer::forceFill()
{
  if(this->finished) { return true; }

  while(this->producer_buffer.size() < BUFFER_SIZE && !(this->array.end()))
  {
    this->producer_buffer.push_back(*(this->array)); ++(this->array);
  }

  if(this->array.end()) { this->finished = true; }
  return this->finished;
}

void
RankArrayBuffer::read()
{
  std::unique_lock<std::mutex> lock(this->mtx);
  this->forceRead();
  this->empty.notify_one();
}

void
RankArrayBuffer::forceRead()
{
  if(this->producer_buffer.empty())
  {
    this->forceFill();
  }
  this->consumer_buffer.swap(this->producer_buffer);
  this->producer_buffer.clear();
}

//------------------------------------------------------------------------------

MergeBuffers::MergeBuffers(size_type expected_size, const MergeParameters& params) :
  parameters(params),
  merge_buffers(params.merge_buffers),
  ra_values(0), ra_bytes(0), final_size(expected_size)
{
}

MergeBuffers::~MergeBuffers()
{
}

void
MergeBuffers::insert(std::vector<edge_type>& pos_buffer, buffer_type& thread_buffer, bool force_merge)
{
  // Compress pos_buffer.
  buffer_type temp_buffer(pos_buffer); pos_buffer.clear();

  // Merge it into thread_buffer.
  thread_buffer = buffer_type(thread_buffer, temp_buffer);
  if(!force_merge && thread_buffer.bytes() < this->parameters.thread_buffer_size) { return; }

  if(Verbosity::level >= Verbosity::FULL)
  {
    std::lock_guard<std::mutex> lock(this->stderr_access);
    std::cerr << "MergeBuffers::insert(): Thread " << omp_get_thread_num() << ": Inserting "
              << thread_buffer.size() << " values to the merge buffers" << std::endl;
  }

  // Merge with the existing merge buffers until we find an empty slot.
  for(size_type i = 0; i < this->merge_buffers.size(); i++)
  {
    {
      std::lock_guard<std::mutex> lock(this->buffer_lock);
      if(this->merge_buffers[i].empty())
      {
        thread_buffer.swap(this->merge_buffers[i]);
        if(Verbosity::level >= Verbosity::FULL)
        {
          std::lock_guard<std::mutex> lock(this->stderr_access);
          std::cerr << "MergeBuffers::insert(): Thread " << omp_get_thread_num()
                    << ": Inserted the values to buffer " << i << std::endl;
        }
        return;
      }
      else
      {
        temp_buffer.swap(this->merge_buffers[i]);
      }
    }
    thread_buffer = buffer_type(thread_buffer, temp_buffer);
  }

  // All slots were full, write the merged buffer to a file.
  this->write(thread_buffer);
}

void
MergeBuffers::flush()
{
  for(size_type i = 1; i < this->merge_buffers.size(); i++)
  {
    this->merge_buffers[i] = buffer_type(this->merge_buffers[i], this->merge_buffers[i - 1]);
  }
  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    std::lock_guard<std::mutex> lock(this->stderr_access);
    std::cerr << "buildRA(): Flushing " << this->merge_buffers.back().size() << " values to disk" << std::endl;
  }
  this->write(this->merge_buffers.back());
}

void
MergeBuffers::write(buffer_type& buffer)
{
  if(buffer.empty()) { return; }

  std::string filename;
  size_type buffer_values = buffer.size(), buffer_bytes = buffer.bytes();
  {
    std::lock_guard<std::mutex> lock(this->ra_lock);
    this->ra.addFile(buffer);
  }
  buffer.write(filename); buffer.clear();

  double ra_done, ra_gb;
  {
    std::lock_guard<std::mutex> lock(this->ra_lock);
    this->ra_values += buffer_values;
    this->ra_bytes += buffer_bytes + sizeof(size_type);
    ra_done = (100.0 * this->ra_values) / this->final_size;
    ra_gb = inGigabytes(this->ra_bytes);
  }

  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    std::lock_guard<std::mutex> lock(this->stderr_access);
    std::cerr << "MergeBuffers::write(): Thread " << omp_get_thread_num()
              << ": Wrote " << buffer_values << " values to the rank array" << std::endl;
    std::cerr << "MergeBuffers::write(): " << ra_done << "% done; RA size " << ra_gb << " GB" << std::endl;
  }
}

//------------------------------------------------------------------------------

} // namespace gbwt
