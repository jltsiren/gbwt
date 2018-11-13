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
  // Use a block size of 4 bytes per element in the producer buffer.
  // This way we should avoid starving the producer thread.
  array.data = sdsl::int_vector_buffer<8>(filename, std::ios::in, 4 * PRODUCER_BUFFER_SIZE);
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

RankArray::RankArray() :
  tournament_tree(1, tree_type(invalid_edge(), 0))
{
}

RankArray::~RankArray()
{
  this->close();
  for(size_type i = 0; i < this->size(); i++) { TempFile::remove(this->filenames[i]); }
}

std::string
RankArray::addFile(GapArray<BlockArray>& array)
{
  std::string filename = TempFile::getName(TEMP_FILE_PREFIX);
  this->filenames.push_back(filename);
  this->value_counts.push_back(array.size());
  return this->filenames.back();
}

void
RankArray::open()
{
  this->close();
  this->inputs = std::vector<array_type>(this->size());
  this->iterators = std::vector<iterator>(this->size());
  this->buffers = std::vector<ProducerBuffer<iterator>*>(this->size());

  for(size_type i = 0; i < this->size(); i++)
  {
    gbwt::open(this->inputs[i], this->filenames[i], this->value_counts[i]);
    this->iterators[i] = iterator(this->inputs[i]);
    this->buffers[i] = new ProducerBuffer<iterator>(this->iterators[i]);
  }

  this->initTree();
}

void
RankArray::close()
{
  for(size_type i = 0; i < this->buffers.size(); i++)
  {
    delete this->buffers[i]; this->buffers[i] = nullptr;
  }
  this->buffers.clear();
  this->iterators.clear();
  this->inputs.clear();
  this->tournament_tree = std::vector<tree_type>(1, tree_type(invalid_edge(), 0));
}

void RankArray::operator++()
{
  // Advance the active file.
  size_type pos = this->tournament_tree.back().second;
  this->buffers[pos]->operator++();
  this->tournament_tree[pos].first = this->buffers[pos]->operator*();

  // Update the tournament tree.
  size_type level_size = this->leaves, level_offset = 0;
  while(level_size > 1)
  {
    size_type next_offset = level_offset + level_size;
    this->tournament_tree[next_offset + pos / 2] = this->smaller(level_offset + pos);
    level_offset = next_offset;
    pos /= 2; level_size /= 2;
  }
}

void
RankArray::initTree()
{
  // The number of leaves must be a power of two.
  this->leaves = 1;
  while(this->leaves < this->size()) { this->leaves *= 2; }

  // Build the tournament tree.
  this->tournament_tree = std::vector<tree_type>(2 * this->leaves - 1, tree_type(invalid_edge(), 0));
  for(size_type i = 0; i < this->size(); i++) { this->tournament_tree[i].first = this->buffers[i]->operator*(); }
  for(size_type i = 0; i < this->leaves; i++) { this->tournament_tree[i].second = i; }
  size_type level_size = this->leaves, level_offset = 0;
  while(level_size > 1)
  {
    size_type next_offset = level_offset + level_size;
    for(size_type i = 0; i < level_size; i += 2)
    {
      this->tournament_tree[next_offset + i / 2] = this->smaller(level_offset + i);
    }
    level_offset = next_offset;
    level_size /= 2;
  }
}

//------------------------------------------------------------------------------

MergeBuffers::MergeBuffers(size_type expected_size, size_type num_threads, const MergeParameters& params) :
  parameters(params),
  pos_buffers(num_threads), thread_buffers(num_threads),
  merge_buffers(params.merge_buffers),
  ra_values(0), ra_bytes(0), final_size(expected_size)
{
}

MergeBuffers::~MergeBuffers()
{
}

void
MergeBuffers::flush()
{
  // First insert all thread-specific buffers.
  #pragma omp parallel for schedule(static, 1)
  for(size_type thread = 0; thread < this->threads(); thread++)
  {
    this->insert(this->pos_buffers[thread], this->thread_buffers[thread], true);
  }

  // Then merge and write the global buffers.
  for(size_type i = 1; i < this->merge_buffers.size(); i++)
  {
    this->merge_buffers[i] = buffer_type(this->merge_buffers[i], this->merge_buffers[i - 1]);
  }
  size_type buffer_values = this->merge_buffers.back().size();
  this->write(this->merge_buffers.back());

  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    std::lock_guard<std::mutex> lock(this->stderr_access);
    std::cerr << "MergeBuffers::flush(): Wrote " << buffer_values << " values to disk" << std::endl;
  }
}

void
MergeBuffers::insert(std::vector<edge_type>& pos_buffer, buffer_type& thread_buffer, bool force_merge)
{
  // Compress the pos_buffer and merge it into thread_buffer.
  if(!(pos_buffer.empty()))
  {
    buffer_type temp_buffer(pos_buffer); pos_buffer.clear();
    thread_buffer = buffer_type(thread_buffer, temp_buffer);
  }

  // Should we merge thread_buffer into the merge buffers?
  if(thread_buffer.empty()) { return; }
  if(!force_merge && thread_buffer.bytes() < this->parameters.threadBufferBytes()) { return; }
  if(Verbosity::level >= Verbosity::FULL)
  {
    std::lock_guard<std::mutex> lock(this->stderr_access);
    std::cerr << "MergeBuffers::insert(): Thread " << omp_get_thread_num() << ": Inserting "
              << thread_buffer.size() << " values to the merge buffers" << std::endl;
  }

  // Merge with the existing merge buffers until we find an empty slot.
  for(size_type i = 0; i < this->merge_buffers.size(); i++)
  {
    buffer_type temp_buffer;
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
MergeBuffers::write(buffer_type& buffer)
{
  if(buffer.empty()) { return; }

  size_type buffer_values = buffer.size(), buffer_bytes = buffer.bytes();
  std::string filename;
  {
    std::lock_guard<std::mutex> lock(this->ra_lock);
    filename = this->ra.addFile(buffer);
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
