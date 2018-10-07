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

#include <gbwt/internal.h>

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

  ByteArray data;
  size_type value_count;

private:
  void copy(const GapArray& source)
  {
    this->data = source.data;
    this->value_count = source.value_count;
  }

  void writeValue(edge_type curr, edge_type& prev)
  {
    ByteCode::write(this->data, curr.first - prev.first);
    if(curr.first != prev.first) { prev.second = 0; } // Node changed, set previous offset to 0.
    ByteCode::write(this->data, curr.second - prev.second);
    prev = curr;
  }
};  // class GapArray

template<> GapArray<BlockArray>::GapArray(std::vector<edge_type>& source);
template<> GapArray<BlockArray>::GapArray(GapArray& a, GapArray& b);
template<> void GapArray<BlockArray>::clear();
template<> void GapArray<BlockArray>::write(const std::string& filename);

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

  GapIterator(const GapIterator& source) :
    array(source.array), pos(source.pos), data_pointer(source.data_pointer), value(source.value)
  {
  }

  edge_type operator*() const { return this->value; }
  void operator++() { this->pos++; this->read(); }
  bool end() const { return (this->pos >= this->array->size()); }

  GapArray<ByteArray>* array;
  size_type pos, data_pointer;
  edge_type value;

private:
  // There is a specialized destructive read() BlockArray.
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
};  // class GapIterator

template<> void GapIterator<BlockArray>::read();

//------------------------------------------------------------------------------

/*
  A RankArray is a number of temporary files containing GapArrays on disk. When the object
  is deleted, the files are also deleted.
*/

class RankArray
{
public:
  typedef GapArray<sdsl::int_vector_buffer<8>> array_type;
  typedef array_type::size_type                size_type;
  typedef array_type::iterator                 iterator;

  const static std::string TEMP_FILE_PREFIX;  // "ranks"

  RankArray();
  ~RankArray();

  // Also clears the array. Do not call after opening the RankArray.
  void addFile(GapArray<BlockArray>& array);

  void open();
  void close();

  // Iterator operations.
  edge_type operator*() const { return *(this->iterators[0]); }
  void operator++() { ++(this->iterators[0]); this->down(0); }
  bool end() const { return (this->empty() || this->iterators[0].end()); }

  size_type size() const { return this->filenames.size(); }
  size_type empty() const { return (this->size() == 0); }

  std::vector<std::string> filenames;
  std::vector<size_type>   value_counts;

  std::vector<array_type> inputs;
  std::vector<iterator>   iterators;

private:
  /*
    Heap operations.
  */
  static size_type parent(size_type i) { return (i - 1) / 2; }
  static size_type left(size_type i) { return 2 * i + 1; }
  static size_type right(size_type i) { return 2 * i + 2; }

  size_type smaller(size_type i, size_type j) const
  {
    return (this->iterators[j].value < this->iterators[i].value ? j : i);
  }

  void down(size_type i)
  {
    while(left(i) < this->size())
    {
      size_type next = this->smaller(i, left(i));
      if(right(i) < this->size()) { next = this->smaller(next, right(i)); }
      if(next == i) { return; }
      std::swap(this->iterators[i], this->iterators[next]);
      i = next;
    }
  }

  void heapify();

  RankArray(const RankArray&) = delete;
  RankArray& operator= (const RankArray&) = delete;
};

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_BWTMERGE_H
