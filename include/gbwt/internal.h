/*
  Copyright (c) 2017, 2018 Jouni Siren
  Copyright (c) 2015, 2016, 2017 Genome Research Ltd.

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

#ifndef GBWT_INTERNAL_H
#define GBWT_INTERNAL_H

#include <array>

#include <gbwt/support.h>

namespace gbwt
{

/*
  support.h: Internal support structures.
*/

//------------------------------------------------------------------------------

/*
  Reads and writes are done in blocks of 1048576 elements to avoid the bug with
  large reads in GCC / macOS.
*/

struct DiskIO
{
  const static size_type block_size;

  template<class Element>
  static bool read(std::istream& in, Element* data, size_type n = 1)
  {
    for(size_type offset = 0; offset < n; offset += block_size)
    {
      size_type bytes = std::min(block_size, n - offset) * sizeof(Element);
      in.read(reinterpret_cast<char*>(data + offset), bytes);
      size_type read_bytes = in.gcount();
      if(read_bytes < bytes) { return false; }
    }
    return true;
  }

  template<class Element>
  static void write(std::ostream& out, const Element* data, size_type n = 1)
  {
    for(size_type offset = 0; offset < n; offset += block_size)
    {
      size_type bytes = std::min(block_size, n - offset) * sizeof(Element);
      out.write(reinterpret_cast<const char*>(data + offset), bytes);
      if(out.fail())
      {
        std::cerr << "DiskIO::write(): Write failed" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  }
};

// Serialize an std::vector of integers.
template<class Element>
size_type
serializeVector(const std::vector<Element>& data, std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "")
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(data));
  size_type written_bytes = 0;

  size_type data_size = data.size();
  written_bytes += sdsl::write_member(data_size, out, child, "size");

  if(data_size > 0)
  {
    sdsl::structure_tree_node* data_node =
      sdsl::structure_tree::add_child(child, "data", sdsl::util::class_name(data[0]));
    DiskIO::write(out, data.data(), data_size);
    sdsl::structure_tree::add_size(data_node, data_size * sizeof(Element));
    written_bytes += data_size * sizeof(Element);
  }

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

// Specialization for a vector of strings.
template<>
size_type serializeVector<std::string>(const std::vector<std::string>& data, std::ostream& out, sdsl::structure_tree_node* v, std::string name);

// Load an std::vector of integers.
template<class Element>
void
loadVector(std::vector<Element>& data, std::istream& in)
{
  size_type data_size = 0;
  sdsl::read_member(data_size, in);

  data.resize(data_size);
  if(data_size > 0)
  {
    DiskIO::read(in, data.data(), data_size);
  }
}

// Specialization for a vector of strings.
template<>
void loadVector<std::string>(std::vector<std::string>& data, std::istream& in);

//------------------------------------------------------------------------------

/*
  Encodes unsigned integers as byte sequences. Each byte contains 7 bits of data
  and one bit telling whether the encoding continues in the next byte. The data is
  stored in LSB order.
*/

struct ByteCode
{
  typedef gbwt::size_type value_type;
  typedef gbwt::byte_type code_type;

  constexpr static size_type DATA_BITS  = 7;
  constexpr static code_type DATA_MASK  = 0x7F;
  constexpr static code_type NEXT_BYTE  = 0x80;

  /*
    Reads the next value and updates i to point to the byte after the value.
  */
  template<class ByteArray>
  static value_type read(ByteArray& array, size_type& i)
  {
    size_type offset = 0;
    value_type res = array[i] & DATA_MASK;
    while(array[i] & NEXT_BYTE)
    {
      i++; offset += DATA_BITS;
      res += ((value_type)(array[i] & DATA_MASK)) << offset;
    }
    i++;
    return res;
  }

  /*
    Encodes the value and stores it in the array using push_back().
  */
  template<class ByteArray>
  static void write(ByteArray& array, value_type value)
  {
    while(value > DATA_MASK)
    {
      array.push_back((value & DATA_MASK) | NEXT_BYTE);
      value >>= DATA_BITS;
    }
    array.push_back(value);
  }
};

//------------------------------------------------------------------------------

/*
  Run-length encoding using ByteCode. Run lengths and alphabet size are assumed to be > 0.
*/

struct Run
{
  typedef ByteCode::value_type value_type;
  typedef ByteCode::code_type  code_type;

  size_type sigma, run_continues;

  explicit Run(size_type alphabet_size);

  /*
    Returns (value, run length) and updates i to point past the run.
  */
  template<class ByteArray>
  run_type read(ByteArray& array, size_type& i)
  {
    run_type run;
    if(this->run_continues == 0)
    {
      run.first = ByteCode::read(array, i);
      run.second = ByteCode::read(array, i) + 1;
    }
    else
    {
      run = this->decodeBasic(array[i]); i++;
      if(run.second >= this->run_continues) { run.second += ByteCode::read(array, i); }
    }
    return run;
  }

  /*
    Encodes the run and stores it in the array using push_back().
  */
  template<class ByteArray>
  void write(ByteArray& array, value_type value, size_type length)
  {
    if(this->run_continues == 0)
    {
      ByteCode::write(array, value);
      ByteCode::write(array, length - 1);
    }
    else if(length < this->run_continues)
    {
      array.push_back(this->encodeBasic(value, length));
    }
    else
    {
      array.push_back(this->encodeBasic(value, this->run_continues));
      ByteCode::write(array, length - this->run_continues);
    }
  }

  template<class ByteArray>
  void write(ByteArray& array, run_type run) { this->write(array, run.first, run.second); }

  code_type encodeBasic(value_type value, size_type length)
  {
    return value + this->sigma * (length - 1);
  }

  run_type decodeBasic(code_type code)
  {
    return run_type(code % this->sigma, code / this->sigma + 1);
  }
};

//------------------------------------------------------------------------------

/*
  A support structure for run-length encoding outrank sequences.
*/

struct RunMerger
{
  size_type              total_size;
  run_type               accumulator;
  std::vector<run_type>  runs;
  std::vector<size_type> counts;

  RunMerger(size_type sigma) : total_size(0), accumulator(0, 0), counts(sigma) {}

  size_type size() const { return this->total_size; }

  void insert(run_type run)
  {
    this->total_size += run.second; counts[run.first] += run.second;
    if(run.first == accumulator.first) { accumulator.second += run.second; }
    else { this->flush(); this->accumulator = run; }
  }

  void insert(rank_type outrank) { this->insert(run_type(outrank, 1)); }

  void flush();

  void addEdge() { this->counts.push_back(0); }
};

//------------------------------------------------------------------------------

/*
  A text iterator corresponding to a sequence. Used for GBWT construction.
*/

struct Sequence
{
  size_type id;
  node_type curr, next;
  size_type offset; // Offset in the current record.
  size_type pos;    // Position in the text or offset in the source record.

  Sequence();

  // Create a sequence starting from text[i].
  Sequence(const text_type& text, size_type i, size_type seq_id);
  Sequence(const vector_type& text, size_type i, size_type seq_id);

  // Create a sequence where endmarker[source_pos] == node.
  Sequence(node_type node, size_type seq_id, size_type source_pos);

  // Sort by reverse prefixes text[..pos+1].
  bool operator<(const Sequence& another) const
  {
    if(this->next != another.next) { return (this->next < another.next); }
    if(this->curr != another.curr) { return (this->curr < another.curr); }
    return (this->offset < another.offset);
  }
};

//------------------------------------------------------------------------------

/*
  Iterators for CompressedRecords.

  - CompressedRecordIterator is the fastest, as it only iterates over the runs.
  - CompressedRecordRankIterator keeps track of the rank for one successor node.
  - CompressedRecordFullIterator is the slowest, as it keeps track of the ranks
    for all successor nodes.
  - CompressedRecordArrayIterator is a faster version of the full iterator for
    records with outdegree <= MAX_OUTDEGREE_FOR_ARRAY.
*/

template<class RankSupport>
struct CompressedRecordGenericIterator
{
  explicit CompressedRecordGenericIterator(const CompressedRecord& source, rank_type outrank = 0) :
    record(source), decoder(source.outdegree()),
    record_offset(0), curr_offset(0), next_offset(0),
    rank_support(source, outrank)
  {
    this->read();
  }

  bool end() const { return (this->curr_offset >= this->record.data_size); }
  void operator++() { this->curr_offset = this->next_offset; this->read(); }

  // Read while offset < i.
  void readUntil(size_type i)
  {
    while(this->offset() < i && !(this->end()))
    {
      this->curr_offset = this->next_offset;
      this->readUnsafe();
    }
  }

  // Read while offset <= i.
  void readPast(size_type i)
  {
    while(this->offset() <= i && !(this->end()))
    {
      this->curr_offset = this->next_offset;
      this->readUnsafe();
    }
  }

  run_type operator*() const { return this->run; }
  const run_type* operator->() const { return &(this->run); }

  // After the current run.
  size_type offset() const { return this->record_offset; }
  size_type rank() const { return this->rank_support.rank(*this); }
  size_type rank(rank_type outrank) const { return this->rank_support.rank(outrank); }
  edge_type edge() const { return this->rank_support.edge(*this); }
  edge_type edge(rank_type outrank) const { return this->rank_support.edge(outrank); }

  // Intended for positions i covered by or after the current run. May advance the iterator.
  size_type rankAt(size_type i) { return this->rank_support.rankAt(*this, i); }
  edge_type edgeAt(size_type i) { return this->rank_support.edgeAt(*this, i); }

  const CompressedRecord& record;
  Run                     decoder;

  size_type               record_offset;
  size_type               curr_offset, next_offset;
  run_type                run;

  RankSupport             rank_support;

private:
  void read()
  {
    if(!(this->end()))
    {
      this->readUnsafe();
    }
  }

  void readUnsafe()
  {
    this->run = this->decoder.read(this->record.body, this->next_offset);
    this->record_offset += this->run.second;
    this->rank_support.handle(this->run);
  }
};

struct DummyRankSupport
{
  typedef CompressedRecordGenericIterator<DummyRankSupport> Iterator;

  DummyRankSupport(const CompressedRecord&, rank_type) {}

  void handle(run_type) {}

  size_type rank(const Iterator&) const { return invalid_offset(); }
  size_type rank(rank_type) const { return invalid_offset(); }
  edge_type edge(const Iterator&) const { return invalid_edge(); }
  edge_type edge(rank_type) const { return invalid_edge(); }

  size_type rankAt(Iterator&, size_type) { return invalid_offset(); }
  edge_type edgeAt(Iterator&, size_type) { return invalid_edge(); }
};

struct OccurrenceCounter
{
  typedef CompressedRecordGenericIterator<OccurrenceCounter> Iterator;

  OccurrenceCounter(const CompressedRecord& source, rank_type outrank) : value(outrank), result(source.offset(outrank)) {}

  void handle(run_type run)
  {
    if(run.first == this->value) { this->result += run.second; }
  }

  size_type rank(const Iterator&) const { return this->result; }
  size_type rank(rank_type) const { return invalid_offset(); }
  edge_type edge(const Iterator&) const { return invalid_edge(); }
  edge_type edge(rank_type) const { return invalid_edge(); }

  size_type rankAt(Iterator& iter, size_type i)
  {
    iter.readUntil(i);
    size_type temp = this->rank(iter);
    if(i < iter.offset() && iter->first == this->value) { temp -= (iter.offset() - i); }
    return temp;
  }

  edge_type edgeAt(Iterator&, size_type) { return invalid_edge(); }

  rank_type               value;
  size_type               result;
};

inline void
copyEdges(const std::vector<edge_type>& from, std::vector<edge_type>& to)
{
  to = from;
}

constexpr size_type MAX_OUTDEGREE_FOR_ARRAY = 4;

inline void
copyEdges(const std::vector<edge_type>& from, std::array<edge_type, MAX_OUTDEGREE_FOR_ARRAY>& to)
{
  for(size_type i = 0; i < from.size(); i++) { to[i] = from[i]; }
}

template<class ArrayType>
struct FullRankSupport
{
  typedef CompressedRecordGenericIterator<FullRankSupport<ArrayType>> Iterator;

  FullRankSupport(const CompressedRecord& source, rank_type)
  {
    copyEdges(source.outgoing, this->ranks);
  }

  void handle(run_type run)
  {
    this->ranks[run.first].second += run.second;
  }

  size_type rank(const Iterator& iter) const { return this->rank(iter->first); }
  size_type rank(rank_type outrank) const { return this->ranks[outrank].second; }
  edge_type edge(const Iterator& iter) const { return this->edge(iter->first); }
  edge_type edge(rank_type outrank) const { return this->ranks[outrank]; }

  size_type rankAt(Iterator& iter, size_type i)
  {
    iter.readPast(i); // We read past offset i to get BWT[i].
    if(iter.end()) { return invalid_offset(); }
    return this->rank(iter) - (iter.offset() - i);
  }

  edge_type edgeAt(Iterator& iter, size_type i)
  {
    iter.readPast(i); // We read past offset i to get BWT[i].
    if(iter.end()) { return invalid_edge(); }
    edge_type temp = this->edge(iter);
    temp.second -= (iter.offset() - i);
    return temp;
  }

  ArrayType ranks;
};

typedef CompressedRecordGenericIterator<DummyRankSupport> CompressedRecordIterator;
typedef CompressedRecordGenericIterator<OccurrenceCounter> CompressedRecordRankIterator;
typedef CompressedRecordGenericIterator<FullRankSupport<std::vector<edge_type>>> CompressedRecordFullIterator;
typedef CompressedRecordGenericIterator<FullRankSupport<std::array<edge_type, MAX_OUTDEGREE_FOR_ARRAY>>> CompressedRecordArrayIterator;

//------------------------------------------------------------------------------

/*
  Iterator for DASamples. The iterator does not care about records. If the record
  for the current sample starts at offset i, the correct sample_type is
  (iter.offset() - i, *iter).
*/

struct SampleIterator
{
  explicit SampleIterator(const DASamples& source) :
    data(source),
    pos(0), sample_offset(0),
    offset_select(&(source.sampled_offsets))
  {
    this->update();
  }

  bool end() const { return (this->pos >= this->data.size()); }
  void operator++() { this->pos++; this->update(); }

  size_type operator*() const { return this->data.array[this->pos]; }

  size_type offset() const { return this->sample_offset; }

  const DASamples& data;
  size_type        pos, sample_offset;

private:
  sdsl::sd_vector<>::select_1_type offset_select;

  void update()
  {
    if(!(this->end()))
    {
      this->sample_offset = this->offset_select(this->pos + 1);
    }
  }
};

/*
  Iterator for sampled ranges in DASamples.
*/

struct SampleRangeIterator
{
  explicit SampleRangeIterator(const DASamples& source) :
    data(source),
    record_id(0), record_rank(0),
    record_start(0), record_limit(0)
  {
    this->advance();
  }

  bool end() const { return (this->record_id >= this->data.records()); }
  void operator++() { this->record_id++; this->record_rank++; this->advance(); }

  size_type record() const { return this->record_id; }
  size_type rank() const { return this->record_rank; }
  size_type start() const { return this->record_start; }
  size_type limit() const { return this->record_limit; }
  size_type length() const { return this->limit() - this->start(); }

  const DASamples& data;
  size_type        record_id, record_rank;
  size_type        record_start, record_limit;

private:
  void advance()
  {
    while(!(this->end()))
    {
      if(this->data.isSampled(this->record_id))
      {
        this->record_start = this->record_limit;
        this->record_limit = this->data.limit(this->record_rank);
        return;
      }
      record_id++;
    }
  }
};

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_INTERNAL_H
