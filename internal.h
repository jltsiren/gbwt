/*
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

#include "support.h"

namespace gbwt
{

/*
  support.h: Internal support structures.
*/

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

  const static size_type DATA_BITS  = 7;
  const static code_type DATA_MASK  = 0x7F;
  const static code_type NEXT_BYTE  = 0x80;

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
  inline void write(ByteArray& array, run_type run) { this->write(array, run.first, run.second); }

  inline code_type encodeBasic(value_type value, size_type length)
  {
    return value + this->sigma * (length - 1);
  }

  inline run_type decodeBasic(code_type code)
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

  inline size_type size() const { return this->total_size; }

  inline void insert(run_type run)
  {
    this->total_size += run.second; counts[run.first] += run.second;
    if(run.first == accumulator.first) { accumulator.second += run.second; }
    else { this->flush(); this->accumulator = run; }
  }

  inline void insert(rank_type outrank) { this->insert(run_type(outrank, 1)); }

  void flush();

  inline void addEdge() { this->counts.push_back(0); }
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

  // Create a sequence where endmarker[source_pos] == node.
  Sequence(node_type node, size_type seq_id, size_type source_pos);

  // Sort by reverse prefixes text[..pos+1].
  inline bool operator<(const Sequence& another) const
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

  FIXME a single iterator with a RankCalculator as a template parameter.
*/

struct CompressedRecordIterator
{
  explicit CompressedRecordIterator(const CompressedRecord& source) :
    record(source), decoder(source.outdegree()),
    record_offset(0), curr_offset(0), next_offset(0)
  {
    this->read();
  }

  inline bool end() const { return (this->curr_offset >= this->record.data_size); }
  inline void operator++() { this->curr_offset = this->next_offset; this->read(); }

  inline run_type operator*() const { return this->run; }
  inline const run_type* operator->() { return &(this->run); }

  // After the current run.
  inline size_type offset() const { return this->record_offset; }

  const CompressedRecord& record;
  Run                     decoder;

  size_type               record_offset;
  size_type               curr_offset, next_offset;
  run_type                run;

private:
  inline void read()
  {
    if(!(this->end()))
    {
      this->run = this->decoder.read(this->record.body, this->next_offset);
      this->record_offset += this->run.second;
    }
  }
};

struct CompressedRecordRankIterator
{
  explicit CompressedRecordRankIterator(const CompressedRecord& source, rank_type outrank) :
    record(source), decoder(source.outdegree()),
    record_offset(0), curr_offset(0), next_offset(0),
    value(outrank), result(source.offset(outrank))
  {
    this->read();
  }

  inline bool end() const { return (this->curr_offset >= this->record.data_size); }
  inline void operator++() { this->curr_offset = this->next_offset; this->read(); }

  inline run_type operator*() const { return this->run; }
  inline const run_type* operator->() { return &(this->run); }

  // After the current run.
  inline size_type offset() const { return this->record_offset; }
  inline size_type rank() const { return this->result; }

  // Intended for positions i covered by the current run.
  inline size_type rankAt(size_type i) const
  {
    size_type temp = this->rank();
    if(i < this->offset() && this->run.first == value) { temp -= (this->offset() - i); }
    return temp;
  }

  const CompressedRecord& record;
  Run                     decoder;

  size_type               record_offset;
  size_type               curr_offset, next_offset;
  run_type                run;

  rank_type               value;
  size_type               result;

private:
  inline void read()
  {
    if(!(this->end()))
    {
      this->run = this->decoder.read(this->record.body, this->next_offset);
      this->record_offset += this->run.second;
      if(this->run.first == value) { this->result += this->run.second; }
    }
  }
};

struct CompressedRecordFullIterator
{
  explicit CompressedRecordFullIterator(const CompressedRecord& source) :
    record(source), decoder(source.outdegree()), ranks(source.outgoing),
    record_offset(0), curr_offset(0), next_offset(0)
  {
    this->read();
  }

  inline bool end() const { return (this->curr_offset >= this->record.data_size); }
  inline void operator++() { this->curr_offset = this->next_offset; this->read(); }

  inline run_type operator*() const { return this->run; }
  inline const run_type* operator->() { return &(this->run); }

  // After the current run.
  inline size_type offset() const { return this->record_offset; }
  inline size_type rank() const { return this->rank(this->run.first); }
  inline size_type rank(rank_type outrank) const { return this->ranks[outrank].second; }
  inline edge_type edge() const { return this->edge(this->run.first); }
  inline edge_type edge(rank_type outrank) const { return this->ranks[outrank]; }

  const CompressedRecord& record;
  Run                     decoder;
  std::vector<edge_type>  ranks;

  size_type               record_offset;
  size_type               curr_offset, next_offset;
  run_type                run;

private:
  inline void read()
  {
    if(!(this->end()))
    {
      this->run = this->decoder.read(this->record.body, this->next_offset);
      this->record_offset += this->run.second;
      this->ranks[this->run.first].second += this->run.second;
    }
  }
};

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_INTERNAL_H
