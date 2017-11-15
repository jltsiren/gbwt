/*
  Copyright (c) 2017 Jouni Siren
  Copyright (c) 2017 Genome Research Ltd.

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

#ifndef GBWT_SUPPORT_H
#define GBWT_SUPPORT_H

#include <gbwt/utils.h>

namespace gbwt
{

/*
  support.h: Public support structures.
*/

//------------------------------------------------------------------------------

/*
  A simple encoding between (node id, orientation) <-> node_type. Not used by GBWT itself.
*/
struct Node
{
  const static node_type REVERSE_MASK = 0x1;
  const static size_type ID_SHIFT     = 1;

  static size_type id(node_type node) { return (node >> ID_SHIFT); }
  static bool is_reverse(node_type node) { return (node & REVERSE_MASK); }
  static node_type encode(size_type node_id, bool reversed) { return ((node_id << ID_SHIFT) | reversed); }
  static node_type reverse(node_type node) { return (node ^ REVERSE_MASK); }
};

/*
  A simple encoding between (path id, orientation) <-> size_type. Not used by GBWT itself.
*/
struct Path
{
  const static size_type REVERSE_MASK = 0x1;
  const static size_type ID_SHIFT     = 1;

  static size_type id(size_type path) { return (path >> ID_SHIFT); }
  static bool is_reverse(size_type path) { return (path & REVERSE_MASK); }
  static size_type encode(size_type path_id, bool reversed) { return ((path_id << ID_SHIFT) | reversed); }
  static size_type reverse(size_type path) { return (path ^ REVERSE_MASK); }
};

//------------------------------------------------------------------------------

/*
  The part of the BWT corresponding to a single node (the suffixes starting with / the
  prefixes ending with that node).

  - Incoming edges are sorted by the source node.
  - Outgoing edges are sorted by the destination node.
  - Sampled sequence ids are sorted by the offset.
*/

struct DynamicRecord
{
  typedef gbwt::size_type size_type;

  size_type                body_size;
  std::vector<edge_type>   incoming, outgoing;
  std::vector<run_type>    body;
  std::vector<sample_type> ids;

//------------------------------------------------------------------------------

  DynamicRecord() : body_size(0) {}

  size_type size() const { return this->body_size; }
  bool empty() const { return (this->size() == 0); }
  size_type indegree() const { return this->incoming.size(); }
  size_type outdegree() const { return this->outgoing.size(); }
  size_type runs() const { return this->body.size(); }
  size_type samples() const { return this->ids.size(); }

  void clear();
  void swap(DynamicRecord& another);

//------------------------------------------------------------------------------

  // Sort the outgoing edges if they are not sorted.
  void recode();

  // Write the compressed representation.
  void writeBWT(std::vector<byte_type>& data) const;

//------------------------------------------------------------------------------

  // Returns (node, LF(i, node)) or invalid_edge() if the offset is invalid.
  edge_type LF(size_type i) const;

  // As above, but also sets 'run_end' to the last offset of the current run.
  edge_type runLF(size_type i, size_type& run_end) const;

  // Returns invalid_offset() if there is no edge to the destination.
  size_type LF(size_type i, node_type to) const;

  // Returns Range::empty_range() if the range is empty or the destination is invalid.
  range_type LF(range_type range, node_type to) const;

  // Returns BWT[i] within the record.
  node_type operator[](size_type i) const;

//------------------------------------------------------------------------------

  // Maps successor nodes to outranks.
  rank_type edgeTo(node_type to) const;

  // These assume that 'outrank' is a valid outgoing edge.
  node_type successor(rank_type outrank) const { return this->outgoing[outrank].first; }
#ifdef GBWT_SAVE_MEMORY
  short_type& offset(rank_type outrank) { return this->outgoing[outrank].second; }
#else
  size_type& offset(rank_type outrank) { return this->outgoing[outrank].second; }
#endif
  size_type offset(rank_type outrank) const { return this->outgoing[outrank].second; }

//------------------------------------------------------------------------------

  // Return the first node >= 'from' with an incoming edge to this node.
  rank_type findFirst(node_type from) const;

  // These assume that 'inrank' is a valid incoming edge.
  node_type predecessor(rank_type inrank) const { return this->incoming[inrank].first; }
#ifdef GBWT_SAVE_MEMORY
  short_type& count(rank_type inrank) { return this->incoming[inrank].second; }
#else
  size_type& count(rank_type inrank) { return this->incoming[inrank].second; }
#endif
  size_type count(rank_type inrank) const { return this->incoming[inrank].second; }

  // Increment the count of the incoming edge from 'from'.
  void increment(node_type from);

  // Add a new incoming edge.
  void addIncoming(edge_type inedge);

//------------------------------------------------------------------------------

  // Returns the first sample at offset >= i or ids.end() if there is no sample.
  std::vector<sample_type>::const_iterator nextSample(size_type i) const;

};  // struct DynamicRecord

std::ostream& operator<<(std::ostream& out, const DynamicRecord& record);

//------------------------------------------------------------------------------

struct CompressedRecord
{
  typedef gbwt::size_type size_type;

  std::vector<edge_type> outgoing;
  const byte_type*       body;
  size_type              data_size;

  CompressedRecord();
  CompressedRecord(const std::vector<byte_type>& source, size_type start, size_type limit);

  size_type size() const; // Expensive.
  bool empty() const { return (this->size() == 0); }
  size_type runs() const; // Expensive.
  size_type outdegree() const { return this->outgoing.size(); }

  // Returns (node, LF(i, node)) or invalid_edge() if the offset is invalid.
  edge_type LF(size_type i) const;

  // As above, but also sets 'run_end' to the last offset of the current run.
  edge_type runLF(size_type i, size_type& run_end) const;

  // Returns invalid_offset() if there is no edge to the destination.
  size_type LF(size_type i, node_type to) const;

  // Returns Range::empty_range() if the range is empty or the destination is invalid.
  range_type LF(range_type range, node_type to) const;

  // Returns BWT[i] within the record.
  node_type operator[](size_type i) const;

  // Maps successor nodes to outranks.
  rank_type edgeTo(node_type to) const;

  // These assume that 'outrank' is a valid outgoing edge.
  node_type successor(rank_type outrank) const { return this->outgoing[outrank].first; }
  size_type offset(rank_type outrank) const { return this->outgoing[outrank].second; }
};

//------------------------------------------------------------------------------

struct RecordArray
{
  typedef gbwt::size_type size_type;

  size_type                        records;
  sdsl::sd_vector<>                index;
  sdsl::sd_vector<>::select_1_type select;
  std::vector<byte_type>           data;

  RecordArray();
  RecordArray(const RecordArray& source);
  RecordArray(RecordArray&& source);
  ~RecordArray();

  explicit RecordArray(const std::vector<DynamicRecord>& bwt);

  // Set the number of records, build the data manually, and give the offsets to build the index.
  explicit RecordArray(size_type array_size);
  void buildIndex(const std::vector<size_type>& offsets);

  void swap(RecordArray& another);
  RecordArray& operator=(const RecordArray& source);
  RecordArray& operator=(RecordArray&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  // 0-based indexing.
  size_type start(size_type record) const { return this->select(record + 1); }
  size_type limit(size_type record) const
  {
    return (record + 1 < this->records ? this->select(record + 2) : this->data.size());
  }

private:
  void copy(const RecordArray& source);
};

//------------------------------------------------------------------------------

struct DASamples
{
  typedef gbwt::size_type size_type;

  // Does node i have samples?
  sdsl::bit_vector                 sampled_records;
  sdsl::bit_vector::rank_1_type    record_rank;

  // Map from record ranks to BWT offsets.
  sdsl::sd_vector<>                bwt_ranges;
  sdsl::sd_vector<>::select_1_type bwt_select;

  // Sampled offsets.
  sdsl::sd_vector<>                sampled_offsets;
  sdsl::sd_vector<>::rank_1_type   sample_rank;

  sdsl::int_vector<0>              array;

  DASamples();
  DASamples(const DASamples& source);
  DASamples(DASamples&& source);
  ~DASamples();

  explicit DASamples(const std::vector<DynamicRecord>& bwt);

  void swap(DASamples& another);
  DASamples& operator=(const DASamples& source);
  DASamples& operator=(DASamples&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  size_type size() const { return this->array.size(); }

  // Returns invalid_sequence() if there is no sample.
  size_type tryLocate(size_type record, size_type offset) const;

  // Returns the first sample at >= offset or invalid_sample() if there is no sample.
  sample_type nextSample(size_type record, size_type offset) const;

private:
  void copy(const DASamples& source);
  void setVectors();
};

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_SUPPORT_H
