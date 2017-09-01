/*
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

#include "utils.h"

namespace gbwt
{

/*
  support.h: Public support structures.
*/

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

  inline size_type size() const { return this->body_size; }
  inline bool empty() const { return (this->size() == 0); }
  inline size_type indegree() const { return this->incoming.size(); }
  inline size_type outdegree() const { return this->outgoing.size(); }
  inline size_type runs() const { return this->body.size(); }
  inline size_type samples() const { return this->ids.size(); }

  void clear();
  void swap(DynamicRecord& another);

//------------------------------------------------------------------------------

  // Sort the outgoing edges if they are not sorted.
  void recode();

//------------------------------------------------------------------------------

  // Returns invalid_offset() if there is no edge to the destination.
  size_type LF(size_type i, node_type to) const;

  // Returns (node, LF(i, node)) or invalid_edge() if the offset is invalid.
  edge_type LF(size_type i) const;

  // Returns Range::empty_range() if the range is empty or the destination is invalid.
  range_type LF(range_type range, node_type to) const;

  // Returns BWT[i] within the record.
  node_type operator[](size_type i) const;

//------------------------------------------------------------------------------

  // Maps successor nodes to outranks.
  rank_type edgeTo(node_type to) const;

  // These assume that 'outrank' is a valid outgoing edge.
  inline node_type successor(rank_type outrank) const { return this->outgoing[outrank].first; }
#ifdef GBWT_SAVE_MEMORY
  inline short_type& offset(rank_type outrank) { return this->outgoing[outrank].second; }
#else
  inline size_type& offset(rank_type outrank) { return this->outgoing[outrank].second; }
#endif
  inline size_type offset(rank_type outrank) const { return this->outgoing[outrank].second; }

//------------------------------------------------------------------------------

  // Return the first node >= 'from' with an incoming edge to this node.
  rank_type findFirst(node_type from) const;

  // These assume that 'inrank' is a valid incoming edge.
  inline node_type predecessor(rank_type inrank) const { return this->incoming[inrank].first; }
#ifdef GBWT_SAVE_MEMORY
  inline short_type& count(rank_type inrank) { return this->incoming[inrank].second; }
#else
  inline size_type& count(rank_type inrank) { return this->incoming[inrank].second; }
#endif
  inline size_type count(rank_type inrank) const { return this->incoming[inrank].second; }

  // Increment the count of the incoming edge from 'from'.
  void increment(node_type from);

  // Add a new incoming edge.
  void addIncoming(edge_type inedge);

//------------------------------------------------------------------------------

};  // struct DynamicRecord

std::ostream& operator<<(std::ostream& out, const DynamicRecord& record);

//------------------------------------------------------------------------------

struct CompressedRecord
{
  typedef gbwt::size_type size_type;

  std::vector<edge_type> outgoing;
  const byte_type*       body;
  size_type              data_size;

  CompressedRecord(const std::vector<byte_type>& source, size_type start, size_type limit);

  size_type size() const; // Expensive.
  inline bool empty() const { return (this->size() == 0); }
  size_type runs() const; // Expensive.
  inline size_type outdegree() const { return this->outgoing.size(); }

  // Returns (node, LF(i, node)) or invalid_edge() if the offset is invalid.
  edge_type LF(size_type i) const;

  // Returns invalid_offset() if there is no edge to the destination.
  size_type LF(size_type i, node_type to) const;

  // Returns Range::empty_range() if the range is empty or the destination is invalid.
  range_type LF(range_type range, node_type to) const;

  // Returns BWT[i] within the record.
  node_type operator[](size_type i) const;

  // Maps successor nodes to outranks.
  rank_type edgeTo(node_type to) const;

  // These assume that 'outrank' is a valid outgoing edge.
  inline node_type successor(rank_type outrank) const { return this->outgoing[outrank].first; }
  inline size_type offset(rank_type outrank) const { return this->outgoing[outrank].second; }
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

  void swap(RecordArray& another);
  RecordArray& operator=(const RecordArray& source);
  RecordArray& operator=(RecordArray&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  // 0-based indexing.
  inline size_type start(size_type record) const { return this->select(record + 1); }
  inline size_type limit(size_type record) const
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
  sdsl::bit_vector                 sampled_nodes;
  sdsl::bit_vector::rank_1_type    node_rank;

  // Map from node ranks to BWT offsets.
  sdsl::sd_vector<>                bwt_ranges;
  sdsl::sd_vector<>::select_1_type bwt_select;

  // Sampled offsets.
  sdsl::sd_vector<>                sampled_offsets;
  sdsl::sd_vector<>::rank_1_type   sample_rank;

  sdsl::int_vector<0>              samples;

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

  inline size_type size() const { return this->samples.size(); }

private:
  void copy(const DASamples& source);
  void setVectors();
};

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_SUPPORT_H
