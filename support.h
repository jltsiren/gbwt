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

#ifndef _GBWT_SUPPORT_H
#define _GBWT_SUPPORT_H

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
  - It is not necessary to have the outgoing edges sorted, but we will sort them anyway
  after each insert.
*/

struct DynamicRecord
{
  typedef gbwt::size_type size_type;

  size_type              body_size;
  std::vector<edge_type> incoming, outgoing;
  std::vector<run_type>  body;

//------------------------------------------------------------------------------

  DynamicRecord() : body_size(0) {}

  inline size_type size() const { return this->body_size; }
  inline bool empty() const { return (this->size() == 0); }
  inline size_type runs() const { return this->body.size(); }
  inline size_type indegree() const { return this->incoming.size(); }
  inline size_type outdegree() const { return this->outgoing.size(); }

  bool operator==(const DynamicRecord& another) const;
  inline bool operator!=(const DynamicRecord& another) const { return !(this->operator==(another)); }

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
};

//------------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& out, const DynamicRecord& record);

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // _GBWT_SUPPORT_H
