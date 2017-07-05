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

  inline size_type size() const { return this->body_size; }
  inline size_type runs() const { return this->body.size(); }
  inline size_type indegree() const { return this->incoming.size(); }
  inline size_type outdegree() const { return this->outgoing.size(); }

//------------------------------------------------------------------------------

  // Sort the outgoing edges if they are not sorted.
  void recode();

  // This assumes that 'outrank' is a valid outgoing edge.
  size_type LF(size_type i, rank_type outrank) const;

//------------------------------------------------------------------------------

  // Map nodes to outranks.
  inline rank_type edgeTo(node_type to) const
  {
    for(rank_type outrank = 0; outrank < this->outdegree(); outrank++)
    {
      if(this->successor(outrank) == to) { return outrank; }
    }
    return this->outdegree();
  }

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
  inline rank_type findFirst(node_type from) const
  {
    for(size_type i = 0; i < this->indegree(); i++)
    {
      if(this->incoming[i].first >= from) { return i; }
    }
    return this->indegree();
  }

  // These assume that 'inrank' is a valid incoming edge.
  inline node_type predecessor(rank_type inrank) const { return this->incoming[inrank].first; }
#ifdef GBWT_SAVE_MEMORY
  inline short_type& count(rank_type inrank) { return this->incoming[inrank].second; }
#else
  inline size_type& count(rank_type inrank) { return this->incoming[inrank].second; }
#endif
  inline size_type count(rank_type inrank) const { return this->incoming[inrank].second; }

  // Increment the count of the incoming edge from 'from'.
  inline void increment(node_type from)
  {
    for(rank_type inrank = 0; inrank < this->indegree(); inrank++)
    {
      if(this->predecessor(inrank) == from) { this->count(inrank)++; return; }
    }
    this->addIncoming(edge_type(from, 1));
  }

  // Add a new incoming edge.
  inline void addIncoming(edge_type inedge)
  {
    this->incoming.push_back(inedge);
    sequentialSort(this->incoming.begin(), this->incoming.end());
  }
};

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // _GBWT_SUPPORT_H
