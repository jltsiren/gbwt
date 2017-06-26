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

#ifndef _GBWT_GBWT_H
#define _GBWT_GBWT_H

#include "support.h"

namespace gbwt
{

/*
  gbwt.h: Main GBWT structures.
*/

//------------------------------------------------------------------------------

struct Sequence
{
  size_type id;
  node_type curr, next;
  size_type offset; // Offset in the current record.
  size_type pos;    // Position in the data.

  Sequence();
  Sequence(const text_type& text, size_type i, size_type seq_id);
};

//------------------------------------------------------------------------------

struct DynamicRecord
{
  typedef gbwt::size_type                 size_type;
  typedef Run::value_type                 rank_type;  // Rank of incoming/outgoing edge.
  typedef Run::run_type                   run_type;
  typedef std::pair<node_type, size_type> edge_type;

  size_type              body_size;
  std::vector<edge_type> incoming, outgoing;
  std::vector<run_type>  body;

  inline size_type size() const { return this->body_size; }
  inline size_type runs() const { return this->body.size(); }
  inline size_type indegree() const { return this->incoming.size(); }
  inline size_type outdegree() const { return this->outgoing.size(); }

  // Map global alphabet to local alphabet.
  inline rank_type edge_to(node_type to) const
  {
    for(rank_type outrank = 0; outrank < this->outdegree(); outrank++)
    {
      if(this->outgoing[outrank].first == to) { return outrank; }
    }
    return this->outdegree();
  }

  // This assumes that 'outrank' is a valid outgoing edge.
  inline size_type offset_in(rank_type outrank) const
  {
    return this->outgoing[outrank].second;
  }

  // This assumes that 'outrank' is a valid outgoing edge.
  size_type LF(size_type i, rank_type outrank) const;

  // Return the first node >= 'from' with an incoming edge to this node.
  inline rank_type find_first(node_type from) const
  {
    for(size_type i = 0; i < this->indegree(); i++)
    {
      if(this->incoming[i].first >= from) { return i; }
    }
    return this->indegree();
  }

  // This assumes that 'inrank' is a valid incoming edge.
  inline node_type predecessor(rank_type inrank) const
  {
    return this->incoming[inrank].first;
  }
};

//------------------------------------------------------------------------------

class DynamicGBWT
{
  typedef DynamicRecord::size_type  size_type;
  typedef DynamicRecord::rank_type  rank_type;
  typedef DynamicRecord::run_type   run_type;
  typedef DynamicRecord::edge_type  edge_type;

  inline size_type size() const { return this->bwt_size; }
  inline size_type sigma() const { return this->bwt.size(); }
  inline size_type count(node_type node) const { return this->bwt[node].size(); }

  void insert(const text_type& text);

  size_type LF(node_type from, size_type i, node_type to) const;

  size_type                  bwt_size;
  std::vector<DynamicRecord> bwt;
}; // class DynamicGBWT

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // _GBWT_GBWT_H
