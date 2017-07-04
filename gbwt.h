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

#include "files.h"
#include "support.h"

namespace gbwt
{

/*
  gbwt.h: Main GBWT structures.

  FIXME Reorganize:
    - rename old support.h to internal.h
    - move Sequence to internal.h, DynamicRecord to support.h
    - rename this file to dynamic_gbwt.h

  FIXME We currently assume that the alphabet is dense. In a multi-chromosome graph, the
  alphabet for an individual chromosome will be dense in range [a,b].
*/

//------------------------------------------------------------------------------

struct Sequence
{
  size_type id;
  node_type curr, next;
  size_type offset; // Offset in the current record.
  size_type pos;    // Position in the text.

  Sequence();
  Sequence(const text_type& text, size_type i, size_type seq_id);

  // Sort by reverse prefixes text[..pos+1].
  inline bool operator<(const Sequence& another) const
  {
    if(this->next != another.next) { return (this->next < another.next); }
    if(this->curr != another.curr) { return (this->curr < another.curr); }
    return (this->offset < another.offset);
  }

  // Do not call if 'next' is an endmarker.
  inline void advance(const text_type& text)
  {
    this->curr = this->next;
    this->pos++;
    this->next = text[this->pos];
  }
};

//------------------------------------------------------------------------------

/*
  The part of the BWT corresponding to a single node (the suffixes starting with / the
  prefixes ending with that node). Incoming edges are sorted by the source node, while
  outgoing edges are not sorted.

  FIXME Implement recoding to make outgoing edges also sorted.
*/

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
  inline rank_type edgeTo(node_type to) const
  {
    for(rank_type outrank = 0; outrank < this->outdegree(); outrank++)
    {
      if(this->outgoing[outrank].first == to) { return outrank; }
    }
    return this->outdegree();
  }

  // These assume that 'outrank' is a valid outgoing edge.
  inline node_type successor(rank_type outrank) const { return this->outgoing[outrank].first; }
  inline size_type& offset(rank_type outrank) { return this->outgoing[outrank].second; }
  inline size_type offset(rank_type outrank) const { return this->outgoing[outrank].second; }

  // This assumes that 'outrank' is a valid outgoing edge.
  size_type LF(size_type i, rank_type outrank) const;

  // Return the first node >= 'from' with an incoming edge to this node.
  inline rank_type findFirst(node_type from) const
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

  // Increment the count of the incoming edge from 'from'.
  inline void increment(node_type from)
  {
    for(size_type i = 0; i < this->indegree(); i++)
    {
      if(this->incoming[i].first == from) { this->incoming[i].second++; return; }
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

class DynamicGBWT
{
public:
  typedef DynamicRecord::size_type size_type;
  typedef DynamicRecord::rank_type rank_type;
  typedef DynamicRecord::run_type  run_type;
  typedef DynamicRecord::edge_type edge_type;

//------------------------------------------------------------------------------

  DynamicGBWT();
  DynamicGBWT(const DynamicGBWT& source);
  DynamicGBWT(DynamicGBWT&& source);
  ~DynamicGBWT();

  void swap(DynamicGBWT& another);
  DynamicGBWT& operator=(const DynamicGBWT& source);
  DynamicGBWT& operator=(DynamicGBWT&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);  // FIXME not tested

  const static std::string EXTENSION; // .gbwt

//------------------------------------------------------------------------------

  /*
    Insert one or more sequences to the GBWT. The text must be a concatenation of sequences,
    each of which ends with an endmarker (0). The new sequences receive identifiers starting
    from this->sequences().
  */
  void insert(const text_type& text);

//------------------------------------------------------------------------------

  inline size_type size() const { return this->header.size; }
  inline size_type sequences() const { return this->header.sequences; }
  inline size_type sigma() const { return this->header.alphabet_size; }
  inline size_type count(node_type node) const { return this->bwt[node].size(); }

//------------------------------------------------------------------------------

  size_type LF(node_type from, size_type i, node_type to) const;

//------------------------------------------------------------------------------

  GBWTHeader                 header;
  std::vector<DynamicRecord> bwt;

//------------------------------------------------------------------------------

private:
  void copy(const DynamicGBWT& source);

//------------------------------------------------------------------------------

}; // class DynamicGBWT

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // _GBWT_GBWT_H
