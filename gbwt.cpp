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

#include "gbwt.h"

namespace gbwt
{

//------------------------------------------------------------------------------

Sequence::Sequence() :
  id(0), curr(ENDMARKER), next(ENDMARKER), offset(0), pos(0)
{
}

Sequence::Sequence(const text_type& text, size_type i, size_type seq_id) :
  id(seq_id), curr(ENDMARKER), next(text[i]), offset(seq_id), pos(i)
{
}

//------------------------------------------------------------------------------

size_type
DynamicRecord::LF(size_type i, rank_type outrank) const
{
  size_type res = this->offset_in(outrank);
  if(i == 0) { return res; }

  size_type j = 0;
  for(run_type run : this->body)
  {
    if(run.first == outrank) { res += run.second; }
    j += run.second;
    if(j + 1 >= i)
    {
      if(run.first == outrank) { res -= j + 1 - i; }
      break;
    }
  }
  return res;
}

//------------------------------------------------------------------------------

size_type
DynamicGBWT::LF(node_type from, size_type i, node_type to) const
{
  if(to >= this->sigma()) { return this->size(); }
  if(from >= this->sigma()) { return this->bwt[to].size(); }

  const DynamicRecord& from_node = this->bwt[from];
  rank_type outrank = from_node.edge_to(to);
  if(outrank < from_node.outdegree())
  {
    return from_node.LF(i, outrank);
  }

  /*
    Edge (from, to) has not been observed. We find the first edge from a node >= 'from' to 'to'.
    If 'inrank' is equal to indegree, all incoming edges are from nodes < 'from'.
    Otherwise the result is the stored offset in the node we found.
  */
  const DynamicRecord& to_node = this->bwt[to];
  rank_type inrank = to_node.find_first(from);
  if(inrank >= to_node.indegree()) { return to_node.size(); }
  const DynamicRecord& next_from = this->bwt[to_node.predecessor(inrank)];
  return next_from.offset_in(next_from.edge_to(to));
}

//------------------------------------------------------------------------------

void
DynamicGBWT::insert(const text_type& text)
{
  if(text.empty()) { return; }
  if(text[text.size() - 1] != ENDMARKER)
  {
    std::cerr << "DynamicGBWT::insert(): The text must end an endmarker" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  for(node_type node : text)
  {
    if(node >= this->sigma())
    {
      std::cerr << "DynamicGBWT::insert(): Cannot insert " << node << " with alphabet size " << this->sigma() << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  std::vector<Sequence> sequences;
  bool insert = true;
  for(size_type i = 0; i < text.size(); i++)
  {
    if(insert)
    {
      Sequence temp(text, i, this->count(text[0]) + sequences.size());
      sequences.push_back(temp); insert = false;
    }
    if(text[i] == ENDMARKER) { insert = true; }
  }

  // Invariant: Sequences are sorted by (curr, offset).
  // FIXME implement
  while(!(sequences.empty()))
  {
    // for each curr
    //   for each sequence
    //     insert next
    //     set offset to rank(next) within the record
    //     increment the count of curr in the incoming edges of next
    // sort by (next, curr, offset) (maintains the invariant)
    // remove the sequences where next is the endmarker
    // for each next
    //   rebuild the offsets of outgoing edges of each predecessor
    //   for each sequence
    //     add the offset of the outgoing edge (curr, next) to offset
    //     advance to next
  }
}

//------------------------------------------------------------------------------

} // namespace gbwt
