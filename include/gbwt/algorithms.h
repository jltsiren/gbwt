/*
  Copyright (c) 2017, 2018 Jouni Siren
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

#ifndef GBWT_ALGORITHMS_H
#define GBWT_ALGORITHMS_H

#include "support.h"

namespace gbwt
{

/*
  algorithms.h: High-level search algorithms.
*/

//------------------------------------------------------------------------------

struct SearchState
{
  node_type  node;
  range_type range;

  SearchState() : node(ENDMARKER), range(Range::empty_range()) {}
  SearchState(node_type node_id, range_type offset_range) : node(node_id), range(offset_range) {}
  SearchState(node_type node_id, size_type sp, size_type ep) : node(node_id), range(sp, ep) {}

  size_type size() const { return Range::length(this->range); }
  bool empty() const { return Range::empty(this->range); }

  bool operator==(SearchState another) const { return (this->node == another.node && this->range == another.range); }
  bool operator!=(SearchState another) const { return (this->node != another.node || this->range != another.range); }
};

std::ostream& operator<<(std::ostream& out, SearchState state);

//------------------------------------------------------------------------------

struct BidirectionalState
{
  SearchState forward, backward;

  BidirectionalState() {}
  BidirectionalState(SearchState forward_state, SearchState backward_state) : forward(forward_state), backward(backward_state) {}

  size_type size() const { return this->forward.size(); }
  bool empty() const { return this->forward.empty(); }

  void flip() { std::swap(this->forward, this->backward); }

  bool operator==(BidirectionalState another) const { return (this->forward == another.forward && this->backward == another.backward); }
  bool operator!=(BidirectionalState another) const { return (this->forward != another.forward || this->backward != another.backward); }
};

std::ostream& operator<<(std::ostream& out, BidirectionalState state);

//------------------------------------------------------------------------------

/*
  If the parameters are invalid or if there are no matches, the search algorithms return an
  empty SearchState.

  Template parameters:
    GBWTType  GBWT or DynamicGBWT
    Iterator  InputIterator
*/

template<class GBWTType>
SearchState
extend(const GBWTType& index, SearchState state, node_type node)
{
  if(state.empty() || !(index.contains(node))) { return SearchState(); }
  state.range = index.LF(state, node);
  state.node = node;
  return state;
}

template<class GBWTType, class Iterator>
SearchState
extend(const GBWTType& index, SearchState state, Iterator begin, Iterator end)
{
  while(begin != end && !(state.empty()))
  {
    state = gbwt::extend(index, state, *begin);
    ++begin;
  }
  return state;
}

template<class GBWTType>
SearchState
find(const GBWTType& index, node_type node)
{
  if(!(index.contains(node))) { return SearchState(); }
  return SearchState(node, 0, index.nodeSize(node) - 1);
}

template<class GBWTType, class Iterator>
SearchState
find(const GBWTType& index, Iterator begin, Iterator end)
{
  if(begin == end) { return SearchState(); }

  SearchState state = gbwt::find(index, *begin);
  ++begin;

  return gbwt::extend(index, state, begin, end);
}

template<class GBWTType>
SearchState
prefix(const GBWTType& index, node_type node)
{
  SearchState state(ENDMARKER, 0, index.sequences() - 1);
  return gbwt::extend(index, state, node);
}

template<class GBWTType, class Iterator>
SearchState
prefix(const GBWTType& index, Iterator begin, Iterator end)
{
  SearchState state(ENDMARKER, 0, index.sequences() - 1);
  return gbwt::extend(index, state, begin, end);
}

//------------------------------------------------------------------------------

/*
  If the parameters are invalid or if there are no matches, the search algorithms return an
  empty BidirectionalState. Bidirectional search expects an index with both orientations of
  each sequence.

  Template parameters:
    GBWTType  GBWT or DynamicGBWT
*/

template<class GBWTType>
BidirectionalState
bdExtendForward(const GBWTType& index, BidirectionalState state, node_type node)
{
  if(state.empty() || !(index.contains(node))) { return BidirectionalState(); }
  size_type reverse_offset = 0;
  state.forward.range = index.bdLF(state.forward, node, reverse_offset);
  state.forward.node = node;
  state.backward.range.first += reverse_offset;
  state.backward.range.second = state.backward.range.first + state.forward.size() - 1;
  return state;
}

template<class GBWTType>
BidirectionalState
bdExtendBackward(const GBWTType& index, BidirectionalState state, node_type node)
{
  state.flip();
  state = bdExtendForward(index, state, Node::reverse(node));
  state.flip();
  return state;
}

template<class GBWTType>
BidirectionalState
bdFind(const GBWTType& index, node_type node)
{
  if(!(index.contains(node))) { return BidirectionalState(); }
  SearchState forward(node, 0, index.nodeSize(node) - 1);
  SearchState backward(Node::reverse(node), forward.range);
  return BidirectionalState(forward, backward);
}

//------------------------------------------------------------------------------

/*
  If the parameters are invalid, the locate algorithms return invalid_sequence() or an
  empty vector.

  Template parameters:
    GBWTType  GBWT or DynamicGBWT
*/

template<class GBWTType>
size_type
locate(const GBWTType& index, edge_type position)
{
  if(!(index.contains(position))) { return invalid_sequence(); }

  // No need to check for invalid_edge(), if the initial position is valid.
  while(true)
  {
    size_type result = index.tryLocate(position);
    if(result != invalid_sequence()) { return result; }
    position = index.LF(position);
  }
}

//------------------------------------------------------------------------------

/*
  If the parameters are invalid, the extraction algorithms return an empty container.

  Template parameters:
    GBWTType  GBWT or DynamicGBWT
*/

template<class GBWTType>
vector_type
extract(const GBWTType& index, edge_type position, size_type max_length)
{
  vector_type result;
  if(position == invalid_edge() || !(index.contains(position))) { return result; }

  // No need to check for invalid_edge(), if the initial position is valid.
  while(position.first != ENDMARKER && result.size() < max_length)
  {
    result.push_back(position.first);
    position = index.LF(position);
  }
  return result;
}

template<class GBWTType>
vector_type
extract(const GBWTType& index, edge_type position)
{
  vector_type result;
  if(position == invalid_edge() || !(index.contains(position))) { return result; }

  // No need to check for invalid_edge(), if the initial position is valid.
  while(position.first != ENDMARKER)
  {
    result.push_back(position.first);
    position = index.LF(position);
  }
  return result;
}

template<class GBWTType>
vector_type
extract(const GBWTType& index, size_type sequence)
{
  if(sequence >= index.sequences()) { return vector_type(); }
  return extract(index, index.start(sequence));
}

//------------------------------------------------------------------------------

/*
  Returns DA samples for the given sample interval in sorted order.

  Template parameters:
    GBWTType  GBWT or DynamicGBWT
*/

template<class GBWTType>
std::vector<std::pair<comp_type, sample_type>>
resample(const GBWTType& index, size_type sample_interval)
{
  if(sample_interval == 0) { sample_interval = std::numeric_limits<size_type>::max(); }
  std::vector<std::pair<comp_type, sample_type>> result;

  //#pragma omp parallel for schedule(dynamic, 1)
  for(size_type sequence = 0; sequence < index.sequences(); sequence++)
  {
    std::vector<edge_type> buffer;
    edge_type curr(ENDMARKER, sequence);
    size_type distance = 0; // From the initial endmarker to the next position.
    do
    {
      edge_type next = index.LF(curr); distance++;
      if(distance % sample_interval == 0 || next.first == ENDMARKER) { buffer.push_back(curr); }
      curr = next;
    }
    while(curr.first != ENDMARKER);
    //#pragma omp critical
    {
      result.reserve(result.size() + buffer.size());
      for(edge_type pos : buffer)
      {
        result.emplace_back(index.toComp(pos.first), sample_type(pos.second, sequence));
      }
    }
  }

  parallelQuickSort(result.begin(), result.end());
  return result;
}

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_ALGORITHMS_H

