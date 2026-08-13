/*
  Copyright (c) 2017, 2018, 2019, 2021 Jouni Siren
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

#ifndef GBWT_CACHED_GBWT_H
#define GBWT_CACHED_GBWT_H

#include "gbwt.h"

namespace gbwt
{

/*
  cached_gbwt.h: Record caching for compressed GBWT. Use a separate CachedGBWT object
  for each thread, as the object is not thread-safe.
*/

//------------------------------------------------------------------------------

class CachedGBWT
{
public:
  typedef GBWT::size_type size_type;

  constexpr static size_type INITIAL_CAPACITY = 256;
  constexpr static size_type SINGLE_CAPACITY  = 2;
  constexpr static double    MAX_LOAD_FACTOR  = 0.77;

//------------------------------------------------------------------------------

  CachedGBWT();
  CachedGBWT(const CachedGBWT& source);
  CachedGBWT(CachedGBWT&& source);
  ~CachedGBWT();

  // Set single_record = true to quickly cache a single record.
  explicit CachedGBWT(const GBWT& gbwt_index, bool single_record = false);

  void swap(CachedGBWT& another);
  CachedGBWT& operator=(const CachedGBWT& source);
  CachedGBWT& operator=(CachedGBWT&& source);

//------------------------------------------------------------------------------

  /*
    Cache interface. Primarily used for accessing the successors of a node.
  */

  size_type cacheSize() const { return this->cached_records.size(); }
  size_type cacheCapacity() const { return this->cache_index.capacity(); }
  void clearCache();

  // Insert the record into the cache if it is not already there. Return the cache offset of the record.
  // Note: This assumes that the node does exist. Use contains() to check.
  size_type findRecord(node_type node) const;

  // Return the outdegree of the cached node.
  size_type outdegree(size_type cache_offset) const { return this->cached_records[cache_offset].outdegree(); }

  // Return the i-th successor node of the cached node.
  node_type successor(size_type cache_offset, size_type i) const { return this->cached_records[cache_offset].successor(i); }

  // Extend the state forward to the i-th successor of the cached node (state.node).
  SearchState cachedExtend(SearchState state, size_type cache_offset, size_type i) const;

  // Extend the state forward to the i-th successor of the cached node (state.forward.node).
  BidirectionalState cachedExtendForward(BidirectionalState state, size_type cache_offset, size_type i) const;

  // Extend the state backward to the i-th successor of the cached node (state.backward.node).
  BidirectionalState cachedExtendBackward(BidirectionalState state, size_type cache_offset, size_type i) const;

//------------------------------------------------------------------------------

  /*
    Low-level interface: Statistics.

    Note: These are simple wrappers.
  */

  size_type size() const { return this->index->size(); }
  bool empty() const { return this->index->empty(); }
  size_type sequences() const { return this->index->sequences(); }
  size_type sigma() const { return this->index->sigma(); }
  size_type effective() const { return this->index->effective(); }

  std::pair<size_type, size_type> runs() const { return this->index->runs(); } // Not cached.
  size_type samples() const { return this->index->samples(); }

  bool bidirectional() const { return this->index->bidirectional(); }

//------------------------------------------------------------------------------

  /*
    High-level interface. The queries check that the parameters are valid. Iterators
    must be InputIterators. On error or failed search, the return values will be the
    following:

    find     empty search state
    prefix   empty search state
    extend   empty search state
    locate   invalid_sequence() or empty vector
    extract  empty vector

    Note: These are all cached.
  */

  SearchState find(node_type node) const { return gbwt::find(*this, node); }

  template<class Iterator>
  SearchState find(Iterator begin, Iterator end) const { return gbwt::find(*this, begin, end); }

  SearchState prefix(node_type node) const { return gbwt::prefix(*this, node); }

  template<class Iterator>
  SearchState prefix(Iterator begin, Iterator end) const { return gbwt::prefix(*this, begin, end); }

  SearchState extend(SearchState state, node_type node) const { return gbwt::extend(*this, state, node); }

  template<class Iterator>
  SearchState extend(SearchState state, Iterator begin, Iterator end) const { return gbwt::extend(*this, state, begin, end); }

  size_type locate(node_type node, size_type i) const { return gbwt::locate(*this, edge_type(node, i)); }
  size_type locate(edge_type position) const { return gbwt::locate(*this, position); }

  std::vector<size_type> locate(node_type node, range_type range) const { return this->locate(SearchState(node, range)); }
  std::vector<size_type> locate(SearchState state) const;

  vector_type extract(size_type sequence) const { return gbwt::extract(*this, sequence); }
  vector_type extract(edge_type position) const { return gbwt::extract(*this, position); }
  vector_type extract(edge_type position, size_type max_length) const { return gbwt::extract(*this, position, max_length); }

//------------------------------------------------------------------------------

  /*
    Bidirectional search interface. The queries check that the parameters are valid.
    On error or failed search, the return value is an empty bidirectional search state.

    Note: These are cached.
  */

  BidirectionalState bdFind(node_type node) const { return gbwt::bdFind(*this, node); }

  BidirectionalState bdExtendForward(BidirectionalState state, node_type node) const { return gbwt::bdExtendForward(*this, state, node); }

  BidirectionalState bdExtendBackward(BidirectionalState state, node_type node) const { return gbwt::bdExtendBackward(*this, state, node); }

//------------------------------------------------------------------------------

  /*
    Low-level interface: Nodes. The interface assumes that node identifiers are valid,
    except in contains() / hasEdge(). This can be checked with contains().

    Note: These are cached whenever they access records and simple wrappers otherwise.
  */

  bool contains(node_type node) const { return this->index->contains(node); }

  bool contains(edge_type position) const
  {
    return (this->contains(position.first) && position.second < this->nodeSize(position.first));
  }

  bool contains(SearchState state) const
  {
    return (this->contains(state.node) && !(state.empty()) && state.range.second < this->nodeSize(state.node));
  }

  bool hasEdge(node_type from, node_type to) const
  {
    return (this->contains(from) && this->record(from).hasEdge(to));
  }

  std::vector<edge_type> edges(node_type from) const
  {
    return this->record(from).outgoing;
  }

  node_type firstNode() const { return this->index->firstNode(); }
  comp_type toComp(node_type node) const { return this->index->toComp(node); }
  node_type toNode(comp_type comp) const { return this->index->toNode(comp); }

  size_type nodeSize(node_type node) const { return this->record(node).size(); }
  bool empty(node_type node) const { return this->index->empty(node); }

//------------------------------------------------------------------------------

  /*
    Low-level interface: Navigation and searching. The interface assumes that node
    identifiers are valid. This can be checked with contains().

    Note: These are cached.
  */

  // On error: invalid_edge().
  edge_type LF(node_type from, size_type i) const
  {
    if(from == ENDMARKER) { return this->endmarker().LF(i); }
    return this->record(from).LF(i);
  }

  // On error: invalid_edge().
  edge_type LF(edge_type position) const
  {
    if(position.first == ENDMARKER) { return this->endmarker().LF(position.second); }
    return this->record(position.first).LF(position.second);
  }

  // Only works in bidirectional indexes. May be slow when the predecessor is the endmarker.
  // On error: invalid_edge().
  edge_type inverseLF(node_type from, size_type i) const;

  // Only works in bidirectional indexes. May be slow when the predecessor is the endmarker.
  // On error: invalid_edge().
  edge_type inverseLF(edge_type position) const
  {
    return this->inverseLF(position.first, position.second);
  }

  // On error: invalid_offset().
  size_type LF(node_type from, size_type i, node_type to) const
  {
    return this->record(from).LF(i, to);
  }

  // On error: invalid_offset().
  size_type LF(edge_type position, node_type to) const
  {
    return this->record(position.first).LF(position.second, to);
  }

  // On error: Range::empty_range().
  range_type LF(node_type from, range_type range, node_type to) const
  {
    return this->record(from).LF(range, to);
  }

  // On error: Range::empty_range().
  range_type LF(SearchState state, node_type to) const
  {
    return this->record(state.node).LF(state.range, to);
  }

  // On error: Range::empty_range().
  range_type bdLF(SearchState state, node_type to, size_type& reverse_offset) const
  {
    return this->record(state.node).bdLF(state.range, to, reverse_offset);
  }

//------------------------------------------------------------------------------

  /*
    Low-level interface: Sequences. The interface assumes that node identifiers are
    valid. This can be checked with contains().

    Note: These are simple wrappers.
  */

  // Starting position of the sequence or invalid_edge() if something fails.
  edge_type start(size_type sequence) const { return this->index->start(sequence); }

  // Returns the sampled document identifier or invalid_sequence() if there is no sample.
  size_type tryLocate(node_type node, size_type i) const { return this->index->tryLocate(node, i); }

  // Returns the sampled document identifier or invalid_sequence() if there is no sample.
  size_type tryLocate(edge_type position) const { return this->index->tryLocate(position); }

//------------------------------------------------------------------------------

  const GBWT* index;

  // Node node_in_cache[i].first is at cached_records[node_in_cache[i].second].
  // NOTE: We want to update the cache in const member functions.
  mutable std::vector<edge_type>        cache_index;
  mutable std::vector<CompressedRecord> cached_records;

//------------------------------------------------------------------------------

/*
  Internal interface. Do not use.
*/

private:
  void copy(const CachedGBWT& source);
  size_type indexOffset(node_type node) const;
  void rehash() const;

public:
  // The reference may be invalid after accessing other records.
  const CompressedRecord& record(node_type node) const { return this->cached_records[this->findRecord(node)]; }

  const DecompressedRecord& endmarker() const { return this->index->endmarker(); }
}; // class CachedGBWT

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_CACHED_GBWT_H
