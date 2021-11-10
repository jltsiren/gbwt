/*
  Copyright (c) 2017, 2018, 2019, 2020, 2021 Jouni Siren
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

#ifndef GBWT_GBWT_H
#define GBWT_GBWT_H

#include <gbwt/algorithms.h>
#include <gbwt/files.h>
#include <gbwt/metadata.h>

namespace gbwt
{

/*
  gbwt.h: Compressed GBWT structures for construction.
*/

//------------------------------------------------------------------------------

// Forward declaration.
class DynamicGBWT;

class GBWT
{
public:
  typedef CompressedRecord::size_type size_type;

//------------------------------------------------------------------------------

  GBWT();
  GBWT(const GBWT& source);
  GBWT(const DynamicGBWT& source);
  GBWT(GBWT&& source);
  ~GBWT();

  // Merge the sources, assuming that node ids do not overlap.
  // Also merges the metadata if all indexes contain it.
  explicit GBWT(const std::vector<GBWT>& sources);

  void swap(GBWT& another);
  GBWT& operator=(const GBWT& source);
  GBWT& operator=(const DynamicGBWT& source);
  GBWT& operator=(GBWT&& source);

  void resample(size_type sample_interval);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  void simple_sds_serialize(std::ostream& out) const;
  void simple_sds_load(std::istream& in);
  size_t simple_sds_size() const;

  const static std::string EXTENSION; // .gbwt

//------------------------------------------------------------------------------

  /*
    Low-level interface: Statistics.
  */

  size_type size() const { return this->header.size; }
  bool empty() const { return (this->size() == 0); }
  size_type sequences() const { return this->header.sequences; }
  size_type sigma() const { return this->header.alphabet_size; }
  size_type effective() const { return this->header.alphabet_size - this->header.offset; }

  std::pair<size_type, size_type> runs() const;
  size_type samples() const { return this->da_samples.size(); }

  bool bidirectional() const { return this->header.get(GBWTHeader::FLAG_BIDIRECTIONAL); }

//------------------------------------------------------------------------------

  /*
    Metadata interface.
  */

  bool hasMetadata() const { return this->header.get(GBWTHeader::FLAG_METADATA); }
  void addMetadata() { this->header.set(GBWTHeader::FLAG_METADATA); }
  void clearMetadata() { this->metadata.clear(); this->header.unset(GBWTHeader::FLAG_METADATA); };

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
  */

  BidirectionalState bdFind(node_type node) const { return gbwt::bdFind(*this, node); }

  BidirectionalState bdExtendForward(BidirectionalState state, node_type node) const { return gbwt::bdExtendForward(*this, state, node); }

  BidirectionalState bdExtendBackward(BidirectionalState state, node_type node) const { return gbwt::bdExtendBackward(*this, state, node); }

//------------------------------------------------------------------------------

  /*
    Low-level interface: Nodes. The interface assumes that node identifiers are valid,
    except in contains() / hasEdge(). This can be checked with contains().
  */

  bool contains(node_type node) const
  {
    return ((node < this->sigma() && node > this->header.offset) || node == ENDMARKER);
  }

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

  node_type firstNode() const { return this->header.offset + 1; }
  comp_type toComp(node_type node) const { return (node == 0 ? node : node - this->header.offset); }
  node_type toNode(comp_type comp) const { return (comp == 0 ? comp : comp + this->header.offset); }

  size_type nodeSize(node_type node) const { return this->bwt.size(this->toComp(node)); }
  bool empty(node_type node) const { return this->bwt.empty(this->toComp(node)); }

//------------------------------------------------------------------------------

  /*
    Low-level interface: Navigation and searching. The interface assumes that node
    identifiers are valid. This can be checked with contains().
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
  */

  // Starting position of the sequence or invalid_edge() if something fails.
  edge_type start(size_type sequence) const { return this->LF(ENDMARKER, sequence); }

  // Returns the sampled document identifier or invalid_sequence() if there is no sample.
  size_type tryLocate(node_type node, size_type i) const
  {
    return this->da_samples.tryLocate(this->toComp(node), i);
  }

  // Returns the sampled document identifier or invalid_sequence() if there is no sample.
  size_type tryLocate(edge_type position) const
  {
    return this->da_samples.tryLocate(this->toComp(position.first), position.second);
  }

//------------------------------------------------------------------------------

  GBWTHeader  header;
  Tags        tags;
  RecordArray bwt;
  DASamples   da_samples;
  Metadata    metadata;

  // Decompress and cache the endmarker, because decompressing it is expensive.
  DecompressedRecord endmarker_record;

//------------------------------------------------------------------------------

/*
  Internal interface. Do not use.
*/

private:
  void copy(const GBWT& source);
  void resetTags();
  void addSource();
  void cacheEndmarker();

public:
  CompressedRecord record(node_type node) const;
  const DecompressedRecord& endmarker() const { return this->endmarker_record; }
}; // class GBWT

//------------------------------------------------------------------------------

void printStatistics(const GBWT& gbwt, const std::string& name, std::ostream& out = std::cout);
std::string indexType(const GBWT&);

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_GBWT_H
