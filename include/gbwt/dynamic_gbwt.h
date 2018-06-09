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

#ifndef GBWT_DYNAMIC_GBWT_H
#define GBWT_DYNAMIC_GBWT_H

#include <thread>

#include <gbwt/gbwt.h>

namespace gbwt
{

/*
  dynamic_gbwt.h: Dynamic GBWT structures for construction.
*/

//------------------------------------------------------------------------------

class GBWTBuilder;

class DynamicGBWT
{
public:
  typedef DynamicRecord::size_type size_type;

  constexpr static size_type INSERT_BATCH_SIZE = 100 * MILLION; // Nodes.
  constexpr static size_type MERGE_BATCH_SIZE = 2000;           // Sequences.
  constexpr static size_type SAMPLE_INTERVAL = 1024;            // Positions in a sequence.

//------------------------------------------------------------------------------

  DynamicGBWT();
  DynamicGBWT(const DynamicGBWT& source);
  DynamicGBWT(DynamicGBWT&& source);
  ~DynamicGBWT();

  void swap(DynamicGBWT& another);
  DynamicGBWT& operator=(const DynamicGBWT& source);
  DynamicGBWT& operator=(DynamicGBWT&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  const static std::string EXTENSION; // .gbwt

//------------------------------------------------------------------------------

  /*
    Insert one or more sequences to the GBWT. The text must be a concatenation of sequences,
    each of which ends with an endmarker (0). The new sequences receive identifiers starting
    from this->sequences().

    Set has_both_orientations to true if the text has both orientations of each sequence.
    Both orientations are required for bidirectional search.
  */
  void insert(const text_type& text, bool has_both_orientations = false);
  void insert(const text_type& text, size_type text_length, bool has_both_orientations = false);
  void insert(const vector_type& text, bool has_both_orientations = false);

  /*
    Use the above to insert the sequences in batches of up to 'batch_size' nodes. Use batch
    size 0 to insert the entire text at once. By default, the sequences are only inserted in
    forward orientation. Set both_orientations = true to insert the reverse complement as well.
    Both orientations are required for bidirectional search.
  */
  void insert(text_buffer_type& text, size_type batch_size = INSERT_BATCH_SIZE, bool both_orientations = false);

  /*
    Insert the sequences from the other GBWT into this. Use batch size 0 to insert all
    sequences at once.
  */
  void merge(const GBWT& source, size_type batch_size = MERGE_BATCH_SIZE);

//------------------------------------------------------------------------------

  /*
    Low-level interface: Statistics.
  */

  size_type size() const { return this->header.size; }
  bool empty() const { return (this->size() == 0); }
  size_type sequences() const { return this->header.sequences; }
  size_type sigma() const { return this->header.alphabet_size; }
  size_type effective() const { return this->header.alphabet_size - this->header.offset; }

  size_type runs() const;     // Expensive.
  size_type samples() const;  // Expensive.

  bool bidirectional() const { return this->header.get(GBWTHeader::FLAG_BIDIRECTIONAL); }

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

  size_type locate(node_type node, size_type i) const { return gbwt::locate(*this, range_type(node, i)); }
  size_type locate(edge_type position) const { return gbwt::locate(*this, position); }

  std::vector<size_type> locate(node_type node, range_type range) const { return this->locate(SearchState(node, range)); }
  std::vector<size_type> locate(SearchState state) const;

  vector_type extract(size_type sequence) const { return gbwt::extract(*this, sequence); }

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
    return ((node < this->sigma() && node > this->header.offset) || node == 0);
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

  comp_type toComp(node_type node) const { return (node == 0 ? node : node - this->header.offset); }
  node_type toNode(comp_type comp) const { return (comp == 0 ? comp : comp + this->header.offset); }

  size_type nodeSize(node_type node) const { return this->record(node).size(); }

  DynamicRecord& record(node_type node)
  {
    return this->bwt[this->toComp(node)];
  }

  const DynamicRecord& record(node_type node) const
  {
    return this->bwt[this->toComp(node)];
  }

//------------------------------------------------------------------------------

  /*
    Low-level interface: Navigation and searching. The interface assumes that node
    identifiers are valid. This can be checked with contains().
  */

  // On error: invalid_edge().
  edge_type LF(node_type from, size_type i) const
  {
    return this->record(from).LF(i);
  }

  // On error: invalid_edge().
  edge_type LF(edge_type position) const
  {
    return this->record(position.first).LF(position.second);
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
  size_type tryLocate(node_type node, size_type i) const;
  size_type tryLocate(edge_type position) const { return this->tryLocate(position.first, position.second); }

//------------------------------------------------------------------------------

  GBWTHeader                 header;
  std::vector<DynamicRecord> bwt;

//------------------------------------------------------------------------------

  // Change offset or alphabet size if the new values are beyond the current values.
  void resize(size_type new_offset, size_type new_sigma);

private:
  void copy(const DynamicGBWT& source);

  /*
    Sort the outgoing edges and change the outranks in the runs accordingly.
    While the GBWT works with any edge order, serialization requires sorted edges,
    as the identifiers of destination nodes are gap-encoded.
  */
  void recode();

  friend class GBWTBuilder;

//------------------------------------------------------------------------------

}; // class DynamicGBWT

void printStatistics(const DynamicGBWT& gbwt, const std::string& name);

//------------------------------------------------------------------------------

class GBWTBuilder
{
public:
  GBWTBuilder(size_type node_width, size_type batch_size = DynamicGBWT::INSERT_BATCH_SIZE);
  ~GBWTBuilder();

  void swapIndex(DynamicGBWT& another_index);

  void insert(const vector_type& sequence, bool both_orientations = false);
  void finish();

  DynamicGBWT index;
  text_type   input_buffer, construction_buffer;
  size_type   input_tail, construction_tail;
  size_type   inserted_sequences, batch_sequences;
  std::thread builder;
  bool        has_both_orientations;

  GBWTBuilder(const GBWTBuilder&) = delete;
  GBWTBuilder& operator= (const GBWTBuilder&) = delete;

private:
  void flush();
}; // class GBWTBuilder

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_DYNAMIC_GBWT_H
