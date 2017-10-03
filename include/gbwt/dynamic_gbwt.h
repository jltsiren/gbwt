/*
  Copyright (c) 2017 Jouni Siren
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

#include <gbwt/gbwt.h>

namespace gbwt
{

/*
  dynamic_gbwt.h: Dynamic GBWT structures for construction.
*/

//------------------------------------------------------------------------------

class DynamicGBWT
{
public:
  typedef DynamicRecord::size_type size_type;

  const static size_type INSERT_BATCH_SIZE = 100 * MILLION; // Nodes.
  const static size_type MERGE_BATCH_SIZE = 2000;           // Sequences.
  const static size_type SAMPLE_INTERVAL = 1024;            // Positions in a sequence.

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
  */
  void insert(const text_type& text);
  void insert(const std::vector<node_type>& text);

  /*
    Use the above to insert the sequences in batches of up to 'batch_size' nodes. Use batch
    size 0 to insert the entire text at once.
  */
  void insert(text_buffer_type& text, size_type batch_size = INSERT_BATCH_SIZE);

  /*
    Insert the sequences from the other GBWT into this. Use batch size 0 to insert all
    sequences at once.

    FIXME Special case when the node ids do not overlap.
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

  std::vector<node_type> extract(size_type sequence) const { return gbwt::extract(*this, sequence); }

//------------------------------------------------------------------------------

  /*
    Low-level interface: Nodes. The interface assumes that node identifiers are valid.
    This can be checked with contains().
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

  comp_type toComp(node_type node) const { return (node == 0 ? node : node - this->header.offset); }

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

//------------------------------------------------------------------------------

}; // class DynamicGBWT

void printStatistics(const DynamicGBWT& gbwt, const std::string& name);

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_DYNAMIC_GBWT_H
