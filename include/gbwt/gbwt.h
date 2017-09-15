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

#ifndef GBWT_GBWT_H
#define GBWT_GBWT_H

#include <gbwt/files.h>
#include <gbwt/support.h>

namespace gbwt
{

/*
  gbwt.h: Compressed GBWT structures for construction.
*/

//------------------------------------------------------------------------------

class GBWT
{
public:
  typedef CompressedRecord::size_type size_type;
  typedef node_type                   comp_type; // Index of a record in this->bwt.

//------------------------------------------------------------------------------

  GBWT();
  GBWT(const GBWT& source);
  GBWT(GBWT&& source);
  ~GBWT();

  void swap(GBWT& another);
  GBWT& operator=(const GBWT& source);
  GBWT& operator=(GBWT&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  const static std::string EXTENSION; // .gbwt

//------------------------------------------------------------------------------

  inline size_type size() const { return this->header.size; }
  inline bool empty() const { return (this->size() == 0); }
  inline size_type sequences() const { return this->header.sequences; }
  inline size_type sigma() const { return this->header.alphabet_size; }
  inline size_type effective() const { return this->header.alphabet_size - this->header.offset; }
  inline size_type count(node_type node) const { return this->record(node).size(); }  // Expensive.

  inline bool contains(node_type node) const
  {
    return ((node < this->sigma() && node > this->header.offset) || node == 0);
  }

  inline comp_type toComp(node_type node) const { return (node == 0 ? node : node - this->header.offset); }

  size_type runs() const;
  inline size_type samples() const { return this->da_samples.size(); }

//------------------------------------------------------------------------------

  /*
    The interface assumes that the node identifiers are valid. They can be checked with
    contains().
  */

  // On error: invalid_edge().
  inline edge_type LF(node_type from, size_type i) const
  {
    return this->record(from).LF(i);
  }

  // On error: invalid_edge().
  inline edge_type LF(edge_type position) const
  {
    return this->record(position.first).LF(position.second);
  }

  // On error: invalid_offset().
  inline size_type LF(node_type from, size_type i, node_type to) const
  {
    return this->record(from).LF(i, to);
  }

  // On error: invalid_offset().
  inline size_type LF(edge_type position, node_type to) const
  {
    return this->record(position.first).LF(position.second, to);
  }

  // On error: Range::empty_range().
  inline range_type LF(node_type from, range_type range, node_type to) const
  {
    return this->record(from).LF(range, to);
  }

  // Returns the sampled document identifier or invalid_sequence() if there is no sample.
  inline size_type tryLocate(node_type node, size_type i) const
  {
    return this->da_samples.tryLocate(this->toComp(node), i);
  }

  // Returns the sampled document identifier or invalid_sequence() if there is no sample.
  inline size_type tryLocate(edge_type position) const
  {
    return this->da_samples.tryLocate(this->toComp(position.first), position.second);
  }

  // On error: invalid_sequence().
  size_type locate(node_type node, size_type i) const;
  inline size_type locate(edge_type position) const { return this->locate(position.first, position.second); }

  // On error: empty vector.
  std::vector<size_type> locate(node_type node, range_type range) const;

//------------------------------------------------------------------------------

  // This returns the compressed record for the given node, assuming that it exists.
  CompressedRecord record(node_type node) const;

//------------------------------------------------------------------------------

  GBWTHeader  header;
  RecordArray bwt;
  DASamples   da_samples;

//------------------------------------------------------------------------------

private:
  void copy(const GBWT& source);

//------------------------------------------------------------------------------

}; // class GBWT

void printStatistics(const GBWT& gbwt, const std::string& name);

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_GBWT_H
