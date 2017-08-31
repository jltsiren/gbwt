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

#ifndef GBWT_DYNAMIC_GBWT_H
#define GBWT_DYNAMIC_GBWT_H

#include "gbwt.h"

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
  typedef node_type                comp_type; // Index of a record in this->bwt.

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

  inline size_type size() const { return this->header.size; }
  inline bool empty() const { return (this->size() == 0); }
  inline size_type sequences() const { return this->header.sequences; }
  inline size_type sigma() const { return this->header.alphabet_size; }
  inline size_type effective() const { return this->header.alphabet_size - this->header.offset; }
  inline size_type count(node_type node) const { return this->record(node).size(); }

  size_type runs() const;
  size_type samples() const;

//------------------------------------------------------------------------------

  // Returns invalid_offset() if the destination is invalid.
  size_type LF(node_type from, size_type i, node_type to) const;

  // Returns invalid_edge() if the node or the offset is invalid.
  edge_type LF(node_type from, size_type i) const;

  // Returns Range::empty_range() if the range is empty or the destination is invalid.
  range_type LF(node_type from, range_type range, node_type to) const;

//------------------------------------------------------------------------------

  /*
    These functions get the BWT record for the given node. Because the alphabet is empty
    in range [1..offset], we cannot simply access the BWT.
  */

  inline DynamicRecord& record(node_type node)
  {
    return this->bwt[node == 0 ? node : node - this->header.offset];
  }

  inline const DynamicRecord& record(node_type node) const
  {
    return this->bwt[node == 0 ? node : node - this->header.offset];
  }

//------------------------------------------------------------------------------

  GBWTHeader                 header;
  std::vector<DynamicRecord> bwt;

//------------------------------------------------------------------------------

private:
  void copy(const DynamicGBWT& source);
  void resize(size_type new_offset, size_type new_sigma); // Change offset or alphabet size.

  /*
    Sort the outgoing edges and change the outranks in the runs accordingly.
    While the GBWT works with any edge order, serialization requires sorted edges,
    as the identifiers of destination nodes are gap-encoded.
  */
  void recode();

  /*
    Insert a batch of sequences with ids (in the current input) starting from 'start_id'.
  */
  void insertBatch(const text_type& text, size_type start_id = 0);

//------------------------------------------------------------------------------

}; // class DynamicGBWT

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_DYNAMIC_GBWT_H
