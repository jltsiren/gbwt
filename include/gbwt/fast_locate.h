/*
  Copyright (c) 2020, 2022 Jouni Siren

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

#ifndef GBWT_FAST_LOCATE_H
#define GBWT_FAST_LOCATE_H

#include "gbwt.h"

namespace gbwt
{

/*
  fast_locate.h: A support structure implementing the r-index locate().
*/

//------------------------------------------------------------------------------

/*
  An optional locate() structure based on the r-index. If the source GBWT is changed
  in any way, this structure must be rebuilt. The implementation is based on the
  simplified version in:

    Gagie, Navarro, and Prezza: Fully Functional Suffix Trees and Optimal Text
    Searching in BWT-Runs Bounded Space. Journal of the ACM, 2020.

  We store samples at the beginning of each run instead of at the end, and derive
  SA[i+1] from SA[i]. Because GBWT is a multi-string BWT, each endmarker is
  considered a separate run. And because we are really indexing the reverse texts,
  sequence offsets tell the distance to the end of the text.

  The implementation assumes that the product of the number of sequences and the length
  of the longest sequence fits into 64 bits.

  Version 1:
  - First version.
*/

class FastLocate
{
public:
  typedef gbwt::size_type size_type; // Needed for SDSL serialization.

  constexpr static size_type NO_POSITION = std::numeric_limits<size_type>::max();
  
//------------------------------------------------------------------------------

  FastLocate();
  FastLocate(const FastLocate& source);
  FastLocate(FastLocate&& source);
  ~FastLocate();

  explicit FastLocate(const GBWT& source);

  void swap(FastLocate& another);
  FastLocate& operator=(const FastLocate& source);
  FastLocate& operator=(FastLocate&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);
  void setGBWT(const GBWT& source) { this->index = &source; }

  const static std::string EXTENSION; // .ri

//------------------------------------------------------------------------------

  struct Header
  {
    std::uint32_t tag;
    std::uint32_t version;
    std::uint64_t max_length; // Length of the longest sequence.
    std::uint64_t flags;

    constexpr static std::uint32_t TAG = 0x6B3741D8;
    constexpr static std::uint32_t VERSION = Version::R_INDEX_VERSION;

    Header();

    size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
    void load(std::istream& in);

    // Throws `sdsl::simple_sds::InvalidData` if the header is invalid.
    void check() const;

    void setVersion() { this->version = VERSION; }
    void set(std::uint64_t flag) { this->flags |= flag; }
    void unset(std::uint64_t flag) { this->flags &= ~flag; }
    bool get(std::uint64_t flag) const { return (this->flags & flag); }
  };

//------------------------------------------------------------------------------

  // Source GBWT.
  const GBWT* index;

  Header header;

  // (sequence id, sequence offset) samples at the start of each run.
  sdsl::int_vector<0> samples;

  // Mark the text positions at the end of a run.
  sdsl::sd_vector<> last;

  // If last[i] = 1, last_to_run[last_rank(i)] is the identifier of the run.
  sdsl::int_vector<0> last_to_run;

  // Run identifier of the first run in each node.
  sdsl::int_vector<0> comp_to_run;

//------------------------------------------------------------------------------

  /*
    Low-level interface: Statistics.
  */

  size_type size() const { return this->samples.size(); }
  bool empty() const { return (this->size() == 0); }

//------------------------------------------------------------------------------

  /*
    High-level interface. The queries check that the parameters are valid. Iterators
    must be InputIterators. On error or failed search, the return value will be empty.
    If the state is non-empty, first will be updated to the packed position corresponding
    to the first occurrence in the range.
  */

  SearchState find(node_type node, size_type& first) const;

  template<class Iterator>
  SearchState find(Iterator begin, Iterator end, size_type& first) const;

  SearchState extend(SearchState state, node_type node, size_type& first) const;

  template<class Iterator>
  SearchState extend(SearchState state, Iterator begin, Iterator end, size_type& first) const;

  std::vector<size_type> locate(SearchState state, size_type first = NO_POSITION) const;

  std::vector<size_type> locate(node_type node, range_type range, size_type first = NO_POSITION) const
  {
    return this->locate(SearchState(node, range), first);
  }

  std::vector<size_type> decompressSA(node_type node) const;

  std::vector<size_type> decompressDA(node_type node) const;

//------------------------------------------------------------------------------

  /*
    Low-level interface. The interface assumes that the arguments are valid. This
    be checked with index->contains(node) and seq_id < index->sequences(). There is
    no check for the offset.
  */

  size_type pack(size_type seq_id, size_type seq_offset) const
  {
    return seq_id * this->header.max_length + seq_offset;
  }

  size_type seqId(size_type offset) const { return offset / this->header.max_length; }
  size_type seqOffset(size_type offset) const { return offset % this->header.max_length; }

  std::pair<size_type, size_type> unpack(size_type offset) const
  {
    return std::make_pair(this->seqId(offset), this->seqOffset(offset));
  }

  size_type locateFirst(node_type node) const
  {
    return this->getSample(node, 0);
  }

  size_type locateNext(size_type prev) const;

//------------------------------------------------------------------------------

  /*
    Internal interface. Do not use.
  */

private:
  void copy(const FastLocate& source);

  size_type globalRunId(node_type node, size_type run_id) const
  {
    return this->comp_to_run[this->index->toComp(node)] + run_id;
  }

  size_type getSample(node_type node, size_type run_id) const
  {
    return this->samples[this->globalRunId(node, run_id)];
  }
}; // class FastLocate

//------------------------------------------------------------------------------

/*
  Template query implementations.
*/

template<class Iterator>
SearchState
FastLocate::find(Iterator begin, Iterator end, size_type& first) const
{
  if(begin == end) { return SearchState(); }

  SearchState state = this->find(*begin, first);
  ++begin;

  return this->extend(state, begin, end, first);
}

template<class Iterator>
SearchState
FastLocate::extend(SearchState state, Iterator begin, Iterator end, size_type& first) const
{
  while(begin != end && !(state.empty()))
  {
    state = this->extend(state, *begin, first);
    ++begin;
  }
  return state;
}

//------------------------------------------------------------------------------

void printStatistics(const FastLocate& index, const std::string& name);
std::string indexType(const FastLocate&);

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_FAST_LOCATE_H
