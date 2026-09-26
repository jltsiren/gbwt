#ifndef GBWT_SEQUENCE_LOCATE_H
#define GBWT_SEQUENCE_LOCATE_H

#include <gbwt/gbwt.h>

#include <functional>

namespace gbwt
{

/*
  sequence_locate.h: An r-index locate() structure with offsets in base pairs.
*/

//------------------------------------------------------------------------------

/*
  A variant of FastLocate where sequence offsets are measured in base pairs instead
  of nodes. Node lengths come from a user-provided function, which must return the
  same length (>= 1) for both orientations of a node. If the source GBWT or the node
  lengths change, this structure must be rebuilt.

  The sequence offset of a visit is the distance in bp from the end of the node to
  the end of the sequence. The last node has offset 0, and the endmarker has offset
  equal to the sequence length. A visit to node v with offset o covers forward
  positions [len - o - length(v), len - o), where len is the sequence length.

  The r-index invariants hold with bp offsets, because sequences in the same run
  traverse identical nodes after each LF step.

  Version 1:
  - First version.
*/

class SequenceLocate
{
public:
  typedef gbwt::size_type size_type; // Needed for SDSL serialization.
  typedef std::function<size_type(node_type)> length_function;

  constexpr static size_type NO_POSITION = std::numeric_limits<size_type>::max();

//------------------------------------------------------------------------------

  SequenceLocate();
  SequenceLocate(const SequenceLocate& source);
  SequenceLocate(SequenceLocate&& source) noexcept;
  ~SequenceLocate();

  SequenceLocate(const GBWT& source, const length_function& node_length);

  void swap(SequenceLocate& another) noexcept;
  SequenceLocate& operator=(const SequenceLocate& source);
  SequenceLocate& operator=(SequenceLocate&& source) noexcept;

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  // Sets the source GBWT and node lengths; must be called after load().
  void setGBWT(const GBWT& source, const length_function& node_length)
  {
    this->index = &source;
    this->length = node_length;
  }

  const static std::string EXTENSION; // .sri

//------------------------------------------------------------------------------

  struct Header
  {
    std::uint32_t tag;
    std::uint32_t version;
    std::uint64_t max_length; // Length of the longest sequence in bp + 1.
    std::uint64_t flags;

    constexpr static std::uint32_t TAG = 0x5E9B10C8;
    constexpr static std::uint32_t VERSION = Version::SEQUENCE_LOCATE_VERSION;

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

  // Source GBWT and node lengths.
  const GBWT* index;
  length_function length;

  Header header;

  // (sequence id, sequence offset) samples at the start of each run.
  sdsl::int_vector<0> samples;

  // Mark the text positions at the end of a run.
  sdsl::sd_vector<> last;

  // If last[i] = 1, last_to_run[last_rank(i)] is the identifier of the run.
  sdsl::int_vector<0> last_to_run;

  // Run identifier of the first run in each node.
  sdsl::int_vector<0> comp_to_run;

  // Length of each sequence in bp.
  sdsl::int_vector<0> sequence_length;

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

  // Returns the distinct sequence ids in the range.
  std::vector<size_type> locate(SearchState state, size_type first = NO_POSITION) const;

  std::vector<size_type> locate(node_type node, range_type range, size_type first = NO_POSITION) const
  {
    return this->locate(SearchState(node, range), first);
  }

  // Returns (sequence id, sequence offset) for each occurrence in the range, in BWT order.
  std::vector<std::pair<size_type, size_type>> locatePositions(SearchState state, size_type first = NO_POSITION) const;

  // As locatePositions(), but returns (sequence id, forward start of state.node).
  std::vector<std::pair<size_type, size_type>> locateForward(SearchState state, size_type first = NO_POSITION) const;

  /*
    Returns (node, offset in node) for the base at forward position `forward_bp`
    in sequence `seq_id`, or invalid_edge() if the position is invalid.

    We start from the nearest tail sample at or before the target on the same
    sequence (the successor in `last`, as offsets decrease along the sequence),
    or from the endmarker if there is none. Then we follow LF until we reach
    the node containing the target.
  */
  edge_type locateSequence(size_type seq_id, size_type forward_bp) const;

  // As above, but also sets `gbwt_pos` to the GBWT position (node, record offset) of the visit, for use with LF.
  edge_type locateSequence(size_type seq_id, size_type forward_bp, edge_type& gbwt_pos) const;

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

  // Length of the sequence in bp.
  size_type sequenceLength(size_type seq_id) const { return this->sequence_length[seq_id]; }

  // Forward start of a visit to `node` with the given sequence offset.
  size_type forwardStart(size_type seq_id, size_type seq_offset, node_type node) const
  {
    return this->sequenceLength(seq_id) - seq_offset - this->length(node);
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
  void copy(const SequenceLocate& source);

  size_type globalRunId(node_type node, size_type run_id) const
  {
    return this->comp_to_run[this->index->toComp(node)] + run_id;
  }

  size_type getSample(node_type node, size_type run_id) const
  {
    return this->samples[this->globalRunId(node, run_id)];
  }

  // Packed position of the first occurrence in the range.
  size_type firstPosition(SearchState state, size_type first) const;

  // GBWT position at the end of the given global run.
  edge_type runTail(size_type run) const;
}; // class SequenceLocate

//------------------------------------------------------------------------------

/*
  Template query implementations.
*/

template<class Iterator>
SearchState
SequenceLocate::find(Iterator begin, Iterator end, size_type& first) const
{
  if(begin == end) { return SearchState(); }

  SearchState state = this->find(*begin, first);
  ++begin;

  return this->extend(state, begin, end, first);
}

template<class Iterator>
SearchState
SequenceLocate::extend(SearchState state, Iterator begin, Iterator end, size_type& first) const
{
  while(begin != end && !(state.empty()))
  {
    state = this->extend(state, *begin, first);
    ++begin;
  }
  return state;
}

//------------------------------------------------------------------------------

void printStatistics(const SequenceLocate& index, const std::string& name);
std::string indexType(const SequenceLocate&);

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_SEQUENCE_LOCATE_H
