/*
  Copyright (c) 2017, 2018, 2019 Jouni Siren
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

#ifndef GBWT_SUPPORT_H
#define GBWT_SUPPORT_H

#include <gbwt/utils.h>

namespace gbwt
{

/*
  support.h: Public support structures.
*/

//------------------------------------------------------------------------------

/*
  A simple encoding between (node id, orientation) <-> node_type.
*/
struct Node
{
  constexpr static node_type REVERSE_MASK = 0x1;
  constexpr static size_type ID_SHIFT     = 1;

  static size_type id(node_type node) { return (node >> ID_SHIFT); }
  static bool is_reverse(node_type node) { return (node & REVERSE_MASK); }
  static node_type encode(size_type node_id, bool reversed) { return ((node_id << ID_SHIFT) | reversed); }
  static node_type reverse(node_type node) { return (node ^ REVERSE_MASK); }
};

/*
  A simple encoding between (path id, orientation) <-> size_type.
*/
struct Path
{
  constexpr static size_type REVERSE_MASK = 0x1;
  constexpr static size_type ID_SHIFT     = 1;

  static size_type id(size_type path) { return (path >> ID_SHIFT); }
  static bool is_reverse(size_type path) { return (path & REVERSE_MASK); }
  static size_type encode(size_type path_id, bool reversed) { return ((path_id << ID_SHIFT) | reversed); }
  static size_type reverse(size_type path) { return (path ^ REVERSE_MASK); }
};

/*
  Create a path traversing the reverse nodes in reverse order:
    - in place
    - appending it to the to the output vector
    - inserting it to the tail of the output text, updating the tail
*/
void reversePath(vector_type& path);
void reversePath(const vector_type& path, vector_type& output);
void reversePath(const vector_type& path, text_type& output, size_type& tail);

//------------------------------------------------------------------------------

/*
  The part of the BWT corresponding to a single node (the suffixes starting with / the
  prefixes ending with that node).

  - Incoming edges are sorted by the source node.
  - Outgoing edges are sorted by the destination node.
  - Sampled sequence ids are sorted by the offset.
*/

rank_type edgeTo(node_type to, const std::vector<edge_type>& outgoing);

struct DynamicRecord
{
  typedef gbwt::size_type size_type;

  size_type                body_size;
  std::vector<edge_type>   incoming, outgoing;
  std::vector<run_type>    body;
  std::vector<sample_type> ids;

//------------------------------------------------------------------------------

  DynamicRecord();

  size_type size() const { return this->body_size; }
  bool empty() const { return (this->size() == 0); }
  size_type indegree() const { return this->incoming.size(); }
  size_type outdegree() const { return this->outgoing.size(); }
  size_type runs() const { return this->body.size(); }
  size_type samples() const { return this->ids.size(); }

  void clear();
  void swap(DynamicRecord& another);

//------------------------------------------------------------------------------

  // Sort the outgoing edges if they are not sorted.
  void recode();

  // Remove outgoing edges that are not used and recode the body.
  void removeUnusedEdges();

  // Write the compressed representation.
  void writeBWT(std::vector<byte_type>& data) const;

//------------------------------------------------------------------------------

  // Returns (node, LF(i, node)) or invalid_edge() if the offset is invalid.
  edge_type LF(size_type i) const;

  // As above, but also sets 'run_end' to the last offset of the current run.
  edge_type runLF(size_type i, size_type& run_end) const;

  // Returns invalid_offset() if there is no edge to the destination.
  size_type LF(size_type i, node_type to) const;

  // Returns Range::empty_range() if the range is empty or the destination is invalid.
  range_type LF(range_type range, node_type to) const;

  // As above, but also returns the number of characters x with
  // Node::reverse(x) < Node::reverse(to) in the range.
  range_type bdLF(range_type range, node_type to, size_type& reverse_offset) const;

  // Returns BWT[i] within the record.
  node_type operator[](size_type i) const;

//------------------------------------------------------------------------------

  bool hasEdge(node_type to) const;

  // Maps successor nodes to outranks.
  rank_type edgeTo(node_type to) const { return gbwt::edgeTo(to, this->outgoing); }

  // This version works when the edges are not sorted.
  rank_type edgeToLinear(node_type to) const;

  // These assume that 'outrank' is a valid outgoing edge.
  node_type successor(rank_type outrank) const { return this->outgoing[outrank].first; }
#ifdef GBWT_SAVE_MEMORY
  short_type& offset(rank_type outrank) { return this->outgoing[outrank].second; }
#else
  size_type& offset(rank_type outrank) { return this->outgoing[outrank].second; }
#endif
  size_type offset(rank_type outrank) const { return this->outgoing[outrank].second; }

//------------------------------------------------------------------------------

  // These assume that 'inrank' is a valid incoming edge.
  node_type predecessor(rank_type inrank) const { return this->incoming[inrank].first; }
#ifdef GBWT_SAVE_MEMORY
  short_type& count(rank_type inrank) { return this->incoming[inrank].second; }
#else
  size_type& count(rank_type inrank) { return this->incoming[inrank].second; }
#endif
  size_type count(rank_type inrank) const { return this->incoming[inrank].second; }

  // The sum of count(inrank) for all 'inrank' with predecessor(inrank) < 'from'.
  size_type countBefore(node_type from) const;

  // The sum of count(inrank) for all 'inrank' with predecessor(inrank) <= 'from'.
  size_type countUntil(node_type from) const;

  // Increment the count of the incoming edge from 'from'.
  void increment(node_type from);

  // Add a new incoming edge.
  void addIncoming(edge_type inedge);

//------------------------------------------------------------------------------

  // Returns the first sample at offset >= i or ids.end() if there is no sample.
  std::vector<sample_type>::const_iterator nextSample(size_type i) const;

};  // struct DynamicRecord

std::ostream& operator<<(std::ostream& out, const DynamicRecord& record);

//------------------------------------------------------------------------------

struct CompressedRecord
{
  typedef gbwt::size_type size_type;

  std::vector<edge_type> outgoing;
  const byte_type*       body;
  size_type              data_size;

  CompressedRecord();
  CompressedRecord(const std::vector<byte_type>& source, size_type start, size_type limit);

  size_type size() const; // Expensive.
  bool empty() const { return (this->size() == 0); }
  size_type runs() const; // Expensive.
  size_type outdegree() const { return this->outgoing.size(); }

  // Returns (node, LF(i, node)) or invalid_edge() if the offset is invalid.
  edge_type LF(size_type i) const;

  // As above, but also sets 'run_end' to the last offset of the current run.
  edge_type runLF(size_type i, size_type& run_end) const;

  // Returns invalid_offset() if there is no edge to the destination.
  size_type LF(size_type i, node_type to) const;

  // Returns Range::empty_range() if the range is empty or the destination is invalid.
  range_type LF(range_type range, node_type to) const;

  // As above, but also returns the number of characters x with
  // Node::reverse(x) < Node::reverse(to) in the range.
  range_type bdLF(range_type range, node_type to, size_type& reverse_offset) const;

  // Returns BWT[i] within the record.
  node_type operator[](size_type i) const;

  bool hasEdge(node_type to) const;

  // Maps successor nodes to outranks.
  rank_type edgeTo(node_type to) const { return gbwt::edgeTo(to, this->outgoing); };

  // These assume that 'outrank' is a valid outgoing edge.
  node_type successor(rank_type outrank) const { return this->outgoing[outrank].first; }
  size_type offset(rank_type outrank) const { return this->outgoing[outrank].second; }
};

//------------------------------------------------------------------------------

/*
  A record decompressed into an edge array. Good for extracting entire paths,
  but no support for searching with LF(i, to), LF(range, to), or bdLF().
*/
struct DecompressedRecord
{
  typedef gbwt::size_type size_type;

  std::vector<edge_type> outgoing;
  std::vector<edge_type> after; // Outgoing edges after this record.
  std::vector<edge_type> body;

  DecompressedRecord();
  explicit DecompressedRecord(const DynamicRecord& source);
  explicit DecompressedRecord(const CompressedRecord& source);

  size_type size() const { return this->body.size(); }
  bool empty() const { return (this->size() == 0); }
  size_type runs() const; // Expensive.
  size_type outdegree() const { return this->outgoing.size(); }

  // Returns (node, LF(i, node)) or invalid_edge() if the offset is invalid.
  edge_type LF(size_type i) const;

  // As above, but also sets 'run_end' to the last offset of the current run.
  edge_type runLF(size_type i, size_type& run_end) const;

  // Returns BWT[i] within the record.
  node_type operator[](size_type i) const;

  bool hasEdge(node_type to) const;

  // Maps successor nodes to outranks.
  rank_type edgeTo(node_type to) const { return gbwt::edgeTo(to, this->outgoing); };

  // These assume that 'outrank' is a valid outgoing edge.
  node_type successor(rank_type outrank) const { return this->outgoing[outrank].first; }
  size_type offset(rank_type outrank) const { return this->outgoing[outrank].second; }
  size_type offsetAfter(rank_type outrank) const { return this->after[outrank].second; }
};

//------------------------------------------------------------------------------

/*
  An iterator over the 1-bits in sdsl::sd_vector<>.
*/
struct SDIterator
{
  const sdsl::sd_vector<>& vector;

  size_type low_offset, high_offset;
  size_type vector_offset;

  SDIterator(const sdsl::sd_vector<>& v, size_type i) :
    vector(v)
  {
    this->select(i);
  }

  size_type operator*() const { return this->vector_offset; }
  size_type rank() const { return this->low_offset; }
  size_type size() const { return this->vector.low.size(); }
  bool end() const { return (this->rank() >= this->size()); }

  void select(size_type i)
  {
    this->low_offset = i - 1;
    this->high_offset = this->vector.high_1_select(i);
    this->vector_offset = this->vector.low[this->low_offset] + ((this->high_offset + 1 - i) << this->vector.wl);
  }

  void operator++()
  {
    this->low_offset++;
    if(this->end()) { return; }
    do
    {
      this->high_offset++;
    }
    while(this->vector.high[this->high_offset] != 1);
    this->vector_offset = this->vector.low[low_offset] + ((this->high_offset - this->low_offset) << this->vector.wl);
  }
};

//------------------------------------------------------------------------------

struct RecordArray
{
  typedef gbwt::size_type size_type;

  size_type                        records;
  sdsl::sd_vector<>                index;
  sdsl::sd_vector<>::select_1_type select;
  std::vector<byte_type>           data;

  RecordArray();
  RecordArray(const RecordArray& source);
  RecordArray(RecordArray&& source);
  ~RecordArray();

  explicit RecordArray(const std::vector<DynamicRecord>& bwt);
  RecordArray(const std::vector<RecordArray const*> sources, const sdsl::int_vector<0>& origins, const std::vector<size_type>& record_offsets);

  // Set the number of records, build the data manually, and give the offsets to build the index.
  explicit RecordArray(size_type array_size);
  void buildIndex(const std::vector<size_type>& offsets);

  void swap(RecordArray& another);
  RecordArray& operator=(const RecordArray& source);
  RecordArray& operator=(RecordArray&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  size_type size() const { return this->records; }
  bool empty() const { return (this->size() == 0); }

  // 0-based indexing.
  size_type start(size_type record) const { return this->select(record + 1); }
  size_type limit(size_type record) const
  {
    return (record + 1 < this->size() ? this->select(record + 2) : this->data.size());
  }

private:
  void copy(const RecordArray& source);
};

//------------------------------------------------------------------------------

struct DASamples
{
  typedef gbwt::size_type size_type;

  // Does node i have samples?
  sdsl::bit_vector                 sampled_records;
  sdsl::bit_vector::rank_1_type    record_rank;

  // Map from record ranks to BWT offsets.
  sdsl::sd_vector<>                bwt_ranges;
  sdsl::sd_vector<>::select_1_type bwt_select;

  // Sampled offsets.
  sdsl::sd_vector<>                sampled_offsets;
  sdsl::sd_vector<>::rank_1_type   sample_rank;

  sdsl::int_vector<0>              array;

  DASamples();
  DASamples(const DASamples& source);
  DASamples(DASamples&& source);
  ~DASamples();

  explicit DASamples(const std::vector<DynamicRecord>& bwt);
  DASamples(const std::vector<DASamples const*> sources, const sdsl::int_vector<0>& origins, const std::vector<size_type>& record_offsets, const std::vector<size_type>& sequence_counts);

  void swap(DASamples& another);
  DASamples& operator=(const DASamples& source);
  DASamples& operator=(DASamples&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  size_type records() const { return this->sampled_records.size(); }
  size_type size() const { return this->array.size(); }

  // Returns invalid_sequence() if there is no sample.
  size_type tryLocate(size_type record, size_type offset) const;

  // Returns the first sample at >= offset or invalid_sample() if there is no sample.
  sample_type nextSample(size_type record, size_type offset) const;

  bool isSampled(size_type record) const { return this->sampled_records[record]; }

  // We assume that 'record' has samples.
  size_type start(size_type record) const { return this->bwt_select(this->record_rank(record) + 1); }

  // Upper bound for the range of a record, given its rank among records with samples.
  size_type limit(size_type rank) const;

private:
  void copy(const DASamples& source);
  void setVectors();
};

//------------------------------------------------------------------------------

struct MergeParameters
{
  constexpr static size_type POS_BUFFER_SIZE = 64; // Megabytes.
  constexpr static size_type THREAD_BUFFER_SIZE = 256; // Megabytes.
  constexpr static size_type MERGE_BUFFERS = 6;
  constexpr static size_type CHUNK_SIZE = 1; // Sequences per thread.
  constexpr static size_type MERGE_JOBS = 4;

  constexpr static size_type MAX_BUFFER_SIZE = 16384; // Megabytes.
  constexpr static size_type MAX_MERGE_BUFFERS = 16;
  constexpr static size_type MAX_MERGE_JOBS = 16;

  MergeParameters();

  void setPosBufferSize(size_type megabytes);
  void setThreadBufferSize(size_type megabytes);
  void setMergeBuffers(size_type n);
  void setChunkSize(size_type n);
  void setMergeJobs(size_type n);

  // These return the sizes in positions/bytes.
  size_type posBufferPositions() const { return (this->pos_buffer_size * MEGABYTE) / sizeof(edge_type); }
  size_type threadBufferBytes() const { return this->thread_buffer_size * MEGABYTE; }

  size_type pos_buffer_size, thread_buffer_size;
  size_type merge_buffers;
  size_type chunk_size;
  size_type merge_jobs;
};

//------------------------------------------------------------------------------

/*
  GBWT metadata. Merging assumes that all metadata records cover the same samples
  and haplotypes.

  Version 0:
  - Preliminary version with counts of samples, haplotypes, and contigs.

  Next version:
  - Haplotype coverage over ranges of node ids (run-length encode with sd_vector).
  - Assign contig ids to ranges of node ids (as above).
  - Contig and sample names.
  - Sequence names as tuples (sample, phase, count) with implied orientation.
*/

struct Metadata
{
  typedef gbwt::size_type size_type;  // Needed for SDSL serialization.

  std::uint32_t tag;
  std::uint32_t version;
  std::uint64_t sample_count;
  std::uint64_t haplotype_count;
  std::uint64_t contig_count;
  std::uint64_t flags;

  constexpr static std::uint32_t TAG = 0x6B375E7A;
  constexpr static std::uint32_t VERSION = Version::METADATA_VERSION;

  // Flag masks for old compatible versions.

  constexpr static std::uint64_t FLAG_MASK = 0x0000;

  Metadata();

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);
  bool check() const;
  bool checkNew() const;

  void setVersion() { this->version = VERSION; }

  void set(std::uint64_t flag) { this->flags |= flag; }
  void unset(std::uint64_t flag) { this->flags &= ~flag; }
  bool get(std::uint64_t flag) const { return (this->flags & flag); }

  void swap(Metadata& another);

  bool operator==(const Metadata& another) const;
  bool operator!=(const Metadata& another) const { return !(this->operator==(another)); }

  size_type samples() const { return this->sample_count; }
  size_type haplotypes() const { return this->haplotype_count; }
  size_type contigs() const { return this->contig_count; }

  void setSamples(size_type n) { this->sample_count = n; }
  void setHaplotypes(size_type n) { this->haplotype_count = n; }
  void setContigs(size_type n) { this->contig_count = n; }

  void merge(const Metadata& source, bool same_samples, bool same_contigs);
  void merge(std::vector<const Metadata*> sources, bool same_samples, bool same_contigs);
  void clear();
};

std::ostream& operator<<(std::ostream& stream, const Metadata& metadata);

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_SUPPORT_H
