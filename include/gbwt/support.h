/*
  Copyright (c) 2017, 2018, 2019, 2020, 2021, 2025 Jouni Siren
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

#include <functional>
#include <map>

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
  std::pair<size_type, size_type> runs() const; // (concrete, logical)
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

  // Returns `offset` such that `LF(offset) == (to, i)`, or `invalid_offset()`
  // if there is no such offset.
  // This can be used for computing inverse LF in a bidirectional GBWT.
  size_type offsetTo(node_type to, size_type i) const;

  // As above, but also reports the closed offset range ('run') and the identifier
  // ('run_id') of the logical run used for computing LF().
  edge_type LF(size_type i, range_type& run, size_type& run_id) const;

  // As above, but also sets 'run_end' to the last offset of the current logical run.
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

  // Returns the size of the record corresponding to the given semiopen interval.
  static size_type recordSize(const std::vector<byte_type>& source, size_type start, size_type limit);

  // Checks whether the record starting at the given position is empty.
  static bool emptyRecord(const std::vector<byte_type>& source, size_type start);

  size_type size() const; // Expensive.
  bool empty() const { return (this->size() == 0); }
  std::pair<size_type, size_type> runs() const; // (concrete, logical)
  size_type outdegree() const { return this->outgoing.size(); }

  // Returns (node, LF(i, node)) or invalid_edge() if the offset is invalid.
  edge_type LF(size_type i) const;

  // Returns the predecessor node for the sequene at offset `i` in the other orientation
  // of this node, or `invalid_node()` if there is no such node.
  // This can be used for computing inverse LF in a bidirectional GBWT.
  node_type predecessorAt(size_type i) const;

  // Returns `offset` such that `LF(offset) == (to, i)`, or `invalid_offset()`
  // if there is no such offset.
  // This can be used for computing inverse LF in a bidirectional GBWT.
  size_type offsetTo(node_type to, size_type i) const;

  // As above, but also reports the closed offset range ('run') and the identifier
  // ('run_id') of the logical run used for computing LF().
  edge_type LF(size_type i, range_type& run, size_type& run_id) const;

  // As above, but also sets 'run_end' to the last offset of the current logical run.
  edge_type runLF(size_type i, size_type& run_end) const;

  // Returns invalid_offset() if there is no edge to the destination.
  size_type LF(size_type i, node_type to) const;

  // Returns Range::empty_range() if the range is empty or the destination is invalid.
  range_type LF(range_type range, node_type to) const;

  // As above, but also sets 'starts_with_to' if the range starts with node 'to', and
  // sets 'first_run' to the run identifier of the first run of to in overlapping
  // with the range.
  range_type LF(range_type range, node_type to, bool& starts_with_to, size_type& first_run) const;

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
  DecompressedRecord(const DecompressedRecord& source);
  DecompressedRecord(DecompressedRecord&& source);
  ~DecompressedRecord();

  explicit DecompressedRecord(const DynamicRecord& source);
  explicit DecompressedRecord(const CompressedRecord& source);

  void swap(DecompressedRecord& another);
  DecompressedRecord& operator=(const DecompressedRecord& source);
  DecompressedRecord& operator=(DecompressedRecord&& source);

  size_type size() const { return this->body.size(); }
  bool empty() const { return (this->size() == 0); }
  std::pair<size_type, size_type> runs() const; // (concrete, logical)
  size_type outdegree() const { return this->outgoing.size(); }

  // Returns (node, LF(i, node)) or invalid_edge() if the offset is invalid.
  // LF(i, run, run_id) is not supported, because there are no good ways of
  // determining the run_id.
  edge_type LF(size_type i) const;

  // As above, but also sets 'run_end' to the last offset of the current logical run.
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

private:
  void copy(const DecompressedRecord& source);
};

//------------------------------------------------------------------------------

struct RecordArray
{
  typedef gbwt::size_type size_type;

  size_type                        records; // This is redundant, as we could use `index.ones()`.
  sdsl::sd_vector<>                index;
  sdsl::sd_vector<>::select_1_type select;
  std::vector<byte_type>           data;

  RecordArray();
  RecordArray(const RecordArray& source);
  RecordArray(RecordArray&& source);
  ~RecordArray();

  explicit RecordArray(const std::vector<DynamicRecord>& bwt);
  RecordArray(const std::vector<RecordArray const*> sources, const sdsl::int_vector<0>& origins);

  // Set the number of records, build the data manually, and give the offsets to build the index.
  explicit RecordArray(size_type array_size);
  void buildIndex(const std::vector<size_type>& offsets);

  void swap(RecordArray& another);
  RecordArray& operator=(const RecordArray& source);
  RecordArray& operator=(RecordArray&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  void simple_sds_serialize(std::ostream& out) const;
  void simple_sds_load(std::istream& in);
  size_t simple_sds_size() const;

  size_type size() const { return this->records; }
  bool empty() const { return (this->size() == 0); }

  size_type size(size_type record) const;
  bool empty(size_type record) const { return CompressedRecord::emptyRecord(this->data, this->select(record + 1)); }

  // Records use 0-based indexing and semiopen ranges [start, limit).
  std::pair<size_type, size_type> getRange(size_type record) const;
  void forEach(std::function<void(size_type, const CompressedRecord&)> iteratee) const;

private:
  void copy(const RecordArray& source);

  // Throws `sdsl::simple_sds::InvalidData` if the checks fail.
  void sanityChecks() const;
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

  sdsl::int_vector<0>              array;

  DASamples();
  DASamples(const DASamples& source);
  DASamples(DASamples&& source);
  ~DASamples();

  explicit DASamples(const std::vector<DynamicRecord>& bwt);

  // Assumes that the samples are in sorted order.
  explicit DASamples(const RecordArray& bwt, const std::vector<std::pair<node_type, sample_type>>& samples);

  DASamples(const std::vector<DASamples const*> sources, const sdsl::int_vector<0>& origins, const std::vector<size_type>& record_offsets, const std::vector<size_type>& sequence_counts);

  void swap(DASamples& another);
  DASamples& operator=(const DASamples& source);
  DASamples& operator=(DASamples&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  void simple_sds_serialize(std::ostream& out) const;
  void simple_sds_load(std::istream& in);
  size_t simple_sds_size() const;

  size_type records() const { return this->sampled_records.size(); }
  size_type size() const { return this->array.size(); }

  // Returns invalid_sequence() if there is no sample.
  size_type tryLocate(size_type record, size_type offset) const;

  // Returns the first sample at >= offset or invalid_sample() if there is no sample.
  sample_type nextSample(size_type record, size_type offset) const;

  bool isSampled(size_type record) const { return this->sampled_records[record]; }

  // We assume that 'record' has samples.
  size_type start(size_type record) const { return this->bwt_select(this->record_rank(record) + 1); }

private:
  void copy(const DASamples& source);
  void setVectors();

  // Throws `sdsl::simple_sds::InvalidData` if the checks fail.
  void sanityChecks() const;
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
  An array of strings stored in a single character vector, with starting offsets
  stored in an integer vector. This can be serialized and loaded much faster than
  an array of actual strings.
*/
class StringArray
{
public:
  typedef gbwt::size_type size_type;

  StringArray() : index(1, 0, 1) {}
  StringArray(const std::vector<std::string>& source);
  StringArray(const std::map<std::string, std::string>& source);
  StringArray(size_type n, const std::function<size_type(size_type)>& length, const std::function<view_type(size_type)>& sequence);
  StringArray(size_type n, const std::function<bool(size_type)>& choose, const std::function<size_type(size_type)>& length, const std::function<view_type(size_type)>& sequence);
  StringArray(size_type n, const std::function<size_type(size_type)>& length, const std::function<std::string(size_type)>& sequence);

  void swap(StringArray& another);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  void simple_sds_serialize(std::ostream& out) const;
  void simple_sds_load(std::istream& in);
  size_t simple_sds_size() const;

  // This version loads each string twice and transforms the second copy.
  // The transform function should not change the length of the string.
  void simple_sds_load_duplicate(std::istream& in, const std::function<void(std::string&)>& transform);

  bool operator==(const StringArray& another) const;
  bool operator!=(const StringArray& another) const;

  size_type size() const { return this->index.size() - 1; }
  bool empty() const { return (this->size() == 0); }
  size_type length() const { return this->strings.size(); }
  size_type length(size_type i) const { return (this->index[i + 1] - this->index[i]); }
  size_type length(size_type start, size_t limit) const { return (this->index[limit] - this->index[start]); }

  std::string str(size_type i) const
  {
    return std::string(this->strings.data() + this->index[i], this->strings.data() + this->index[i + 1]);
  }

  view_type view(size_type i) const
  {
    return view_type(this->strings.data() + this->index[i], this->length(i));
  }

  void remove(size_type i);

  sdsl::int_vector<0> index;
  std::vector<char>   strings;

private:
  // Throws `sdsl::simple_sds::InvalidData` if the checks fail.
  void sanityChecks() const;
};

//------------------------------------------------------------------------------

class Dictionary
{
public:
  typedef gbwt::size_type size_type;

  StringArray         strings;
  sdsl::int_vector<0> sorted_ids; // String ids in sorted order.

  Dictionary();
  Dictionary(const Dictionary& source);
  Dictionary(Dictionary&& source);
  ~Dictionary();

  explicit Dictionary(const std::vector<std::string>& source);

  // Create a dictionary that contains all strings from the first and all missing strings
  // from the second dictionary.
  Dictionary(const Dictionary& first, const Dictionary& second);

  void swap(Dictionary& another);
  Dictionary& operator=(const Dictionary& source);
  Dictionary& operator=(Dictionary&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  // Load the version of Dictionary used in Metadata version 1.
  void load_v1(std::istream& in);

  void simple_sds_serialize(std::ostream& out) const;
  void simple_sds_load(std::istream& in);
  size_t simple_sds_size() const;

  bool operator==(const Dictionary& another) const;
  bool operator!=(const Dictionary& another) const { return !(this->operator==(another)); }

  void clear();

  size_type size() const { return this->strings.size(); }
  bool empty() const { return (this->size() == 0); }
  size_type length() const { return this->strings.length(); }

  // Return key i or an empty string if there is no such key.
  std::string operator[](size_type i) const 
  {
    return (i >= this->size() ? std::string() : this->strings.str(i));
  }

  // Returns size() if not found.
  size_type find(const std::string& s) const { return this->find(str_to_view(s)); }
  size_type find(view_type view) const;

  // Removes key i.
  void remove(size_type i);

  // The same as `*this = Dictionary(*this, source)`.
  void append(const Dictionary& source);

  bool hasDuplicates() const;

private:
  void copy(const Dictionary& source);

  // Throws `sdsl::simple_sds::InvalidData` if the checks fail.
  void sanityChecks() const;

  void sortKeys();

  // Indexes in sorted_ids.
  bool smaller_by_order(size_type left, size_type right) const;
  bool smaller_by_order(size_type left, view_type right) const;
  bool smaller_by_order(view_type left, size_type right) const;

  // Indexes in offsets.
  bool smaller_by_id(size_type left, size_type right) const;
  bool smaller_by_id(size_type left, view_type right) const;
  bool smaller_by_id(view_type left, size_type right) const;
};

//------------------------------------------------------------------------------

class Tags
{
public:
  typedef gbwt::size_type size_type;

  std::map<std::string, std::string> tags;

  Tags();
  Tags(const Tags& source);
  Tags(Tags&& source);
  ~Tags();

  void swap(Tags& another);
  Tags& operator=(const Tags& source);
  Tags& operator=(Tags&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  void simple_sds_serialize(std::ostream& out) const;
  void simple_sds_load(std::istream& in);
  size_t simple_sds_size() const;

  bool operator==(const Tags& another) const;
  bool operator!=(const Tags& another) const { return !(this->operator==(another)); }

  void set(const std::string& key, const std::string& value);

  // Returns an empty string if the key does not exist.
  std::string get(const std::string& key) const;

  bool contains(const std::string& key) const;

  void clear();

  size_type size() const { return this->tags.size(); }
  bool empty() const { return (this->size() == 0); }

private:
  void copy(const Tags& source);
  void build(const StringArray& source);
  static std::string normalize(const std::string& key);
};

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_SUPPORT_H
