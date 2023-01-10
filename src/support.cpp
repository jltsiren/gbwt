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

#include <gbwt/internal.h>

#include <cctype>
#include <unordered_set>

namespace gbwt {

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr node_type Node::REVERSE_MASK;
constexpr size_type Node::ID_SHIFT;

constexpr node_type Path::REVERSE_MASK;
constexpr size_type Path::ID_SHIFT;

constexpr size_type MergeParameters::POS_BUFFER_SIZE;
constexpr size_type MergeParameters::THREAD_BUFFER_SIZE;
constexpr size_type MergeParameters::MERGE_BUFFERS;
constexpr size_type MergeParameters::CHUNK_SIZE;
constexpr size_type MergeParameters::MERGE_JOBS;
constexpr size_type MergeParameters::MAX_BUFFER_SIZE;
constexpr size_type MergeParameters::MAX_MERGE_BUFFERS;
constexpr size_type MergeParameters::MAX_MERGE_JOBS;

//------------------------------------------------------------------------------

void reversePath(vector_type &path) {
  std::reverse(path.begin(), path.end());
  for (auto &node : path) {
    node = Node::reverse(node);
  }
}

void reversePath(const vector_type &path, vector_type &output) {
  for (auto iter = path.rbegin(); iter != path.rend(); ++iter) {
    output.push_back(Node::reverse(*iter));
  }
}

void reversePath(const vector_type &path, text_type &output, size_type &tail) {
  for (auto iter = path.rbegin(); iter != path.rend(); ++iter) {
    output[tail] = Node::reverse(*iter);
    tail++;
  }
}

//------------------------------------------------------------------------------

rank_type edgeTo(node_type to, const std::vector<edge_type> &outgoing) {
  rank_type low = 0, high = outgoing.size();
  while (low < high) {
    rank_type mid = low + (high - low) / 2;
    if (outgoing[mid].first == to) {
      return mid;
    }
    if (outgoing[mid].first > to) {
      high = mid;
    } else {
      low = mid + 1;
    }
  }
  return outgoing.size();
}

DynamicRecord::DynamicRecord() : body_size(0) {}

DynamicRecord::DynamicRecord(const DynamicRecord &source) {
  body_size = source.body_size;
  incoming.assign(source.incoming.begin(), source.incoming.end());
  outgoing.assign(source.outgoing.begin(), source.outgoing.end());
  body.assign(source.body.begin(), source.body.end());
  ids.assign(source.ids.begin(), source.ids.end());
  //
  outgoing_offset_map = source.outgoing_offset_map;
  incoming_offset_map = source.incoming_offset_map;
  sample_incoming_offset_map = source.sample_incoming_offset_map;
}

DynamicRecord &DynamicRecord::operator=(const DynamicRecord &source) {
  body_size = source.body_size;
  incoming.assign(source.incoming.begin(), source.incoming.end());
  outgoing.assign(source.outgoing.begin(), source.outgoing.end());
  body.assign(source.body.begin(), source.body.end());
  ids.assign(source.ids.begin(), source.ids.end());
  //
  outgoing_offset_map = source.outgoing_offset_map;
  incoming_offset_map = source.incoming_offset_map;
  sample_incoming_offset_map = source.sample_incoming_offset_map;
  return *this;
}

void DynamicRecord::clear() {
  DynamicRecord temp;
  this->swap(temp);
}

void DynamicRecord::swap(DynamicRecord &another) {
  if (this != &another) {
    std::swap(this->body_size, another.body_size);
    this->incoming.swap(another.incoming);
    this->outgoing.swap(another.outgoing);
    this->body.swap(another.body);
    this->ids.swap(another.ids);
    //
    this->outgoing_offset_map.swap(another.outgoing_offset_map);
    this->incoming_offset_map.swap(another.incoming_offset_map);
    this->sample_incoming_offset_map.swap(another.sample_incoming_offset_map);
  }
}

std::pair<size_type, size_type> DynamicRecord::runs() const {
  std::pair<size_type, size_type> result(this->body.size(), 0);
  for (run_type run : this->body) {
    result.second += (this->successor(run.first) == ENDMARKER ? run.second : 1);
  }
  return result;
}

//------------------------------------------------------------------------------

std::uint64_t DynamicRecord::getBodyOffset(std::uint64_t &income_id) {
  std::uint64_t accumulate = 0;
  std::uint64_t body_offset = 0;
  for (auto &e:this->incoming) {
    if (this->incoming_offset_map.count(e.first)==0) {
      income_id = e.first;
      this->incoming_offset_map[income_id] = 0;
      break;
    } else if ((e.second-this->incoming_offset_map[e.first])>0) {
      income_id = e.first;
      accumulate += this->incoming_offset_map[e.first];
      break;
    } else {
      accumulate += e.second;
    }
  }
  body_offset = accumulate;
  return body_offset;
}

void DynamicRecord::updateBodyOffset(const std::uint64_t &income_id) {
  ++this->incoming_offset_map[income_id];
}

std::uint64_t DynamicRecord::getSampleOffset(const size_type &income_id) {
  if (this->sample_incoming_offset_map.count(income_id)==0) {
    std::uint64_t accumulate = 0;
    for (auto &e:incoming) {
      if (e.first>income_id)
        break;
      else {
        if (this->sample_incoming_offset_map.count(e.first)!=0)
          accumulate+=this->sample_incoming_offset_map[e.first];
      }
    } this->sample_incoming_offset_map[income_id] = accumulate;
  }
  return this->sample_incoming_offset_map[income_id];
}

void DynamicRecord::updateSampleOffset(const size_type &income_id) {
  ++this->sample_incoming_offset_map[income_id];
  for(auto &e:incoming) {
    if (e.first>income_id && this->sample_incoming_offset_map.count(e.first)!=0)
      ++sample_incoming_offset_map[e.first];
  }
}

//------------------------------------------------------------------------------

void DynamicRecord::recode() {
  if (this->empty()) {
    return;
  }

  bool sorted = true;
  for (rank_type outrank = 1; outrank < this->outdegree(); outrank++) {
    if (this->successor(outrank) < this->successor(outrank - 1)) {
      sorted = false;
      break;
    }
  }
  if (sorted) {
    return;
  }

  for (run_type &run : this->body) {
    run.first = this->successor(run.first);
  }
  sequentialSort(this->outgoing.begin(), this->outgoing.end());
  for (run_type &run : this->body) {
    run.first = this->edgeTo(run.first);
  }
}

void DynamicRecord::removeUnusedEdges() {
  // Determine which edges are used and replace the outranks with node
  // identifiers.
  std::vector<bool> used(this->outdegree(), false);
  for (run_type &run : this->body) {
    used[run.first] = true;
    run.first = this->successor(run.first);
  }

  // Remove unused edges.
  size_type tail = 0;
  for (size_type i = 0; i < this->outdegree(); i++) {
    this->outgoing[tail] = this->outgoing[i];
    if (used[i]) {
      tail++;
    }
  }
  this->outgoing.resize(tail);

  // Recode the body.
  for (run_type &run : this->body) {
    run.first = this->edgeTo(run.first);
  }
}

void DynamicRecord::writeBWT(std::vector<byte_type> &data) const {
  // Write the outgoing edges.
  ByteCode::write(data, this->outdegree());
  node_type prev = 0;
  for (edge_type outedge : this->outgoing) {
    ByteCode::write(data, outedge.first - prev);
    prev = outedge.first;
    ByteCode::write(data, outedge.second);
  }

  // Write the body.
  if (this->outdegree() > 0) {
    Run encoder(this->outdegree());
    for (run_type run : this->body) {
      encoder.write(data, run);
    }
  }
}
//------------------------------------------------------------------------------

template <class Array>
edge_type LFLoop(Array &result, const std::vector<edge_type> &body, size_type i,
                 range_type &run, size_type &run_id) {
  rank_type last_edge = 0;
  size_type offset = 0;
  size_type runs_seen = 0;
  for (size_type j = 0; j < body.size(); j++) {
    run_type curr = body[j];
    last_edge = curr.first;
    result[curr.first].second += curr.second;
    offset += curr.second;
    runs_seen += (result[curr.first].first == ENDMARKER ? curr.second : 1);
    if (offset > i) {
      if (result[curr.first].first == ENDMARKER) {
        run.first = run.second = i;
        run_id = runs_seen - (offset - i);
      } else {
        run.first = offset - curr.second;
        run.second = offset - 1;
        run_id = runs_seen - 1;
      }
      break;
    }
  }

  result[last_edge].second -= (offset - i);
  return result[last_edge];
}

edge_type DynamicRecord::LF(size_type i) const {
  range_type run(0, 0);
  size_type run_id = 0;
  return this->LF(i, run, run_id);
}

size_type DynamicRecord::offsetTo(node_type to, size_type i) const {
  size_type outrank = this->edgeTo(to);
  if (outrank >= this->outdegree() || this->offset(outrank) > i) {
    return invalid_offset();
  }

  // Find the occurrence of `to` of rank `i`.
  size_type rank = this->offset(outrank), offset = 0;
  for (run_type run : this->body) {
    offset += run.second;
    if (run.first != outrank) {
      continue;
    }
    rank += run.second;
    if (rank > i) {
      return offset - (rank - i);
    }
  }
  return invalid_offset();
}

edge_type DynamicRecord::LF(size_type i, range_type &run,
                            size_type &run_id) const {
  if (i >= this->size()) {
    return invalid_edge();
  }

  if (this->outdegree() <= MAX_OUTDEGREE_FOR_ARRAY) {
    edge_type result[MAX_OUTDEGREE_FOR_ARRAY];
    for (size_type i = 0; i < this->outdegree(); i++) {
      result[i] = this->outgoing[i];
    }
    return LFLoop(result, this->body, i, run, run_id);
  } else {
    std::vector<edge_type> result(this->outgoing);
    return LFLoop(result, this->body, i, run, run_id);
  }
}

edge_type DynamicRecord::runLF(size_type i, size_type &run_end) const {
  range_type run(0, 0);
  size_type run_id = 0;
  edge_type result = this->LF(i, run, run_id);
  run_end = run.second;
  return result;
}

// run is *(--iter); offset and result are for the beginning of the run at iter.
size_type LFLoop(std::vector<run_type>::const_iterator &iter,
                 std::vector<run_type>::const_iterator end, size_type i,
                 rank_type outrank, run_type &run, size_type &offset,
                 size_type &result) {
  while (iter != end && offset < i) {
    run = *iter;
    ++iter;
    offset += run.second;
    if (run.first == outrank) {
      result += run.second;
    }
  }
  return result - (run.first == outrank ? offset - i : 0);
}

size_type DynamicRecord::LF(size_type i, node_type to) const {
  rank_type outrank = this->edgeTo(to);
  if (outrank >= this->outdegree()) {
    return invalid_offset();
  }

  std::vector<run_type>::const_iterator iter = this->body.begin();
  run_type run(0, 0);
  size_type offset = 0, result = this->offset(outrank);

  return LFLoop(iter, this->body.end(), i, outrank, run, offset, result);
}

range_type DynamicRecord::LF(range_type range, node_type to) const {
  if (Range::empty(range)) {
    return Range::empty_range();
  }

  rank_type outrank = this->edgeTo(to);
  if (outrank >= this->outdegree()) {
    return Range::empty_range();
  }

  std::vector<run_type>::const_iterator iter = this->body.begin();
  run_type run(0, 0);
  size_type offset = 0, result = this->offset(outrank);

  // [LF(range.first, to), LF(range.second + 1, to) - 1].
  range.first =
      LFLoop(iter, this->body.end(), range.first, outrank, run, offset, result);
  range.second = LFLoop(iter, this->body.end(), range.second + 1, outrank, run,
                        offset, result) -
                 1;
  return range;
}

range_type DynamicRecord::bdLF(range_type range, node_type to,
                               size_type &reverse_offset) const {
  if (Range::empty(range)) {
    return Range::empty_range();
  }

  rank_type outrank = this->edgeTo(to);
  if (outrank >= this->outdegree()) {
    return Range::empty_range();
  }

  // sp = LF(range.first, to)
  std::vector<run_type>::const_iterator iter = this->body.begin();
  run_type run(0, 0);
  size_type offset = 0, result = this->offset(outrank);
  size_type sp =
      LFLoop(iter, this->body.end(), range.first, outrank, run, offset, result);

  /*
    Count the number of occurrences of nodes x in the query range, where
    Node::reverse(x) < Node::reverse(to), and store it in reverse_offset.

    1. In the easy case, there are no edges to Node::reverse(to), so we only
    compute the occurrences < outrank.
    2. If there are edges to Node::reverse(to) and to is in forward orientation,
    we count the occurrences <= reverse_rank except those of outrank.
    3. If there are edges to Node::reverse(to), and to is in reverse
    orientation, we count the occurrences < reverse_rank < outrank.
  */
  rank_type reverse_rank = this->edgeTo(Node::reverse(to));
  bool subtract_equal = false;
  if (reverse_rank >= this->outdegree()) {
    reverse_rank = outrank;
  } else if (!Node::is_reverse(to)) {
    reverse_rank++;
    subtract_equal = true;
  }

  // Previous run may go past range.first.
  size_type equal = (run.first == outrank ? offset - range.first : 0);
  reverse_offset = (run.first < reverse_rank ? offset - range.first : 0);

  // ep + 1 = LF(range.second + 1, to)
  range.second++;
  while (iter != this->body.end() && offset < range.second) {
    run = *iter;
    ++iter;
    offset += run.second;
    if (run.first == outrank) {
      equal += run.second;
    }
    if (run.first < reverse_rank) {
      reverse_offset += run.second;
    }
  }

  // Last run may go past range.second.
  if (run.first == outrank) {
    equal -= (offset - range.second);
  }
  if (run.first < reverse_rank) {
    reverse_offset -= (offset - range.second);
  }

  if (subtract_equal) {
    reverse_offset -= equal;
  }
  return range_type(sp, sp + equal - 1);
}

node_type DynamicRecord::operator[](size_type i) const {
  if (i >= this->size()) {
    return ENDMARKER;
  }

  size_type offset = 0;
  for (run_type run : this->body) {
    offset += run.second;
    if (offset > i) {
      return this->successor(run.first);
    }
  }

  return ENDMARKER;
}

//------------------------------------------------------------------------------

bool DynamicRecord::hasEdge(node_type to) const {
  for (rank_type outrank = 0; outrank < this->outdegree(); outrank++) {
    if (this->successor(outrank) == to) {
      return true;
    }
  }
  return false;
}

rank_type DynamicRecord::edgeToLinear(node_type to) const {
  for (rank_type outrank = 0; outrank < this->outdegree(); outrank++) {
    if (this->successor(outrank) == to) {
      return outrank;
    }
  }
  return this->outdegree();
}

//------------------------------------------------------------------------------

size_type DynamicRecord::countBefore(node_type from) const {
  size_type result = 0;
  for (rank_type inrank = 0;
       inrank < this->indegree() && this->predecessor(inrank) < from;
       inrank++) {
    result += this->count(inrank);
  }
  return result;
}

size_type DynamicRecord::countUntil(node_type from) const {
  size_type result = 0;
  for (rank_type inrank = 0;
       inrank < this->indegree() && this->predecessor(inrank) <= from;
       inrank++) {
    result += this->count(inrank);
  }
  return result;
}

void DynamicRecord::increment(node_type from) {
  for (rank_type inrank = 0; inrank < this->indegree(); inrank++) {
    if (this->predecessor(inrank) == from) {
      this->count(inrank)++;
      return;
    }
  }
  this->addIncoming(edge_type(from, 1));
}

void DynamicRecord::addIncoming(edge_type inedge) {
  this->incoming.push_back(inedge);
  sequentialSort(this->incoming.begin(), this->incoming.end());
}

//------------------------------------------------------------------------------

std::vector<sample_type>::const_iterator
DynamicRecord::nextSample(size_type i) const {
  std::vector<sample_type>::const_iterator curr = this->ids.begin();
  while (curr != this->ids.end() && curr->first < i) {
    ++curr;
  }
  return curr;
}

//------------------------------------------------------------------------------

std::ostream &operator<<(std::ostream &out, const DynamicRecord &record) {
  std::pair<size_type, size_type> run_counts = record.runs();

  out << "(size " << record.size() << ", " << run_counts.first << " concrete / "
      << run_counts.second << " logical runs, "
      << "indegree " << record.indegree() << ", outdegree "
      << record.outdegree() << ", incoming = " << record.incoming
      << ", outgoing = " << record.outgoing << ", body = " << record.body
      << ", ids = " << record.ids << ")";

  return out;
}

//------------------------------------------------------------------------------

CompressedRecord::CompressedRecord() : outgoing(), body(0), data_size(0) {}

CompressedRecord::CompressedRecord(const std::vector<byte_type> &source,
                                   size_type start, size_type limit) {
  this->outgoing.resize(ByteCode::read(source, start));
  node_type prev = 0;
  for (edge_type &outedge : this->outgoing) {
    outedge.first = ByteCode::read(source, start) + prev;
    prev = outedge.first;
    outedge.second = ByteCode::read(source, start);
  }

  this->body = source.data() + start;
  this->data_size = limit - start;
}

size_type CompressedRecord::recordSize(const std::vector<byte_type> &source,
                                       size_type start, size_type limit) {
  size_type sigma = ByteCode::read(source, start);
  if (sigma == 0) {
    return 0;
  }
  for (size_type i = 0; i < sigma; i++) {
    ByteCode::read(source, start); // Destination node.
    ByteCode::read(source, start); // Rank within the node.
  }

  size_type result = 0;
  Run decoder(sigma);
  while (start < limit) {
    result += decoder.read(source, start).second;
  }

  return result;
}

bool CompressedRecord::emptyRecord(const std::vector<byte_type> &source,
                                   size_type start) {
  return (ByteCode::read(source, start) == 0);
}

size_type CompressedRecord::size() const {
  size_type result = 0;
  if (this->outdegree() > 0) {
    for (CompressedRecordIterator iter(*this); !(iter.end()); ++iter) {
      result += iter->second;
    }
  }
  return result;
}

std::pair<size_type, size_type> CompressedRecord::runs() const {
  std::pair<size_type, size_type> result(0, 0);
  if (this->outdegree() > 0) {
    for (CompressedRecordIterator iter(*this); !(iter.end()); ++iter) {
      result.first++;
      result.second +=
          (this->successor(iter->first) == ENDMARKER ? iter->second : 1);
    }
  }
  return result;
}

edge_type CompressedRecord::LF(size_type i) const {
  range_type run(0, 0);
  size_type run_id = 0;
  return this->LF(i, run, run_id);
}

node_type CompressedRecord::predecessorAt(size_type i) const {
  // 1. Determine the number of sequences going to each successor node.
  std::vector<edge_type> edges = this->outgoing;
  for (edge_type &edge : edges) {
    edge.second = 0;
  }
  for (CompressedRecordIterator iter(*this); !(iter.end()); ++iter) {
    edges[iter->first].second += iter->second;
  }

  // 2. Flip the successors to make them predecessors of the other orientation
  // of this node.
  for (edge_type &edge : edges) {
    if (edge.first != ENDMARKER) {
      edge.first = Node::reverse(edge.first);
    }
  }

  // 3. If the list of predecessors has both orientations of the same node, swap
  // them to make the list sorted again.
  for (size_type rank = 1; rank < edges.size(); rank++) {
    if (Node::id(edges[rank - 1].first) == Node::id(edges[rank].first)) {
      std::swap(edges[rank - 1], edges[rank]);
    }
  }

  // 4. Find the predecessor, if it exists.
  size_type offset = 0;
  for (edge_type edge : edges) {
    offset += edge.second;
    if (offset > i) {
      return edge.first;
    }
  }
  return invalid_node();
}

size_type CompressedRecord::offsetTo(node_type to, size_type i) const {
  size_type outrank = this->edgeTo(to);
  if (outrank >= this->outdegree() || this->offset(outrank) > i) {
    return invalid_offset();
  }

  // Find the occurrence of `to` of rank `i`.
  for (CompressedRecordRankIterator iter(*this, outrank); !(iter.end());
       ++iter) {
    if (iter->first == outrank && iter.rank() > i) {
      return iter.offset() - (iter.rank() - i);
    }
  }
  return invalid_offset();
}

template <class Iterator>
edge_type LFLoop(const CompressedRecord &record, size_type i, range_type &run,
                 size_type &run_id) {
  Iterator iter(record);
  edge_type result = iter.edgeAt(i);
  if (result.first == ENDMARKER) {
    run.first = run.second = i;
    run_id = iter.runId() - (iter.offset() - 1 - i);
  } else {
    run.first = iter.offset() - iter->second;
    run.second = iter.offset() - 1;
    run_id = iter.runId();
  }
  return result;
}

edge_type CompressedRecord::LF(size_type i, range_type &run,
                               size_type &run_id) const {
  if (this->outdegree() == 0) {
    return invalid_edge();
  }

  if (this->outdegree() <= MAX_OUTDEGREE_FOR_ARRAY) {
    return LFLoop<CompressedRecordArrayIterator>(*this, i, run, run_id);
  } else {
    return LFLoop<CompressedRecordFullIterator>(*this, i, run, run_id);
  }
}

edge_type CompressedRecord::runLF(size_type i, size_type &run_end) const {
  range_type run(0, 0);
  size_type run_id = 0;
  edge_type result = this->LF(i, run, run_id);
  run_end = run.second;
  return result;
}

size_type CompressedRecord::LF(size_type i, node_type to) const {
  size_type outrank = this->edgeTo(to);
  if (outrank >= this->outdegree()) {
    return invalid_offset();
  }
  CompressedRecordRankIterator iter(*this, outrank);
  return iter.rankAt(i);
}

range_type CompressedRecord::LF(range_type range, node_type to) const {
  if (Range::empty(range)) {
    return Range::empty_range();
  }

  size_type outrank = this->edgeTo(to);
  if (outrank >= this->outdegree()) {
    return Range::empty_range();
  }
  CompressedRecordRankIterator iter(*this, outrank);
  range.first = iter.rankAt(range.first);
  range.second = iter.rankAt(range.second + 1) - 1;

  return range;
}

range_type CompressedRecord::LF(range_type range, node_type to,
                                bool &starts_with_to,
                                size_type &first_run) const {
  starts_with_to = false;
  first_run = invalid_offset();
  if (Range::empty(range)) {
    return Range::empty_range();
  }

  size_type outrank = this->edgeTo(to);
  if (outrank >= this->outdegree()) {
    return Range::empty_range();
  }
  CompressedRecordRankIterator iter(*this, outrank);

  // Range start. If the run that reaches range.first overlaps with the query
  // range, it may be a run of to.
  size_type start_offset = range.first;
  bool first_run_found =
      false; // Have we seen the first run covering the query range?
  range.first = iter.rankAt(range.first);
  if (iter.offset() > start_offset) {
    first_run_found = true;
    if (iter->first == outrank) {
      starts_with_to = true;
      first_run = iter.runId();
    }
  }

  // Range end. Determine the first run of to, if we have not seen it already.
  while (!(iter.end()) && iter.offset() <= range.second) {
    ++iter;
    if (iter->first == outrank && first_run == invalid_offset()) {
      if (!first_run_found) {
        starts_with_to = true;
      }
      first_run = iter.runId();
    }
    first_run_found = true;
  }
  range.second = iter.rankAt(range.second + 1) - 1;

  if (Range::empty(range)) {
    first_run = invalid_offset();
  }
  return range;
}

range_type CompressedRecord::bdLF(range_type range, node_type to,
                                  size_type &reverse_offset) const {
  if (Range::empty(range)) {
    return Range::empty_range();
  }

  size_type outrank = this->edgeTo(to);
  if (outrank >= this->outdegree()) {
    return Range::empty_range();
  }

  CompressedRecordRankIterator iter(*this, outrank);
  size_type sp = iter.rankAt(range.first);

  /*
    Count the number of occurrences of nodes x in the query range, where
    Node::reverse(x) < Node::reverse(to), and store it in reverse_offset.

    1. In the easy case, there are no edges to Node::reverse(to), so we only
    compute the occurrences < outrank.
    2. If there are edges to Node::reverse(to) and to is in forward orientation,
    we count the occurrences <= reverse_rank except those of outrank.
    3. If there are edges to Node::reverse(to), and to is in reverse
    orientation, we count the occurrences < reverse_rank < outrank.
  */
  rank_type reverse_rank = this->edgeTo(Node::reverse(to));
  if (reverse_rank >= this->outdegree()) {
    reverse_rank = outrank;
  } else if (!Node::is_reverse(to)) {
    reverse_rank++;
  }

  // Previous run may go past range.first.
  if (iter->first < reverse_rank && iter->first != outrank) {
    reverse_offset = iter.offset() - range.first;
  } else {
    reverse_offset = 0;
  }

  range.second++; // We compute rank at range.second + 1.
  while (!(iter.end()) && iter.offset() < range.second) {
    ++iter;
    if (iter->first < reverse_rank && iter->first != outrank) {
      reverse_offset += iter->second;
    }
  }

  // Last run may go past range.second.
  if (iter->first < reverse_rank && iter->first != outrank) {
    reverse_offset -= (iter.offset() - range.second);
  }

  return range_type(sp, iter.rankAt(range.second) - 1);
}

node_type CompressedRecord::operator[](size_type i) const {
  if (this->outdegree() == 0) {
    return ENDMARKER;
  }

  for (CompressedRecordIterator iter(*this); !(iter.end()); ++iter) {
    if (iter.offset() > i) {
      return this->successor(iter->first);
    }
  }
  return ENDMARKER;
}

bool CompressedRecord::hasEdge(node_type to) const {
  for (rank_type outrank = 0; outrank < this->outdegree(); outrank++) {
    if (this->successor(outrank) == to) {
      return true;
    }
  }
  return false;
}

//------------------------------------------------------------------------------

DecompressedRecord::DecompressedRecord() : outgoing(), after(), body() {}

DecompressedRecord::DecompressedRecord(const DecompressedRecord &source) {
  this->copy(source);
}

DecompressedRecord::DecompressedRecord(DecompressedRecord &&source) {
  *this = std::move(source);
}

DecompressedRecord::~DecompressedRecord() {}

DecompressedRecord::DecompressedRecord(const DynamicRecord &source)
    : outgoing(source.outgoing), after(source.outgoing), body() {
  this->body.reserve(source.size());
  for (run_type run : source.body) {
    for (size_type i = 0; i < run.second; i++) {
      this->body.push_back(this->after[run.first]);
      this->after[run.first].second++;
    }
  }
}

DecompressedRecord::DecompressedRecord(const CompressedRecord &source)
    : outgoing(source.outgoing), after(source.outgoing), body() {
  this->body.reserve(source.size());
  for (CompressedRecordIterator iter(source); !(iter.end()); ++iter) {
    for (size_type i = 0; i < iter->second; i++) {
      this->body.push_back(this->after[iter->first]);
      this->after[iter->first].second++;
    }
  }
}

void DecompressedRecord::swap(DecompressedRecord &another) {
  if (this != &another) {
    this->outgoing.swap(another.outgoing);
    this->after.swap(another.after);
    this->body.swap(another.body);
  }
}

DecompressedRecord &DecompressedRecord::
operator=(const DecompressedRecord &source) {
  if (this != &source) {
    this->copy(source);
  }
  return *this;
}

DecompressedRecord &DecompressedRecord::operator=(DecompressedRecord &&source) {
  if (this != &source) {
    this->outgoing = std::move(source.outgoing);
    this->after = std::move(source.after);
    this->body = std::move(source.body);
  }
  return *this;
}

void DecompressedRecord::copy(const DecompressedRecord &source) {
  this->outgoing = source.outgoing;
  this->after = source.after;
  this->body = source.body;
}

std::pair<size_type, size_type> DecompressedRecord::runs() const {
  std::pair<size_type, size_type> result(0, 0);
  if (this->empty()) {
    return result;
  }

  node_type prev = invalid_node();
  for (edge_type edge : body) {
    if (edge.first != prev) {
      result.first++;
      result.second++;
      prev = edge.first;
    } else if (edge.first == ENDMARKER) {
      result.second++;
    }
  }

  return result;
}

edge_type DecompressedRecord::LF(size_type i) const {
  if (i >= this->size()) {
    return invalid_edge();
  }
  return this->body[i];
}

edge_type DecompressedRecord::runLF(size_type i, size_type &run_end) const {
  if (i >= this->size()) {
    return invalid_edge();
  }

  run_end = i;
  if (this->body[i].first != ENDMARKER) {
    while (run_end + 1 < this->size() &&
           this->body[run_end + 1].first == this->body[i].first) {
      run_end++;
    }
  }

  return this->body[i];
}

node_type DecompressedRecord::operator[](size_type i) const {
  if (i >= this->size()) {
    return ENDMARKER;
  }
  return this->body[i].first;
}

bool DecompressedRecord::hasEdge(node_type to) const {
  for (rank_type outrank = 0; outrank < this->outdegree(); outrank++) {
    if (this->successor(outrank) == to) {
      return true;
    }
  }
  return false;
}

//------------------------------------------------------------------------------

RecordArray::RecordArray() : records(0) {}

RecordArray::RecordArray(const RecordArray &source) { this->copy(source); }

RecordArray::RecordArray(RecordArray &&source) { *this = std::move(source); }

RecordArray::~RecordArray() {}

RecordArray::RecordArray(const std::vector<DynamicRecord> &bwt)
    : records(bwt.size()) {
  // Find the starting offsets and compress the BWT.
  std::vector<size_type> offsets(bwt.size());
  for (size_type i = 0; i < bwt.size(); i++) {
    offsets[i] = this->data.size();
    bwt[i].writeBWT(this->data);
  }

  this->buildIndex(offsets);
}

RecordArray::RecordArray(const std::vector<RecordArray const *> sources,
                         const sdsl::int_vector<0> &origins)
    : records(origins.size()) {
  size_type data_size = 0;
  for (auto source : sources) {
    data_size += source->data.size();
  }

  // Initialize the iterators.
  std::vector<sdsl::sd_vector<>::one_iterator> iters;
  iters.reserve(sources.size());
  for (size_type i = 0; i < sources.size(); i++) {
    iters.emplace_back(sources[i]->index.one_begin());
  }

  // Merge the endmarkers.
  {
    DynamicRecord merged;
    for (size_type i = 0; i < sources.size(); i++) {
      if (sources[i]->empty()) {
        continue;
      }
      size_type start = iters[i]->second;
      ++(iters[i]);
      CompressedRecord record(sources[i]->data, start, iters[i]->second);
      for (CompressedRecordIterator iter(record); !(iter.end()); ++iter) {
        run_type run = *iter;
        run.first += merged.outdegree();
        merged.body.push_back(run);
        merged.body_size += run.second;
      }
      for (edge_type outedge : record.outgoing) {
        merged.outgoing.push_back(outedge);
      }
    }
    merged.recode();
    merged.writeBWT(this->data);
  }

  // Merge the BWTs.
  this->data.reserve(data_size + this->data.size());
  std::vector<size_type> offsets(origins.size(), 0);
  std::unordered_set<size_type> active; // These sources cover the current comp.
  for (comp_type comp = 1; comp < origins.size(); comp++) {
    offsets[comp] = this->data.size();
    size_type origin = origins[comp];
    if (origin < sources.size()) {
      active.insert(origin);
    } else {
      this->data.push_back(0);
    } // Empty record.
    std::vector<size_type> to_remove;
    // Advance active iterators and copy the record from the origin.
    for (size_type source : active) {
      if (source == origin) {
        size_type start = iters[source]->second;
        ++(iters[source]);
        size_type limit = iters[source]->second;
        for (size_type i = start; i < limit; i++) {
          this->data.push_back(sources[source]->data[i]);
        }
      } else {
        ++(iters[source]);
      }
      if (iters[source] == sources[source]->index.one_end()) {
        to_remove.push_back(source);
      }
    }
    for (size_type source : to_remove) {
      active.erase(source);
    }
  }

  // Build the index for the BWT.
  this->buildIndex(offsets);
}

RecordArray::RecordArray(size_type array_size) : records(array_size) {}

void RecordArray::buildIndex(const std::vector<size_type> &offsets) {
  sdsl::sd_vector_builder builder(this->data.size(), offsets.size());
  for (size_type offset : offsets) {
    builder.set_unsafe(offset);
  }
  this->index = sdsl::sd_vector<>(builder);
  sdsl::util::init_support(this->select, &(this->index));
}

void RecordArray::swap(RecordArray &another) {
  if (this != &another) {
    std::swap(this->records, another.records), this->index.swap(another.index);
    sdsl::util::swap_support(this->select, another.select, &(this->index),
                             &(another.index));
    this->data.swap(another.data);
  }
}

RecordArray &RecordArray::operator=(const RecordArray &source) {
  if (this != &source) {
    this->copy(source);
  }
  return *this;
}

RecordArray &RecordArray::operator=(RecordArray &&source) {
  if (this != &source) {
    this->records = std::move(source.records);
    this->index = std::move(source.index);
    this->select = std::move(source.select);
    this->select.set_vector(&(this->index));
    this->data = std::move(source.data);
  }
  return *this;
}

size_type RecordArray::serialize(std::ostream &out,
                                 sdsl::structure_tree_node *v,
                                 std::string name) const {
  sdsl::structure_tree_node *child =
      sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += sdsl::write_member(this->records, out, child, "records");
  written_bytes += this->index.serialize(out, child, "index");
  written_bytes += this->select.serialize(out, child, "select");

  // Serialize the data.
  size_type data_bytes = this->data.size() * sizeof(byte_type);
  sdsl::structure_tree_node *data_node = sdsl::structure_tree::add_child(
      child, "data", "std::vector<gbwt::byte_type>");
  if (this->data.size() > 0) {
    DiskIO::write(out, this->data.data(), this->data.size());
  }
  sdsl::structure_tree::add_size(data_node, data_bytes);
  written_bytes += data_bytes;

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void RecordArray::load(std::istream &in) {
  sdsl::read_member(this->records, in);

  // Read the record index.
  this->index.load(in);
  this->select.load(in, &(this->index));

  // Read the data.
  this->data.resize(this->index.size());
  if (this->data.size() > 0) {
    DiskIO::read(in, this->data.data(), this->data.size());
  }

  this->sanityChecks();
}

void RecordArray::simple_sds_serialize(std::ostream &out) const {
  this->index.simple_sds_serialize(out);
  sdsl::simple_sds::serialize_vector(this->data, out);
}

void RecordArray::simple_sds_load(std::istream &in) {
  this->index.simple_sds_load(in);
  this->records = this->index.ones();
  sdsl::util::init_support(this->select, &(this->index));
  this->data = sdsl::simple_sds::load_vector<byte_type>(in);
  this->sanityChecks();
}

size_t RecordArray::simple_sds_size() const {
  return this->index.simple_sds_size() +
         sdsl::simple_sds::vector_size(this->data);
}

size_type RecordArray::size(size_type record) const {
  std::pair<size_type, size_type> range = this->getRange(record);
  return CompressedRecord::recordSize(this->data, range.first, range.second);
}

std::pair<size_type, size_type> RecordArray::getRange(size_type record) const {
  auto iter = this->index.select_iter(record + 1);
  size_type start = iter->second;
  ++iter;
  return std::pair<size_type, size_type>(start, iter->second);
}

void RecordArray::forEach(
    std::function<void(size_type, const CompressedRecord &)> iteratee) const {
  if (this->empty()) {
    return;
  }

  auto iter = this->index.one_begin();
  while (iter != this->index.one_end()) {
    size_type start = iter->second;
    ++iter;
    size_type limit = iter->second;
    CompressedRecord record(this->data, start, limit);
    iteratee(iter->first - 1, record);
  }
}

void RecordArray::copy(const RecordArray &source) {
  this->records = source.records;
  this->index = source.index;
  this->select = source.select;
  this->select.set_vector(&(this->index));
  this->data = source.data;
}

void RecordArray::sanityChecks() const {
  if (this->index.size() != this->data.size()) {
    throw sdsl::simple_sds::InvalidData(
        "RecordArray: Index / data size mismatch");
  }
}

//------------------------------------------------------------------------------

DASamples::DASamples() {}

DASamples::DASamples(const DASamples &source) { this->copy(source); }

DASamples::DASamples(DASamples &&source) { *this = std::move(source); }

DASamples::~DASamples() {}

DASamples::DASamples(const std::vector<DynamicRecord> &bwt) {
  // Determine the statistics and mark the sampled nodes.
  size_type record_count = 0, bwt_offsets = 0, sample_count = 0;
  this->sampled_records = sdsl::bit_vector(bwt.size(), 0);
  for (size_type i = 0; i < bwt.size(); i++) {
    if (bwt[i].samples() > 0) {
      record_count++;
      bwt_offsets += bwt[i].size();
      sample_count += bwt[i].samples();
      this->sampled_records[i] = 1;
    }
  }
  sdsl::util::init_support(this->record_rank, &(this->sampled_records));

  // Build the bitvectors over BWT offsets.
  sdsl::sd_vector_builder range_builder(bwt_offsets, record_count);
  sdsl::sd_vector_builder offset_builder(bwt_offsets, sample_count);
  size_type offset = 0, max_sample = 0;
  for (const DynamicRecord &record : bwt) {
    if (record.samples() > 0) {
      range_builder.set_unsafe(offset);
      for (sample_type sample : record.ids) {
        offset_builder.set_unsafe(offset + sample.first);
        max_sample =
            std::max(max_sample, static_cast<size_type>(sample.second));
      }
      offset += record.size();
    }
  }
  this->bwt_ranges = sdsl::sd_vector<>(range_builder);
  sdsl::util::init_support(this->bwt_select, &(this->bwt_ranges));
  this->sampled_offsets = sdsl::sd_vector<>(offset_builder);

  // Store the samples.
  this->array =
      sdsl::int_vector<0>(sample_count, 0, sdsl::bits::length(max_sample));
  size_type curr = 0;
  for (const DynamicRecord &record : bwt) {
    if (record.samples() > 0) {
      for (sample_type sample : record.ids) {
        this->array[curr] = sample.second;
        curr++;
      }
    }
  }
}

DASamples::DASamples(
    const RecordArray &bwt,
    const std::vector<std::pair<comp_type, sample_type>> &samples) {
  // Determine the statistics and mark the sampled nodes.
  size_type record_count = 0, bwt_offsets = 0;
  this->sampled_records = sdsl::bit_vector(bwt.size(), 0);
  comp_type prev = invalid_comp();
  for (auto sample : samples) {
    if (sample.first != prev) {
      record_count++;
      bwt_offsets += bwt.size(sample.first);
      this->sampled_records[sample.first] = 1;
      prev = sample.first;
    }
  }
  sdsl::util::init_support(this->record_rank, &(this->sampled_records));

  // Build the bitvectors over BWT offsets.
  sdsl::sd_vector_builder range_builder(bwt_offsets, record_count);
  sdsl::sd_vector_builder offset_builder(bwt_offsets, samples.size());
  size_type offset = 0, max_sample = 0;
  prev = invalid_comp();
  for (auto sample : samples) {
    if (sample.first != prev) {
      if (prev != invalid_comp()) {
        offset += bwt.size(prev);
      }
      range_builder.set_unsafe(offset);
      prev = sample.first;
    }
    offset_builder.set_unsafe(offset + sample.second.first);
    max_sample =
        std::max(max_sample, static_cast<size_type>(sample.second.second));
  }
  this->bwt_ranges = sdsl::sd_vector<>(range_builder);
  sdsl::util::init_support(this->bwt_select, &(this->bwt_ranges));
  this->sampled_offsets = sdsl::sd_vector<>(offset_builder);

  // Store the samples.
  this->array =
      sdsl::int_vector<0>(samples.size(), 0, sdsl::bits::length(max_sample));
  for (size_type i = 0; i < samples.size(); i++) {
    this->array[i] = samples[i].second.second;
  }
}

DASamples::DASamples(const std::vector<DASamples const *> sources,
                     const sdsl::int_vector<0> &origins,
                     const std::vector<size_type> &record_offsets,
                     const std::vector<size_type> &sequence_counts) {
  // Compute statistics and build iterators over the sources.
  size_type sample_count = 0, total_sequences = 0;
  std::vector<size_type> sequence_offsets(sources.size(), 0);
  std::vector<SampleIterator> sample_iterators;
  std::vector<SampleRangeIterator> range_iterators;
  for (size_type i = 0; i < sources.size(); i++) {
    sample_count += sources[i]->size();
    sequence_offsets[i] = total_sequences;
    total_sequences += sequence_counts[i];
    sample_iterators.push_back(SampleIterator(*(sources[i])));
    range_iterators.push_back(SampleRangeIterator(*(sources[i])));
  }

  // Compute statistics over the records and mark the sampled nodes.
  // Note that the endmarker requires special treatment.
  size_type record_count = 0, bwt_offsets = 0;
  this->sampled_records = sdsl::bit_vector(origins.size(), 0);
  bool sample_endmarker = false;
  for (size_type origin = 0; origin < sources.size(); origin++) {
    if (sources[origin]->isSampled(ENDMARKER)) {
      sample_endmarker = true;
      ++range_iterators[origin];
    }
  }
  if (sample_endmarker) {
    record_count++;
    bwt_offsets += total_sequences;
    this->sampled_records[ENDMARKER] = 1;
  }
  for (size_type i = 1; i < origins.size(); i++) {
    size_type origin = origins[i];
    if (origin >= sources.size()) {
      continue;
    } // No record.
    if (sources[origin]->isSampled(i - record_offsets[origin])) {
      record_count++;
      bwt_offsets += range_iterators[origin].length();
      this->sampled_records[i] = 1;
      ++range_iterators[origin];
    }
  }
  sdsl::util::init_support(this->record_rank, &(this->sampled_records));

  // Reset the range iterators.
  range_iterators.clear();
  for (size_type i = 0; i < sources.size(); i++) {
    range_iterators.push_back(SampleRangeIterator(*(sources[i])));
  }

  // Build the bitvectors over BWT offsets and store the samples.
  // The endmarker requires special treatment again.
  sdsl::sd_vector_builder range_builder(bwt_offsets, record_count);
  sdsl::sd_vector_builder offset_builder(bwt_offsets, sample_count);
  this->array = sdsl::int_vector<0>(sample_count, 0,
                                    sdsl::bits::length(total_sequences - 1));
  size_type record_start = 0, curr = 0;
  if (sample_endmarker) {
    range_builder.set_unsafe(record_start);
    for (size_type origin = 0; origin < sources.size(); origin++) {
      if (!(sources[origin]->isSampled(ENDMARKER))) {
        continue;
      }
      while (!(sample_iterators[origin].end()) &&
             sample_iterators[origin].offset() <
                 range_iterators[origin].limit()) {
        offset_builder.set_unsafe((sample_iterators[origin]).offset() +
                                  sequence_offsets[origin]);
        this->array[curr] =
            *(sample_iterators[origin]) + sequence_offsets[origin];
        curr++;
        ++sample_iterators[origin];
      }
      ++range_iterators[origin];
    }
    record_start += total_sequences;
  }
  for (size_type i = 1; i < origins.size(); i++) {
    if (!(this->isSampled(i))) {
      continue;
    }
    size_type origin = origins[i];
    range_builder.set_unsafe(record_start);
    while (!(sample_iterators[origin].end()) &&
           sample_iterators[origin].offset() <
               range_iterators[origin].limit()) {
      offset_builder.set_unsafe((sample_iterators[origin].offset() -
                                 range_iterators[origin].start()) +
                                record_start);
      this->array[curr] =
          *(sample_iterators[origin]) + sequence_offsets[origin];
      curr++;
      ++sample_iterators[origin];
    }
    record_start += range_iterators[origin].length();
    ++range_iterators[origin];
  }
  this->bwt_ranges = sdsl::sd_vector<>(range_builder);
  sdsl::util::init_support(this->bwt_select, &(this->bwt_ranges));
  this->sampled_offsets = sdsl::sd_vector<>(offset_builder);
}

void DASamples::swap(DASamples &another) {
  if (this != &another) {
    this->sampled_records.swap(another.sampled_records);
    sdsl::util::swap_support(this->record_rank, another.record_rank,
                             &(this->sampled_records),
                             &(another.sampled_records));

    this->bwt_ranges.swap(another.bwt_ranges);
    sdsl::util::swap_support(this->bwt_select, another.bwt_select,
                             &(this->bwt_ranges), &(another.bwt_ranges));

    this->sampled_offsets.swap(another.sampled_offsets);

    this->array.swap(another.array);
  }
}

DASamples &DASamples::operator=(const DASamples &source) {
  if (this != &source) {
    this->copy(source);
  }
  return *this;
}

DASamples &DASamples::operator=(DASamples &&source) {
  if (this != &source) {
    this->sampled_records = std::move(source.sampled_records);
    this->record_rank = std::move(source.record_rank);

    this->bwt_ranges = std::move(source.bwt_ranges);
    this->bwt_select = std::move(source.bwt_select);

    this->sampled_offsets = std::move(source.sampled_offsets);

    this->array = std::move(source.array);

    this->setVectors();
  }
  return *this;
}

size_type DASamples::serialize(std::ostream &out, sdsl::structure_tree_node *v,
                               std::string name) const {
  sdsl::structure_tree_node *child =
      sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes +=
      this->sampled_records.serialize(out, child, "sampled_records");
  written_bytes += this->record_rank.serialize(out, child, "record_rank");

  written_bytes += this->bwt_ranges.serialize(out, child, "bwt_ranges");
  written_bytes += this->bwt_select.serialize(out, child, "bwt_select");

  written_bytes +=
      this->sampled_offsets.serialize(out, child, "sampled_offsets");

  written_bytes += this->array.serialize(out, child, "array");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void DASamples::load(std::istream &in) {
  this->sampled_records.load(in);
  this->record_rank.load(in, &(this->sampled_records));

  this->bwt_ranges.load(in);
  this->bwt_select.load(in, &(this->bwt_ranges));

  this->sampled_offsets.load(in);

  this->array.load(in);

  this->sanityChecks();
}

void DASamples::simple_sds_serialize(std::ostream &out) const {
  this->sampled_records.simple_sds_serialize(out);
  this->bwt_ranges.simple_sds_serialize(out);
  this->sampled_offsets.simple_sds_serialize(out);
  this->array.simple_sds_serialize(out);
}

void DASamples::simple_sds_load(std::istream &in) {
  this->sampled_records.simple_sds_load(in);
  sdsl::util::init_support(this->record_rank, &(this->sampled_records));

  this->bwt_ranges.simple_sds_load(in);
  sdsl::util::init_support(this->bwt_select, &(this->bwt_ranges));

  this->sampled_offsets.simple_sds_load(in);
  this->array.simple_sds_load(in);

  this->sanityChecks();
}

size_t DASamples::simple_sds_size() const {
  return this->sampled_records.simple_sds_size() +
         this->bwt_ranges.simple_sds_size() +
         this->sampled_offsets.simple_sds_size() +
         this->array.simple_sds_size();
}

void DASamples::copy(const DASamples &source) {
  this->sampled_records = source.sampled_records;
  this->record_rank = source.record_rank;

  this->bwt_ranges = source.bwt_ranges;
  this->bwt_select = source.bwt_select;

  this->sampled_offsets = source.sampled_offsets;

  this->array = source.array;

  this->setVectors();
}

void DASamples::setVectors() {
  this->record_rank.set_vector(&(this->sampled_records));
  this->bwt_select.set_vector(&(this->bwt_ranges));
}

void DASamples::sanityChecks() const {
  if (this->record_rank(this->sampled_records.size()) !=
      this->bwt_ranges.ones()) {
    throw sdsl::simple_sds::InvalidData(
        "DASamples: Sampled record / BWT range count mismatch");
  }
  if (this->bwt_ranges.size() != this->sampled_offsets.size()) {
    throw sdsl::simple_sds::InvalidData(
        "DASamples: BWT range / sampled offsets size mismatch");
  }
  if (this->sampled_offsets.ones() != this->array.size()) {
    throw sdsl::simple_sds::InvalidData(
        "DASamples: Sampled offset / sample count mismatch");
  }
}

size_type DASamples::tryLocate(size_type record, size_type offset) const {
  if (!(this->isSampled(record))) {
    return invalid_sequence();
  }

  size_type record_start = this->start(record);
  auto iter = this->sampled_offsets.predecessor(record_start + offset);
  if (iter->second == record_start + offset) {
    return this->array[iter->first];
  }
  return invalid_sequence();
}

sample_type DASamples::nextSample(size_type record, size_type offset) const {
  if (!(this->isSampled(record))) {
    return invalid_sample();
  }

  size_type record_start = this->start(record);
  auto iter = this->sampled_offsets.successor(record_start + offset);
  if (iter->first < this->array.size()) {
    return sample_type(iter->second - record_start, this->array[iter->first]);
  }
  return invalid_sample();
}

//------------------------------------------------------------------------------

MergeParameters::MergeParameters()
    : pos_buffer_size(POS_BUFFER_SIZE), thread_buffer_size(THREAD_BUFFER_SIZE),
      merge_buffers(MERGE_BUFFERS), chunk_size(CHUNK_SIZE),
      merge_jobs(MERGE_JOBS) {}

void MergeParameters::setPosBufferSize(size_type megabytes) {
  this->pos_buffer_size = Range::bound(megabytes, 1, MAX_BUFFER_SIZE);
}

void MergeParameters::setThreadBufferSize(size_type megabytes) {
  this->thread_buffer_size = Range::bound(megabytes, 1, MAX_BUFFER_SIZE);
}

void MergeParameters::setMergeBuffers(size_type n) {
  this->merge_buffers = Range::bound(n, 1, MAX_MERGE_BUFFERS);
}

void MergeParameters::setChunkSize(size_type n) {
  this->chunk_size = std::max(n, static_cast<size_type>(1));
}

void MergeParameters::setMergeJobs(size_type n) {
  this->merge_jobs = Range::bound(n, 1, MAX_MERGE_JOBS);
}

//------------------------------------------------------------------------------

StringArray::StringArray(const std::vector<std::string> &source) {
  *this = StringArray(
      source.size(), [](size_type) -> bool { return true; },
      [&source](size_type i) -> size_type { return source[i].length(); },
      [&source](size_type i) -> view_type { return str_to_view(source[i]); });
}

StringArray::StringArray(const std::map<std::string, std::string> &source) {
  std::vector<std::string> linearized;
  for (auto iter = source.begin(); iter != source.end(); ++iter) {
    linearized.push_back(iter->first);
    linearized.push_back(iter->second);
  }
  *this = StringArray(linearized);
}

StringArray::StringArray(size_type n,
                         const std::function<size_type(size_type)> &length,
                         const std::function<view_type(size_type)> &sequence) {
  *this =
      StringArray(n, [](size_type) -> bool { return true; }, length, sequence);
}

StringArray::StringArray(size_type n,
                         const std::function<bool(size_type)> &choose,
                         const std::function<size_type(size_type)> &length,
                         const std::function<view_type(size_type)> &sequence) {
  size_type chosen = 0, total_length = 0;
  for (size_type i = 0; i < n; i++) {
    if (choose(i)) {
      chosen++;
      total_length += length(i);
    }
  }
  this->index =
      sdsl::int_vector<0>(chosen + 1, 0, sdsl::bits::length(total_length));
  this->strings.reserve(total_length);

  size_type curr = 0, total = 0;
  for (size_type i = 0; i < n; i++) {
    if (!choose(i)) {
      continue;
    }
    view_type view = sequence(i);
    this->index[curr] = total;
    curr++;
    this->strings.insert(this->strings.end(), view.first,
                         view.first + view.second);
    total += view.second;
  }
  this->index[chosen] = total;
}

// This has a separate implementation, because we cannot take a view of a
// temporary string.
StringArray::StringArray(
    size_type n, const std::function<size_type(size_type)> &length,
    const std::function<std::string(size_type)> &sequence) {
  size_type total_length = 0;
  for (size_type i = 0; i < n; i++) {
    total_length += length(i);
  }
  this->index = sdsl::int_vector<0>(n + 1, 0, sdsl::bits::length(total_length));
  this->strings.reserve(total_length);

  size_type total = 0;
  for (size_type i = 0; i < n; i++) {
    std::string str = sequence(i);
    this->index[i] = total;
    this->strings.insert(this->strings.end(), str.begin(), str.end());
    total += str.length();
  }
  this->index[n] = total;
}

void StringArray::swap(StringArray &another) {
  this->index.swap(another.index);
  this->strings.swap(another.strings);
}

size_type StringArray::serialize(std::ostream &out,
                                 sdsl::structure_tree_node *v,
                                 std::string name) const {
  sdsl::structure_tree_node *child =
      sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  // SDSL serialization format stores the strings before the index in order to
  // be compatible with the old implementation in GBWTGraph.
  written_bytes += serializeVector(this->strings, out, child, "strings");
  this->index.serialize(out, child, "index");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void StringArray::load(std::istream &in) {
  loadVector(this->strings, in);
  this->index.load(in);
  this->sanityChecks();
}

void determine_alphabet(const std::vector<char> &strings,
                        std::vector<std::uint8_t> &char_to_comp,
                        sdsl::int_vector<8> &comp_to_char, size_type &width) {
  char_to_comp = std::vector<std::uint8_t>(256, 0);
  for (char c : strings) {
    char_to_comp[static_cast<std::uint8_t>(c)] = 1;
  }

  size_type sigma = 0;
  for (auto c : char_to_comp) {
    sigma += c;
  }
  width = sdsl::bits::length(std::max(sigma, size_type(1)) - 1);

  comp_to_char = sdsl::int_vector<8>(sigma, 0);
  for (size_type i = 0, found = 0; i < char_to_comp.size(); i++) {
    if (char_to_comp[i] != 0) {
      comp_to_char[found] = i;
      char_to_comp[i] = found;
      found++;
    }
  }
}

void StringArray::simple_sds_serialize(std::ostream &out) const {
  // Compress the index without the past-the-end sentinel.
  {
    sdsl::sd_vector<> v(this->index.begin(), this->index.end() - 1);
    v.simple_sds_serialize(out);
  }

  // Determine and serialize the alphabet.
  std::vector<std::uint8_t> char_to_comp;
  sdsl::int_vector<8> comp_to_char;
  size_type width = 0;
  determine_alphabet(this->strings, char_to_comp, comp_to_char, width);
  comp_to_char.simple_sds_serialize(out);

  // Compress the strings.
  {
    sdsl::int_vector<> compressed(this->strings.size(), 0, width);
    for (size_type i = 0; i < this->strings.size(); i++) {
      compressed[i] = char_to_comp[static_cast<uint8_t>(this->strings[i])];
    }
    compressed.simple_sds_serialize(out);
  }
}

void StringArray::simple_sds_load(std::istream &in) {
  // Load the index. We cannot compress it yet because we do not know the width
  // of the sentinel value `strings.size()`.
  sdsl::sd_vector<> v;
  v.simple_sds_load(in);

  // Load the alphabet.
  sdsl::int_vector<8> comp_to_char;
  comp_to_char.simple_sds_load(in);

  // Decompress the strings.
  {
    sdsl::int_vector<> compressed;
    compressed.simple_sds_load(in);
    this->strings = std::vector<char>();
    this->strings.reserve(compressed.size());
    for (auto c : compressed) {
      this->strings.push_back(comp_to_char[c]);
    }
  }

  // Decompress the index.
  this->index = sdsl::int_vector<>(v.ones() + 1, 0,
                                   sdsl::bits::length(this->strings.size()));
  size_t i = 0;
  for (auto iter = v.one_begin(); iter != v.one_end(); ++iter, i++) {
    this->index[i] = iter->second;
  }
  this->index[this->size()] = this->strings.size();

  this->sanityChecks();
}

size_t StringArray::simple_sds_size() const {
  size_t result = 0;

  // Compress the index without the past-the-end sentinel.
  size_t universe = (this->empty() ? 0 : this->index[this->size() - 1] + 1);
  result += sdsl::sd_vector<>::simple_sds_size(universe, this->size());

  // Determine the alphabet.
  std::vector<std::uint8_t> char_to_comp;
  sdsl::int_vector<8> comp_to_char;
  size_type width = 0;
  determine_alphabet(this->strings, char_to_comp, comp_to_char, width);
  result += comp_to_char.simple_sds_size();

  // Compress the strings.
  result += sdsl::int_vector<>::simple_sds_size(this->strings.size(), width);

  return result;
}

void StringArray::remove(size_type i) {
  if (i >= this->size()) {
    return;
  }

  // Update strings.
  size_type tail = this->index[i];
  size_type diff = this->index[i + 1] - tail;
  while (tail + diff < this->strings.size()) {
    this->strings[tail] = this->strings[tail + diff];
    tail++;
  }
  this->strings.resize(tail);

  // Update index.
  for (size_type j = i; j + 1 < this->index.size(); j++) {
    this->index[j] = this->index[j + 1] - diff;
  }
  this->index.resize(this->index.size() - 1);
  sdsl::util::bit_compress(this->index);
}

bool StringArray::operator==(const StringArray &another) const {
  return (this->index == another.index && this->strings == another.strings);
}

bool StringArray::operator!=(const StringArray &another) const {
  return !(this->operator==(another));
}

void StringArray::sanityChecks() const {
  if (this->index.size() == 0 || this->index[0] != 0 ||
      this->index[this->index.size() - 1] != this->strings.size()) {
    throw sdsl::simple_sds::InvalidData(
        "StringArray: Offsets and strings do not match");
  }
}

//------------------------------------------------------------------------------

Dictionary::Dictionary() {}

Dictionary::Dictionary(const Dictionary &source) { this->copy(source); }

Dictionary::Dictionary(Dictionary &&source) { *this = std::move(source); }

Dictionary::~Dictionary() {}

Dictionary::Dictionary(const std::vector<std::string> &source)
    : strings(source) {
  // Build sorted_ids and check for duplicates.
  this->sortKeys();
  if (this->hasDuplicates() && Verbosity::level >= Verbosity::FULL) {
    std::cerr << "Dictionary::Dictionary(): Warning: The dictionary contains "
                 "duplicate strings"
              << std::endl;
  }
}

Dictionary::Dictionary(const Dictionary &first, const Dictionary &second) {
  if (first.empty()) {
    *this = second;
    return;
  } else if (second.empty()) {
    *this = first;
    return;
  }

  // Add all strings from the first and missing strings from the second.
  this->strings = StringArray(
      first.size() + second.size(),
      [&](size_type i) -> bool {
        return (i < first.size() || first.find(second.strings.view(
                                        i - first.size())) >= first.size());
      },
      [&](size_type i) -> size_type {
        return (i < first.size() ? first.strings.length(i)
                                 : second.strings.length(i - first.size()));
      },
      [&](size_type i) -> view_type {
        return (i < first.size() ? first.strings.view(i)
                                 : second.strings.view(i - first.size()));
      });

  // Build sorted_ids and check for duplicates.
  this->sortKeys();
  if (this->hasDuplicates() && Verbosity::level >= Verbosity::FULL) {
    std::cerr << "Dictionary::Dictionary(): Warning: The dictionary contains "
                 "duplicate strings"
              << std::endl;
  }
}

void Dictionary::swap(Dictionary &another) {
  if (this != &another) {
    this->strings.swap(another.strings);
    this->sorted_ids.swap(another.sorted_ids);
  }
}

Dictionary &Dictionary::operator=(const Dictionary &source) {
  if (this != &source) {
    this->copy(source);
  }
  return *this;
}

Dictionary &Dictionary::operator=(Dictionary &&source) {
  if (this != &source) {
    this->strings = std::move(source.strings);
    this->sorted_ids = std::move(source.sorted_ids);
  }
  return *this;
}

size_type Dictionary::serialize(std::ostream &out, sdsl::structure_tree_node *v,
                                std::string name) const {
  sdsl::structure_tree_node *child =
      sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->strings.serialize(out, child, "strings");
  written_bytes += this->sorted_ids.serialize(out, child, "sorted_ids");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void Dictionary::load(std::istream &in) {
  this->strings.load(in);
  this->sorted_ids.load(in);
  this->sanityChecks();
}

void Dictionary::load_v1(std::istream &in) {
  sdsl::int_vector<0> offsets;
  offsets.load(in);
  this->sorted_ids.load(in);
  std::vector<char> data;
  loadVector(data, in);
  this->strings = StringArray(
      offsets.size() - 1,
      [&](size_type i) -> size_type { return offsets[i + 1] - offsets[i]; },
      [&](size_type i) -> view_type {
        return view_type(data.data() + offsets[i], offsets[i + 1] - offsets[i]);
      });
  this->sanityChecks();
}

void Dictionary::simple_sds_serialize(std::ostream &out) const {
  this->strings.simple_sds_serialize(out);
  this->sorted_ids.simple_sds_serialize(out);
}

void Dictionary::simple_sds_load(std::istream &in) {
  this->strings.simple_sds_load(in);
  this->sorted_ids.simple_sds_load(in);
  this->sanityChecks();
}

size_t Dictionary::simple_sds_size() const {
  return this->strings.simple_sds_size() + this->sorted_ids.simple_sds_size();
}

void Dictionary::copy(const Dictionary &source) {
  this->strings = source.strings;
  this->sorted_ids = source.sorted_ids;
}

void Dictionary::sanityChecks() const {
  if (this->sorted_ids.size() != this->strings.size()) {
    throw sdsl::simple_sds::InvalidData(
        "Dictionary: Size mismatch between strings and sorted ids");
  }
}

void Dictionary::sortKeys() {
  if (this->empty()) {
    this->sorted_ids = sdsl::int_vector<0>();
    return;
  }

  this->sorted_ids = sdsl::int_vector<0>(this->size(), 0,
                                         sdsl::bits::length(this->size() - 1));
  for (size_type i = 0; i < this->size(); i++) {
    this->sorted_ids[i] = i;
  }
  sequentialSort(this->sorted_ids.begin(), this->sorted_ids.end(),
                 [this](size_type a, size_type b) -> bool {
                   return this->smaller_by_id(a, b);
                 });
}

bool Dictionary::operator==(const Dictionary &another) const {
  return (this->strings == another.strings &&
          this->sorted_ids == another.sorted_ids);
}

void Dictionary::clear() { *this = Dictionary(); }

size_type Dictionary::find(view_type view) const {
  size_type low = 0, high = this->size();
  while (low < high) {
    size_type mid = low + (high - low) / 2;
    if (this->smaller_by_order(view, mid)) {
      high = mid;
    } else if (this->smaller_by_order(mid, view)) {
      low = mid + 1;
    } else {
      return this->sorted_ids[mid];
    }
  }
  return this->size();
}

void Dictionary::remove(size_type i) {
  if (i >= this->size()) {
    return;
  }

  // Update strings.
  this->strings.remove(i);

  // Update sorted_ids.
  size_type diff = 0;
  for (size_type j = 0; j + 1 < this->sorted_ids.size(); j++) {
    if (this->sorted_ids[j] == i) {
      diff = 1;
    }
    this->sorted_ids[j] = this->sorted_ids[j + diff];
    if (this->sorted_ids[j] > i) {
      this->sorted_ids[j]--;
    }
  }
  this->sorted_ids.resize(this->sorted_ids.size() - 1);
}

void Dictionary::append(const Dictionary &source) {
  *this = Dictionary(*this, source);
}

bool Dictionary::hasDuplicates() const {
  for (size_type i = 0; i + 1 < this->size(); i++) {
    if (!(this->smaller_by_order(i, i + 1))) {
      return true;
    }
  }
  return false;
}

bool stringCompare(const char *a_pos, const char *a_lim, const char *b_pos,
                   const char *b_lim) {
  while (a_pos != a_lim && b_pos != b_lim) {
    if (*a_pos != *b_pos) {
      return (*a_pos < *b_pos);
    }
    ++a_pos;
    ++b_pos;
  }
  return (a_pos == a_lim && b_pos != b_lim);
}

bool Dictionary::smaller_by_order(size_type a, size_type b) const {
  view_type first = this->strings.view(this->sorted_ids[a]);
  view_type second = this->strings.view(this->sorted_ids[b]);
  return stringCompare(first.first, first.first + first.second, second.first,
                       second.first + second.second);
}

bool Dictionary::smaller_by_order(size_type a, view_type b) const {
  view_type first = this->strings.view(this->sorted_ids[a]);
  return stringCompare(first.first, first.first + first.second, b.first,
                       b.first + b.second);
}

bool Dictionary::smaller_by_order(view_type a, size_type b) const {
  view_type second = this->strings.view(this->sorted_ids[b]);
  return stringCompare(a.first, a.first + a.second, second.first,
                       second.first + second.second);
}

bool Dictionary::smaller_by_id(size_type a, size_type b) const {
  view_type first = this->strings.view(a);
  view_type second = this->strings.view(b);
  return stringCompare(first.first, first.first + first.second, second.first,
                       second.first + second.second);
}

bool Dictionary::smaller_by_id(size_type a, view_type b) const {
  view_type first = this->strings.view(a);
  return stringCompare(first.first, first.first + first.second, b.first,
                       b.first + b.second);
}

bool Dictionary::smaller_by_id(view_type a, size_type b) const {
  view_type second = this->strings.view(b);
  return stringCompare(a.first, a.first + a.second, second.first,
                       second.first + second.second);
}

//------------------------------------------------------------------------------

Tags::Tags() {}

Tags::Tags(const Tags &source) { this->copy(source); }

Tags::Tags(Tags &&source) { *this = std::move(source); }

Tags::~Tags() {}

void Tags::swap(Tags &another) {
  if (this != &another) {
    this->tags.swap(another.tags);
  }
}

Tags &Tags::operator=(const Tags &source) {
  if (this != &source) {
    this->copy(source);
  }
  return *this;
}

Tags &Tags::operator=(Tags &&source) {
  if (this != &source) {
    this->tags = std::move(source.tags);
  }
  return *this;
}

size_type Tags::serialize(std::ostream &out, sdsl::structure_tree_node *v,
                          std::string name) const {
  sdsl::structure_tree_node *child =
      sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  {
    StringArray linearized(this->tags);
    written_bytes += linearized.serialize(out, child, "tags");
  }

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void Tags::load(std::istream &in) {
  StringArray linearized;
  linearized.load(in);
  this->build(linearized);
}

void Tags::simple_sds_serialize(std::ostream &out) const {
  StringArray linearized(this->tags);
  linearized.simple_sds_serialize(out);
}

void Tags::simple_sds_load(std::istream &in) {
  StringArray linearized;
  linearized.simple_sds_load(in);
  this->build(linearized);
}

size_t Tags::simple_sds_size() const {
  StringArray linearized(this->tags);
  return linearized.simple_sds_size();
}

void Tags::copy(const Tags &source) { this->tags = source.tags; }

void Tags::build(const StringArray &source) {
  if (source.size() % 2 != 0) {
    throw sdsl::simple_sds::InvalidData("Tags: Key without a value");
  }

  this->tags.clear();
  for (size_type i = 0; i < source.size(); i += 2) {
    std::string key = source.str(i);
    for (auto iter = key.begin(); iter != key.end(); ++iter) {
      *iter = std::tolower(*iter);
    }
    this->tags[key] = source.str(i + 1);
  }

  if (this->tags.size() != source.size() / 2) {
    throw sdsl::simple_sds::InvalidData("Tags: Duplicate keys");
  }
}

bool Tags::operator==(const Tags &another) const {
  return (this->tags == another.tags);
}

void Tags::set(const std::string &key, const std::string &value) {
  this->tags[normalize(key)] = value;
}

std::string Tags::get(const std::string &key) const {
  auto iter = this->tags.find(normalize(key));
  if (iter == this->tags.end()) {
    return std::string();
  }
  return iter->second;
}

bool Tags::contains(const std::string &key) const {
  return (this->tags.find(normalize(key)) != this->tags.end());
}

void Tags::clear() { *this = Tags(); }

std::string Tags::normalize(const std::string &key) {
  std::string normalized = key;
  for (auto iter = normalized.begin(); iter != normalized.end(); ++iter) {
    *iter = std::tolower(*iter);
  }
  return normalized;
}

//------------------------------------------------------------------------------

} // namespace gbwt
