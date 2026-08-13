/*
  Copyright (c) 2017, 2018, 2019, 2020, 2021, 2023 Jouni Siren
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

#include "absl/log/absl_log.h"
#include <gbwt/dynamic_gbwt.h>
#include <gbwt/bwtmerge.h>

#include <memory>
#include <unordered_map>

namespace gbwt
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_type DynamicGBWT::INSERT_BATCH_SIZE;
constexpr size_type DynamicGBWT::MIN_SEQUENCES_PER_BATCH;
constexpr size_type DynamicGBWT::REMOVE_CHUNK_SIZE;
constexpr size_type DynamicGBWT::MERGE_BATCH_SIZE;
constexpr size_type DynamicGBWT::SAMPLE_INTERVAL;

//------------------------------------------------------------------------------

// Other class variables.

const std::string DynamicGBWT::EXTENSION = ".gbwt";

//------------------------------------------------------------------------------

DynamicGBWT::DynamicGBWT()
{
  this->addSource();
}

DynamicGBWT::DynamicGBWT(const DynamicGBWT& source)
{
  this->copy(source);
}

DynamicGBWT::DynamicGBWT(const GBWT& source)
{
  this->copy(source);
}

DynamicGBWT::DynamicGBWT(DynamicGBWT&& source)
{
  *this = std::move(source);
}

DynamicGBWT::~DynamicGBWT()
{
}

void
DynamicGBWT::swap(DynamicGBWT& another)
{
  if(this != &another)
  {
    this->header.swap(another.header);
    this->tags.swap(another.tags);
    this->bwt.swap(another.bwt);
    this->metadata.swap(another.metadata);
  }
}

DynamicGBWT&
DynamicGBWT::operator=(const DynamicGBWT& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

DynamicGBWT&
DynamicGBWT::operator=(const GBWT& source)
{
  this->copy(source);
  return *this;
}

DynamicGBWT&
DynamicGBWT::operator=(DynamicGBWT&& source)
{
  if(this != &source)
  {
    this->header = std::move(source.header);
    this->tags = std::move(source.tags);
    this->bwt = std::move(source.bwt);
    this->metadata = std::move(source.metadata);
  }
  return *this;
}

void
DynamicGBWT::resample(size_type sample_interval)
{
  // Delete old samples to save memory.
  for(DynamicRecord& record : this->bwt) { record.ids = std::vector<sample_type>(); }
  std::vector<std::pair<node_type, sample_type>> samples = gbwt::resample(*this, sample_interval);
  for(auto sample : samples) { this->bwt[sample.first].ids.push_back(sample.second); }
}

size_type
DynamicGBWT::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->header.serialize(out, child, "header");
  written_bytes += this->tags.serialize(out, child, "tags");

  {
    RecordArray array(this->bwt);
    written_bytes += array.serialize(out, child, "bwt");
  }

  {
    DASamples compressed_samples(this->bwt);
    written_bytes += compressed_samples.serialize(out, child, "da_samples");
  }

  if(this->hasMetadata())
  {
    written_bytes += this->metadata.serialize(out, child, "metadata");
  }

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
DynamicGBWT::load(std::istream& in)
{
  // Read the header.
  GBWTHeader h = sdsl::simple_sds::load_value<GBWTHeader>(in);
  h.check();
  bool simple_sds = h.get(GBWTHeader::FLAG_SIMPLE_SDS);
  bool has_tags = h.version >= 5; // FIXME Replace with symbolic constant.
  h.unset(GBWTHeader::FLAG_SIMPLE_SDS); // We only set this flag in the serialized header.
  h.setVersion(); // Update to the current version.
  this->header = h;

  // Read the tags and set the source to Version::SOURCE_VALUE.
  bool source_is_self = false; // We know how to read DASamples.
  if(has_tags)
  {
    if(simple_sds) { this->tags.simple_sds_load(in); }
    else { this->tags.load(in); }
    if(this->tags.get(Version::SOURCE_KEY) == Version::SOURCE_VALUE) { source_is_self = true; }
    this->addSource();
  }
  else { this->resetTags(); }

  // Read and decompress the BWT.
  this->bwt.resize(this->effective());
  {
    RecordArray array;
    if(simple_sds) { array.simple_sds_load(in); }
    else { array.load(in); }
    if(array.size() != this->effective())
    {
      ABSL_LOG(FATAL) << "DynamicGBWT: BWT record count / alphabet size mismatch";
    }
    array.forEach([&](size_type comp, const CompressedRecord& record)
    {
      DynamicRecord& current = this->bwt[comp];
      current.clear();
      current.outgoing = record.outgoing;
      if(current.outdegree() > 0)
      {
        for(CompressedRecordIterator iter(record); !(iter.end()); ++iter)
        {
          current.body.push_back(*iter);
          current.body_size += iter->second;
        }
      }
    });
  }

  // Read and decompress the samples.
  {
    DASamples samples;
    bool found_samples = false;
    if(simple_sds)
    {
      // The samples may be absent or from an unknown source.
      if(source_is_self) { found_samples = sdsl::simple_sds::load_option(samples, in); }
      else { sdsl::simple_sds::skip_option(in); }
    }
    else { samples.load(in); found_samples = true; }
    if(found_samples)
    {
      if(samples.records() != this->effective())
      {
        ABSL_LOG(FATAL) << "DynamicGBWT: Sample record count / alphabet size mismatch";
      }
      SampleIterator sample_iter(samples);
      for(SampleRangeIterator range_iter(samples); !(range_iter.end()); ++range_iter)
      {
        DynamicRecord& current = this->bwt[range_iter.record()];
        while(!(sample_iter.end()) && sample_iter.offset() < range_iter.limit())
        {
          current.ids.push_back(sample_type(sample_iter.offset() - range_iter.start(), *sample_iter));
          ++sample_iter;
        }
      }
    }
    else { this->resample(SAMPLE_INTERVAL); }
  }

  // Read the metadata.
  if(simple_sds)
  {
    bool loaded_metadata = sdsl::simple_sds::load_option(this->metadata, in);
    if(loaded_metadata != this->hasMetadata())
    {
      ABSL_LOG(FATAL) << "DynamicGBWT: Invalid metadata flag in the header";
    }
  }
  else if(this->hasMetadata()) { this->metadata.load(in); }
  if(this->hasMetadata() && this->metadata.hasPathNames())
  {
    size_type expected_paths = (this->bidirectional() ? this->sequences() / 2 : this->sequences());
    if(this->metadata.paths() != expected_paths)
    {
      ABSL_LOG(FATAL) << "DynamicGBWT: Path name / sequence count mismatch";
    }
  }

  // Rebuild the incoming edges.
  this->rebuildIncoming();
}

void
DynamicGBWT::simple_sds_serialize(std::ostream& out) const
{
  GBWTHeader h = this->header;
  h.set(GBWTHeader::FLAG_SIMPLE_SDS); // We only set this flag in the serialized header.
  sdsl::simple_sds::serialize_value(h, out);

  this->tags.simple_sds_serialize(out);
  {
    RecordArray array(this->bwt);
    array.simple_sds_serialize(out);
  }
  {
    DASamples compressed_samples(this->bwt);
    sdsl::simple_sds::serialize_option(compressed_samples, out);
  }
  if(this->hasMetadata()) { sdsl::simple_sds::serialize_option(this->metadata, out); }
  else { sdsl::simple_sds::empty_option(out); }
}

void
DynamicGBWT::simple_sds_load(std::istream& in)
{
  // The same load() function can handle the SDSL and simple-sds formats.
  this->load(in);
}

size_t
DynamicGBWT::simple_sds_size() const
{
  size_t result = sdsl::simple_sds::value_size(this->header);
  result += this->tags.simple_sds_size();
  {
    RecordArray array(this->bwt);
    result += array.simple_sds_size();
  }
  {
    DASamples compressed_samples(this->bwt);
    result += sdsl::simple_sds::option_size(compressed_samples);
  }
  if(this->hasMetadata()) { result += sdsl::simple_sds::option_size(this->metadata); }
  else { result += sdsl::simple_sds::empty_option_size(); }
  return result;
}

void
DynamicGBWT::copy(const DynamicGBWT& source)
{
  this->header = source.header;
  this->tags = source.tags;
  this->bwt = source.bwt;
  this->metadata = source.metadata;
}

void
DynamicGBWT::copy(const GBWT& source)
{
  this->header = source.header;
  this->tags = source.tags;

  // Decompress the BWT.
  this->bwt.resize(this->effective());
  source.bwt.forEach([&](size_type comp, const CompressedRecord& record)
  {
    DynamicRecord& current = this->bwt[comp];
    current.clear();
    current.outgoing = record.outgoing;
    if(current.outdegree() > 0)
    {
      for(CompressedRecordIterator iter(record); !(iter.end()); ++iter)
      {
        current.body.push_back(*iter);
        current.body_size += iter->second;
      }
    }
  });

  // Decompress the samples.
  SampleIterator sample_iter(source.da_samples);
  for(SampleRangeIterator range_iter(source.da_samples); !(range_iter.end()); ++range_iter)
  {
    DynamicRecord& current = this->bwt[range_iter.record()];
    while(!(sample_iter.end()) && sample_iter.offset() < range_iter.limit())
    {
      current.ids.push_back(sample_type(sample_iter.offset() - range_iter.start(), *sample_iter));
      ++sample_iter;
    }
  }

  // Rebuild the incoming edges.
  this->rebuildIncoming();

  this->metadata = source.metadata;
}

void
DynamicGBWT::resetTags()
{
  this->tags.clear();
  this->addSource();
}

void
DynamicGBWT::addSource()
{
  this->tags.set(Version::SOURCE_KEY, Version::SOURCE_VALUE);
}

//------------------------------------------------------------------------------

std::pair<size_type, size_type>
DynamicGBWT::runs() const
{
  std::pair<size_type, size_type> result(0, 0);
  for(const DynamicRecord& node : this->bwt)
  {
    std::pair<size_type, size_type> temp = node.runs();
    result.first += temp.first; result.second += temp.second;
  }
  return result;
}

size_type
DynamicGBWT::samples() const
{
  size_type total = 0;
  for(const DynamicRecord& node : this->bwt) { total += node.samples(); }
  return total;
}

edge_type
DynamicGBWT::inverseLF(node_type from, size_type i) const
{
  if(!(this->bidirectional()) || from == ENDMARKER) { return invalid_edge(); }

  // Find the predecessor node id.
  DynamicRecord curr = this->record(from);
  node_type predecessor = invalid_node();
  size_type curr_offset = 0;
  for(edge_type inedge : curr.incoming)
  {
    curr_offset += inedge.second;
    if(curr_offset > i) { predecessor = inedge.first; break; }
  }
  if(predecessor == invalid_node()) { return invalid_edge(); }

  // Determine the offset.
  DynamicRecord pred_record = this->record(predecessor);
  size_type offset = pred_record.offsetTo(from, i);
  if(offset == invalid_offset()) { return invalid_edge(); }

  return edge_type(predecessor, offset);
}

size_type
DynamicGBWT::fullLF(node_type from, size_type i, node_type to) const
{
  if(to == ENDMARKER) { return invalid_offset(); }
  else if(!(this->contains(to))) { return 0; }
  else if(!(this->hasEdge(from, to))) { return this->record(to).countBefore(from); }
  else { return this->record(from).LF(i, to); }
}

//------------------------------------------------------------------------------

void
DynamicGBWT::resize(size_type new_offset, size_type new_sigma)
{
  /*
    Do not set the new offset, if we already have a smaller real offset or the
    new offset is not a real one.
  */
  if((this->sigma() > 1 && new_offset > this->header.offset) || new_sigma <= 1)
  {
    new_offset = this->header.offset;
  }
  if(this->sigma() > new_sigma) { new_sigma = this->sigma(); }
  if(new_offset > 0 && new_offset >= new_sigma)
  {
    std::cerr << "DynamicGBWT::resize(): Cannot set offset " << new_offset << " with alphabet size " << new_sigma << std::endl;
    std::exit(EXIT_FAILURE);
  }

  this->forceResize(new_offset, new_sigma);
}

void
DynamicGBWT::resize()
{
  if(this->effective() == 0) { return; }

  // Determine the new effective alphabet.
  size_type new_offset = this->header.offset, new_sigma = this->sigma();
  if(this->empty())
  {
    if(Verbosity::level >= Verbosity::FULL)
    {
      std::cerr << "DynamicGBWT::resize(): The index is empty" << std::endl;
    }
    new_offset = 0; new_sigma = 0;
  }
  else
  {
    size_type head_empty = 0;
    while(head_empty + 1 < this->effective() && this->bwt[head_empty + 1].empty()) { head_empty++; }
    if(head_empty + 1 >= this->effective())
    {
      if(Verbosity::level >= Verbosity::FULL)
      {
        std::cerr << "DynamicGBWT::resize(): All nodes apart from the endmarker are empty" << std::endl;
      }
      new_offset = 0; new_sigma = 1;
    }
    else
    {
      new_offset += head_empty;
      for(comp_type tail = this->effective() - 1; this->bwt[tail].empty(); tail--) { new_sigma--; }
    }
  }

  this->forceResize(new_offset, new_sigma);
}

void
DynamicGBWT::forceResize(size_type new_offset, size_type new_sigma)
{
  if(new_offset == this->header.offset && new_sigma == this->sigma()) { return; }

  if(Verbosity::level >= Verbosity::FULL)
  {
    if(new_offset != this->header.offset)
    {
      std::cerr << "DynamicGBWT::resize(): Changing alphabet offset to " << new_offset << std::endl;
    }
    if(new_sigma != this->sigma())
    {
      std::cerr << "DynamicGBWT::resize(): Changing alphabet size to " << new_sigma << std::endl;
    }
  }

  std::vector<DynamicRecord> new_bwt(new_sigma - new_offset);
  if(this->effective() > 0 && new_bwt.size() > 0) { new_bwt.front().swap(this->bwt.front()); }
  comp_type first_comp = 1 + (new_offset > this->header.offset ? new_offset - this->header.offset : 0);
  comp_type comp_tail = std::min(first_comp + new_bwt.size() - 1, this->effective());
  for(comp_type comp = first_comp; comp < comp_tail; comp++)
  {
    new_bwt[comp + this->header.offset - new_offset].swap(this->bwt[comp]);
  }
  this->bwt.swap(new_bwt);
  this->header.offset = new_offset; this->header.alphabet_size = new_sigma;
}

void
DynamicGBWT::recode()
{
  if(Verbosity::level >= Verbosity::FULL)
  {
    std::cerr << "DynamicGBWT::recode(): Sorting the outgoing edges" << std::endl;
  }

  for(comp_type comp = 0; comp < this->effective(); comp++) { this->bwt[comp].recode(); }
}

void
DynamicGBWT::rebuildIncoming()
{
  // Clear the existing incoming edges.
  for(comp_type comp = 0; comp < this->effective(); comp++) { this->bwt[comp].incoming.clear(); }

  // Rebuild them from record bodies and outgoing edges.
  for(comp_type comp = 0; comp < this->effective(); comp++)
  {
    DynamicRecord& current = this->bwt[comp];
    std::vector<size_type> counts(current.outdegree());
    for(run_type run : current.body) { counts[run.first] += run.second; }
    for(rank_type outrank = 0; outrank < current.outdegree(); outrank++)
    {
      if(current.successor(outrank) != ENDMARKER)
      {
        DynamicRecord& successor = this->record(current.successor(outrank));
        successor.addIncoming(edge_type(this->toNode(comp), counts[outrank]));
      }
    }
  }
}

void
DynamicGBWT::rebuildOutgoing()
{
  for(comp_type comp = 0; comp < this->effective(); comp++)
  {
    node_type curr = this->toNode(comp);
    DynamicRecord& current = this->record(curr);
    size_type offset = 0;
    for(rank_type inrank = 0; inrank < current.indegree(); inrank++)
    {
      DynamicRecord& predecessor = this->record(current.predecessor(inrank));
      rank_type outrank = predecessor.edgeTo(curr);
      if(outrank >= predecessor.outdegree())
      {
        std::cerr << "DynamicGBWT::rebuildOutgoing(): Outgoing edge from " << current.predecessor(inrank)
                  << " to " << curr << " not found" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      predecessor.offset(outrank) = offset;
      offset += current.count(inrank);
    }
  }
}

//------------------------------------------------------------------------------

/*
  Support functions for index construction.
*/

void
swapBody(DynamicRecord& record, RunMerger& merger)
{
  merger.flush();
  merger.runs.swap(record.body);
  std::swap(merger.total_size, record.body_size);
}

/*
  Process ranges of sequences sharing the same 'curr' node.
  - Add the outgoing edge (curr, next) if necessary.
  - Add sample (offset, id) if iteration % sample_interval == 0 or next == ENDMARKER.
  - Insert the 'next' node into position 'offset' in the body.
  - Set 'offset' to rank(next) within the record.
  - Update the predecessor count of 'curr' in the incoming edges of 'next'.

  Note: Iteration is the distance from the initial endmarker to the next position.
  We do not maintain incoming edges to the endmarker, because it can be expensive
  and because searching with the endmarker does not work in a multi-string BWT.
*/

void
updateRecords(DynamicGBWT& gbwt, std::vector<Sequence>& seqs, size_type iteration, size_type sample_interval,
  std::unique_ptr<std::unordered_map<node_type, size_type>>& endmarker_edges)
{
  for(size_type i = 0; i < seqs.size(); )
  {
    node_type curr = seqs[i].curr;
    DynamicRecord& current = gbwt.record(curr);
    RunMerger new_body(current.outdegree());
    std::vector<sample_type> new_samples;
    std::vector<run_type>::iterator iter = current.body.begin();
    std::vector<sample_type>::iterator sample_iter = current.ids.begin();
    size_type insert_count = 0;
    while(i < seqs.size() && seqs[i].curr == curr)
    {
      // Determine the edge to the next node or add it if it does not exist.
      rank_type outrank = 0;
      if(curr == ENDMARKER)
      {
        auto iter = endmarker_edges->find(seqs[i].next);
        if(iter != endmarker_edges->end()) { outrank = iter->second; }
        else
        {
          outrank = current.outdegree();
          current.outgoing.push_back(edge_type(seqs[i].next, 0));
          new_body.addEdge();
          (*endmarker_edges)[seqs[i].next] = outrank;
        }
      }
      else
      {
        outrank = current.edgeToLinear(seqs[i].next);
        if(outrank >= current.outdegree())
        {
          current.outgoing.push_back(edge_type(seqs[i].next, 0));
          new_body.addEdge();
        }
      }
      // Add old runs until 'offset'.
      while(new_body.size() < seqs[i].offset)
      {
        if(iter->second <= seqs[i].offset - new_body.size()) { new_body.insert(*iter); ++iter; }
        else
        {
          run_type temp(iter->first, seqs[i].offset - new_body.size());
          new_body.insert(temp);
          iter->second -= temp.second;
        }
      }
      // Add old samples until 'offset'.
      while(sample_iter != current.ids.end() && sample_iter->first + insert_count < seqs[i].offset)
      {
        new_samples.push_back(sample_type(sample_iter->first + insert_count, sample_iter->second));
        ++sample_iter;
      }
      if(iteration % sample_interval == 0 || seqs[i].next == ENDMARKER)  // Sample sequence id.
      {
        new_samples.push_back(sample_type(seqs[i].offset, seqs[i].id));
      }
      seqs[i].offset = new_body.counts[outrank]; // rank(next) within the record.
      new_body.insert(outrank); insert_count++;
      if(seqs[i].next != ENDMARKER)  // The endmarker does not have incoming edges.
      {
        gbwt.record(seqs[i].next).increment(curr);
      }
      i++;
    }
    while(iter != current.body.end()) // Add the rest of the old body.
    {
      new_body.insert(*iter); ++iter;
    }
    while(sample_iter != current.ids.end()) // Add the rest of the old samples.
    {
      new_samples.push_back(sample_type(sample_iter->first + insert_count, sample_iter->second));
      ++sample_iter;
    }
    swapBody(current, new_body);
    current.ids.swap(new_samples);
  }
  gbwt.header.size += seqs.size();
}

/*
  Compute the source offset for each sequence at the next position, assuming that the
  records have been sorted by the node at the current position.
*/

void
nextPosition(std::vector<Sequence>& seqs, const text_type&)
{
  for(Sequence& seq : seqs) { seq.pos++; }
}

void
nextPosition(std::vector<Sequence>& seqs, const vector_type&)
{
  for(Sequence& seq : seqs) { seq.pos++; }
}

void
nextPosition(std::vector<Sequence>& seqs, const GBWT& source)
{
  for(size_type i = 0; i < seqs.size(); )
  {
    node_type curr = seqs[i].curr;
    const CompressedRecord current = source.record(curr);
    CompressedRecordFullIterator iter(current);
    while(i < seqs.size() && seqs[i].curr == curr)
    {
      seqs[i].pos = iter.rankAt(seqs[i].pos);
      i++;
    }
  }
}

void
nextPosition(std::vector<Sequence>& seqs, const DynamicGBWT& source)
{
  for(size_type i = 0; i < seqs.size(); )
  {
    node_type curr = seqs[i].curr;
    const DynamicRecord& current = source.record(curr);
    std::vector<run_type>::const_iterator iter = current.body.begin();
    std::vector<edge_type> result(current.outgoing);
    size_type record_offset = iter->second; result[iter->first].second += iter->second;
    while(i < seqs.size() && seqs[i].curr == curr)
    {
      while(record_offset <= seqs[i].pos)
      {
        ++iter; record_offset += iter->second;
        result[iter->first].second += iter->second;
      }
      seqs[i].pos = result[iter->first].second - (record_offset - seqs[i].pos);
      i++;
    }
  }
}

/*
  Sort the sequences for the next iteration and remove the ones that have reached the endmarker.
  Note that sorting by (next, curr, offset) now is equivalent to sorting by (curr, offset) in the
  next interation.
*/

void
sortSequences(std::vector<Sequence>& seqs)
{
  sequentialSort(seqs.begin(), seqs.end());
  size_type head = 0;
  while(head < seqs.size() && seqs[head].next == ENDMARKER) { head++; }
  if(head > 0)
  {
    for(size_type j = 0; head + j < seqs.size(); j++) { seqs[j] = seqs[head + j]; }
    seqs.resize(seqs.size() - head);
  }
}

/*
  Rebuild the edge offsets in the outgoing edges to each 'next' node. The offsets will be
  valid after the insertions in the next iteration.

  Then add the rebuilt edge offsets to sequence offsets, which have been rank(next)
  within the current record until now.
*/

void
rebuildOffsets(DynamicGBWT& gbwt, std::vector<Sequence>& seqs,
  std::unique_ptr<std::unordered_map<node_type, size_type>>& endmarker_edges)
{
  node_type next = gbwt.sigma();
  for(const Sequence& seq : seqs)
  {
    if(seq.next == next) { continue; }
    next = seq.next;
    size_type offset = 0;
    for(edge_type inedge : gbwt.record(next).incoming)
    {
      DynamicRecord& predecessor = gbwt.record(inedge.first);
      rank_type outrank = (inedge.first == ENDMARKER ? (*endmarker_edges)[next] : predecessor.edgeToLinear(next));
      predecessor.offset(outrank) = offset;
      offset += inedge.second;
    }
  }

  for(Sequence& seq : seqs)
  {
    const DynamicRecord& current = gbwt.record(seq.curr);
    rank_type outrank = (seq.curr == ENDMARKER ? (*endmarker_edges)[seq.next] : current.edgeToLinear(seq.next));
    seq.offset += current.offset(outrank);
  }
}

/*
  Move each sequence to the next position, assuming that the source offset has been
  computed earlier and that the sequences have been sorted by the node at the next
  position.
*/

void
advancePosition(std::vector<Sequence>& seqs, const text_type& text)
{
  for(Sequence& seq : seqs) { seq.curr = seq.next; seq.next = text[seq.pos]; }
}

void
advancePosition(std::vector<Sequence>& seqs, const vector_type& text)
{
  for(Sequence& seq : seqs) { seq.curr = seq.next; seq.next = text[seq.pos]; }
}

void
advancePosition(std::vector<Sequence>& seqs, const GBWT& source)
{
  // FIXME We could optimize further by storing the next position.
  for(size_type i = 0; i < seqs.size(); )
  {
    node_type curr = seqs[i].next;
    const CompressedRecord current = source.record(curr);
    CompressedRecordIterator iter(current);
    while(i < seqs.size() && seqs[i].next == curr)
    {
      seqs[i].curr = seqs[i].next;
      while(iter.offset() <= seqs[i].pos) { ++iter; }
      seqs[i].next = current.successor(iter->first);
      i++;
    }
  }
}

void
advancePosition(std::vector<Sequence>& seqs, const DynamicGBWT& source)
{
  // FIXME We could optimize further by storing the next position.
  for(size_type i = 0; i < seqs.size(); )
  {
    node_type curr = seqs[i].next;
    const DynamicRecord& current = source.record(curr);
    std::vector<run_type>::const_iterator iter = current.body.begin();
    size_type offset = iter->second;
    while(i < seqs.size() && seqs[i].next == curr)
    {
      seqs[i].curr = seqs[i].next;
      while(offset <= seqs[i].pos) { ++iter; offset += iter->second; }
      seqs[i].next = current.successor(iter->first);
      i++;
    }
  }
}

/*
  Insert the sequences from the source to the GBWT. Maintains an invariant that
  the sequences are sorted by (curr, offset).

  Update the bidirectional flag in the header before calling this.
*/

template<class Source>
size_type
insert(DynamicGBWT& gbwt, std::vector<Sequence>& seqs, const Source& source, size_type sample_interval)
{
  // Sanity check sample interval only here.
  if(sample_interval == 0) { sample_interval = std::numeric_limits<size_type>::max(); }

  // The outgoing edges are not expected to be sorted during construction. As the endmarker
  // may have millions of outgoing edges, we need a faster way of mapping destination nodes
  // to edges.
  std::unique_ptr<std::unordered_map<node_type, size_type>> endmarker_edges(new std::unordered_map<node_type, size_type>);
  if(gbwt.sigma() > 0)
  {
    const DynamicRecord& endmarker = gbwt.record(ENDMARKER);
    for(rank_type outrank = 0; outrank < endmarker.outdegree(); outrank++)
    {
      (*endmarker_edges)[endmarker.successor(outrank)] = outrank;
    }
  }

  for(size_type iterations = 1; ; iterations++)
  {
    updateRecords(gbwt, seqs, iterations, sample_interval, endmarker_edges); // Insert the next nodes into the GBWT.
    nextPosition(seqs, source); // Determine the next position for each sequence.
    sortSequences(seqs); // Sort for the next iteration and remove the ones that have finished.
    if(seqs.empty()) { return iterations; }
    rebuildOffsets(gbwt, seqs, endmarker_edges); // Rebuild offsets in outgoing edges and sequences.
    advancePosition(seqs, source); // Move the sequences to the next position.
  }
}

//------------------------------------------------------------------------------

/*
  Insert a batch of sequences with ids (in the current insertion) starting from 'start_id'.
  The template parameter should be an integer vector. Because resizing text_type always
  causes a reallocation, 'text_length' is used to pass the actual length of the text.
  This function assumes that text.size() >= text_length.

  insertBatch() leaves the index in a state where the outgoing edges are not necessarily
  sorted. Subsequent insertBatch() calls are possible, but the index cannot be used
  without calling recode().

  Flag has_both_orientations tells whether the text contains both orientations of each
  sequence.
*/

template<class IntegerVector>
void
insertBatch(DynamicGBWT& index, const IntegerVector& text, size_type text_length,
            size_type start_id, bool has_both_orientations, size_type sample_interval)
{
  double start = readTimer();

  if(text_length == 0) { return; }
  if(text[text_length - 1] != ENDMARKER)
  {
    std::cerr << "insertBatch(): The text must end with an endmarker" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if(!has_both_orientations) { index.header.unset(GBWTHeader::FLAG_BIDIRECTIONAL); }

  /*
    Find the start of each sequence and initialize the sequence objects at the endmarker node.
    Increase alphabet size and decrease offset if necessary.
  */
  bool seq_start = true;
  node_type min_node = (index.empty() ? std::numeric_limits<node_type>::max() : index.header.offset + 1);
  node_type max_node = (index.empty() ? 0 : index.sigma() - 1);
  std::vector<Sequence> seqs;
  for(size_type i = 0; i < text_length; i++)
  {
    if(seq_start)
    {
      seqs.push_back(Sequence(text, i, index.sequences()));
      seq_start = false; index.header.sequences++;
    }
    if(text[i] == ENDMARKER) { seq_start = true; }
    else { min_node = std::min(static_cast<node_type>(text[i]), min_node); }
    max_node = std::max(static_cast<node_type>(text[i]), max_node);
  }
  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    std::cerr << "insertBatch(): Inserting sequences " << start_id << " to " << (start_id + seqs.size() - 1) << std::endl;
  }
  if(max_node == 0) { min_node = 1; } // No real nodes, setting offset to 0.
  index.resize(min_node - 1, max_node + 1);

  // Insert the sequences.
  size_type iterations = gbwt::insert(index, seqs, text, sample_interval);
  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    double seconds = readTimer() - start;
    std::cerr << "insertBatch(): " << iterations << " iterations in " << seconds << " seconds" << std::endl;
  }
}

//------------------------------------------------------------------------------

void
DynamicGBWT::insert(const text_type& text, bool has_both_orientations, size_type sample_interval)
{
  if(text.empty())
  {
    if(Verbosity::level >= Verbosity::FULL)
    {
      std::cerr << "DynamicGBWT::insert(): The input text is empty" << std::endl;
    }
    return;
  }
  gbwt::insertBatch(*this, text, text.size(), 0, has_both_orientations, sample_interval);
  this->recode();
}

void
DynamicGBWT::insert(const text_type& text, size_type text_length, bool has_both_orientations, size_type sample_interval)
{
  if(text_length == 0)
  {
    if(Verbosity::level >= Verbosity::FULL)
    {
      std::cerr << "DynamicGBWT::insert(): The input text is empty" << std::endl;
    }
    return;
  }
  if(text_length > text.size())
  {
    std::cerr << "DynamicGBWT::insert(): Specified text length is larger than container size" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  gbwt::insertBatch(*this, text, text_length, 0, has_both_orientations, sample_interval);
  this->recode();
}

void
DynamicGBWT::insert(const vector_type& text, bool has_both_orientations, size_type sample_interval)
{
  if(text.empty())
  {
    if(Verbosity::level >= Verbosity::FULL)
    {
      std::cerr << "DynamicGBWT::insert(): The input text is empty" << std::endl;
    }
    return;
  }
  gbwt::insertBatch(*this, text, text.size(), 0, has_both_orientations, sample_interval);
  this->recode();
}

void
DynamicGBWT::insert(text_buffer_type& text, size_type batch_size, bool both_orientations, size_type sample_interval)
{
  double start = readTimer();

  if(text.size() == 0)
  {
    if(Verbosity::level >= Verbosity::FULL)
    {
      std::cerr << "DynamicGBWT::insert(): The input text is empty" << std::endl;
    }
    return;
  }
  if(batch_size == 0) { batch_size = text.size(); }
  size_type old_sequences = this->sequences();

  // Create a builder using this index.
  GBWTBuilder builder(text.width(), batch_size, sample_interval);
  builder.swapIndex(*this);

  // Insert all sequences.
  vector_type sequence;
  for(size_type node : text)
  {
    if(node == ENDMARKER) { builder.insert(sequence, both_orientations); sequence.clear(); }
    else { sequence.push_back(node); }
  }
  if(!(sequence.empty())) { builder.insert(sequence, both_orientations); sequence.clear(); }

  // Finish the construction and get the index contents back.
  builder.finish();
  builder.swapIndex(*this);

  if(Verbosity::level >= Verbosity::BASIC)
  {
    double seconds = readTimer() - start;
    std::cerr << "DynamicGBWT::insert(): Inserted " << (this->sequences() - old_sequences)
              << " sequences of total length " << text.size()
              << " in " << seconds << " seconds" << std::endl;
  }
}

//------------------------------------------------------------------------------

void
removePositions(DynamicRecord& record, node_type node, const std::vector<edge_type>& ra, size_type& rank_offset)
{
  if(ra[rank_offset].first != node) { return; }

  size_type body_offset = 0, body_tail = 0;
  size_type bwt_offset = (record.body.empty() ? 0 : record.body.front().second);
  size_type sample_offset = 0, sample_tail = 0;
  size_type removed = 0;

  // Process the removed positions in this record.
  while(rank_offset < ra.size() && ra[rank_offset].first == node)
  {
    // Copy runs that are before the next removed position. Then adjust the run
    // covering the removed position and skip it if it becomes empty.
    // We can safely assume that the run covering the position exists.
    while(bwt_offset <= ra[rank_offset].second)
    {
      record.body[body_tail] = record.body[body_offset];
      body_tail++; body_offset++;
      bwt_offset += record.body[body_offset].second;
    }
    record.body[body_offset].second--;
    if(record.body[body_offset].second == 0)
    {
      body_offset++;
      if(body_offset < record.body.size()) { bwt_offset += record.body[body_offset].second; }
    }
    record.body_size--;

    // Update the samples that are before the next removed position.
    // Then delete the sample that may cover the position.
    while(sample_offset < record.ids.size() && record.ids[sample_offset].first < ra[rank_offset].second)
    {
      record.ids[sample_tail].first = record.ids[sample_offset].first - removed;
      record.ids[sample_tail].second = record.ids[sample_offset].second;
      sample_tail++; sample_offset++;
    }
    if(sample_offset < record.ids.size() && record.ids[sample_offset].first == ra[rank_offset].second)
    {
      sample_offset++;
    }

    rank_offset++;
    removed++;
  }

  // Process the tail of the record and resize the arrays.
  while(body_offset < record.body.size())
  {
    record.body[body_tail] = record.body[body_offset];
    body_tail++; body_offset++;
  }
  record.body.resize(body_tail);
  while(sample_offset < record.ids.size())
  {
    record.ids[sample_tail].first = record.ids[sample_offset].first - removed;
    record.ids[sample_tail].second = record.ids[sample_offset].second;
    sample_tail++; sample_offset++;
  }
  record.ids.resize(sample_tail);

  // Remove unused outgoing edges.
  record.removeUnusedEdges();
}

void
updateSamples(DynamicRecord& record, const sdsl::bit_vector::rank_1_type& remaining_rank)
{
  for(sample_type& sample : record.ids)
  {
    sample.second = remaining_rank(sample.second);
  }
}

size_type
DynamicGBWT::remove(const std::vector<size_type>& seq_ids, size_type chunk_size)
{
  double start = readTimer();

  if(seq_ids.empty())
  {
    if(Verbosity::level >= Verbosity::FULL)
    {
      std::cerr << "DynamicGBWT::remove(): No sequences to remove" << std::endl;
    }
    return 0;
  }

  std::vector<size_type> to_remove;
  for(size_type seq_id : seq_ids)
  {
    if(this->bidirectional())
    {
      to_remove.push_back(Path::encode(seq_id, false));
      to_remove.push_back(Path::encode(seq_id, true));
    }
    else { to_remove.push_back(seq_id); }
  }
  removeDuplicates(to_remove, false);
  for(size_type seq_id : to_remove)
  {
    if(seq_id >= this->sequences())
    {
      std::cerr << "DynamicGBWT::remove(): Invalid sequence id: " << seq_id << std::endl;
      return 0;
    }
  }

  // Build a bitvector of the sequences that will remain.
  sdsl::bit_vector remaining(this->sequences(), 1);
  for(size_type seq_id : to_remove) { remaining[seq_id] = 0; }
  sdsl::bit_vector::rank_1_type remaining_rank;
  sdsl::util::init_support(remaining_rank, &remaining);  

  // Build the rank array.
  double ra_start = readTimer();
  std::vector<edge_type> ra;
  int threads = 1; //omp_get_max_threads();
  std::vector<std::vector<edge_type>> buffers(threads);
  DecompressedRecord fast_endmarker = this->endmarker();
  //#pragma omp parallel for schedule(dynamic, chunk_size)
  for(size_type i = 0; i < to_remove.size(); i++)
  {
    int thread_id = 0;//omp_get_thread_num();
    std::vector<edge_type>& buffer = buffers[thread_id];
    buffer.push_back(edge_type(ENDMARKER, to_remove[i]));
    edge_type pos = fast_endmarker.LF(to_remove[i]);
    while(pos.first != ENDMARKER)
    {
      buffer.push_back(pos);
      pos = this->LF(pos);
    }
    //#pragma omp critical
    {
      ra.insert(ra.end(), buffer.begin(), buffer.end());
    }
    buffer.clear();
  }
  parallelQuickSort(ra.begin(), ra.end());
  if(Verbosity::level >= Verbosity::BASIC)
  {
    double seconds = readTimer() - ra_start;
    std::cerr << "DynamicGBWT::remove(): Rank array built in " << seconds << " seconds" << std::endl;
  }

  // Remove the sequences.
  double update_start = readTimer();
  size_type offset = 0; // Current offset in ra.
  for(comp_type comp = 0; comp < this->effective(); comp++)
  {
    node_type node = this->toNode(comp);
    if(offset < ra.size()) { removePositions(this->record(node), node, ra, offset); }
    updateSamples(this->record(node), remaining_rank);
  }
  if(Verbosity::level >= Verbosity::BASIC)
  {
    double seconds = readTimer() - update_start;
    std::cerr << "DynamicGBWT::remove(): Records updated in " << seconds << " seconds" << std::endl;
  }

  // Update the header.
  this->header.sequences -= to_remove.size();
  this->header.size -= ra.size();

  // Remove empty nodes from both ends and rebuild the edges from the BWT.
  this->resize();
  if(Verbosity::level >= Verbosity::FULL)
  {
    std::cerr << "DynamicGBWT::remove(): Rebuilding the edges" << std::endl;
  }
  this->rebuildIncoming();
  this->rebuildOutgoing();

  if(Verbosity::level >= Verbosity::BASIC)
  {
    double seconds = readTimer() - start;
    std::cerr << "DynamicGBWT::remove(): Removed " << to_remove.size() << " sequences of total length "
              << ra.size() << " in " << seconds << " seconds" << std::endl;
  }

  return ra.size();
}

size_type
DynamicGBWT::remove(size_type seq_id, size_type chunk_size)
{
  return this->remove(std::vector<size_type>(seq_id, 1), chunk_size);
}  

//------------------------------------------------------------------------------

template<class GBWTType>
void
mergeMetadata(DynamicGBWT& index, const GBWTType& source, bool index_was_empty)
{
  // One of the GBWTs is empty.
  if(source.empty()) { return; }
  if(index_was_empty)
  {
    if(source.hasMetadata())
    {
      index.addMetadata();
      index.metadata = source.metadata;
    }
    return;
  }

  // Non-empty GBWTs.
  if(index.hasMetadata() && source.hasMetadata())
  {
    index.metadata.merge(source.metadata, false, true); // Different samples, same contigs.
  }
  else if(index.hasMetadata())
  {
    if(Verbosity::level >= Verbosity::BASIC)
    {
      std::cerr << "DynamicGBWT::merge(): Clearing metadata: no metadata in the other GBWT" << std::endl;
    }
    index.clearMetadata();
  }
}

void
DynamicGBWT::merge(const GBWT& source, size_type batch_size, size_type sample_interval)
{
  double start = readTimer();

  if(source.empty())
  {
    if(Verbosity::level >= Verbosity::FULL)
    {
      std::cerr << "DynamicGBWT::merge(): The input GBWT is empty" << std::endl;
    }
    return;
  }
  bool index_was_empty = this->empty();

  // The merged index is bidirectional only if both indexes are bidirectional.
  if(!(source.bidirectional())) { this->header.unset(GBWTHeader::FLAG_BIDIRECTIONAL); }

  // Increase alphabet size and decrease offset if necessary.
  if(batch_size == 0) { batch_size = source.sequences(); }
  this->resize(source.header.offset, source.sigma());

  // Insert the sequences in batches.
  const DecompressedRecord& endmarker = source.endmarker();
  size_type source_id = 0;
  while(source_id < source.sequences())
  {
    double batch_start = readTimer();
    size_type limit = std::min(source_id + batch_size, source.sequences());
    std::vector<Sequence> seqs; seqs.reserve(limit - source_id);
    while(source_id < limit)  // Create the new sequence iterators.
    {
      seqs.emplace_back(endmarker[source_id], this->sequences(), source_id);
      this->header.sequences++; source_id++;
    }
    if(Verbosity::level >= Verbosity::EXTENDED)
    {
      std::cerr << "DynamicGBWT::merge(): Inserting sequences " << (source_id - seqs.size())
                << " to " << (source_id - 1) << std::endl;
    }
    size_type iterations = gbwt::insert(*this, seqs, source, sample_interval);
    if(Verbosity::level >= Verbosity::EXTENDED)
    {
      double seconds = readTimer() - batch_start;
      std::cerr << "DynamicGBWT::merge(): " << iterations << " iterations in " << seconds << " seconds" << std::endl;
    }
  }

  // Finally sort the outgoing edges.
  this->recode();

  // Merge the metadata.
  mergeMetadata(*this, source, index_was_empty);

  if(Verbosity::level >= Verbosity::BASIC)
  {
    double seconds = readTimer() - start;
    std::cerr << "DynamicGBWT::merge(): Inserted " << source.sequences() << " sequences of total length "
              << source.size() << " in " << seconds << " seconds" << std::endl;
  }
}

//------------------------------------------------------------------------------

/*
  Builds the rank array of 'right' relative to 'left'. When the rank array RA is sorted, we
  know that there will be RA[i] characters from left_BWT before right_BWT[i].

  We assume that 'buffers' has been set to use the same number of threads as OpenMP. The
  sequences from 'right' will get identifiers after those from 'left'.
*/

void
buildRA(const DynamicGBWT& left, const DynamicGBWT& right, MergeBuffers& buffers)
{
  DecompressedRecord right_endmarker = right.endmarker();

  //#pragma omp parallel for schedule(dynamic, buffers.parameters.chunk_size)
  for(size_type sequence = 0; sequence < right.sequences(); sequence++)
  {
    // The new sequence will be after all existing sequences in 'left'.
    size_type thread = 0;//omp_get_thread_num();
    buffers.insert(edge_type(ENDMARKER, left.sequences()), thread);

    // Computing LF() at the endmarker can be expensive, so we do it using the incoming
    // edges at the destination node instead. If the sequence is empty, countUntil().
    // will return 0, because the endmarker does not have incoming edges.
    edge_type right_pos = right_endmarker.LF(sequence);
    edge_type left_pos(right_pos.first, left.record(right_pos.first).countUntil(ENDMARKER));

    // Main loop. We can assume that the positions are always valid.
    while(right_pos.first != ENDMARKER)
    {
      buffers.insert(left_pos, thread);
      right_pos = right.LF(right_pos);
      left_pos = edge_type(right_pos.first, left.fullLF(left_pos, right_pos.first));
    }
  }

  buffers.flush();
}

/*
  Merge the outgoing edges. We can ignore the offsets for now.
*/

std::vector<edge_type>
mergeOutgoing(const std::vector<edge_type>& left, const std::vector<edge_type>& right)
{
  std::vector<edge_type> result;

  auto left_iter = left.begin();
  auto right_iter = right.begin();
  while(left_iter != left.end() && right_iter != right.end())
  {
    if(left_iter->first == right_iter->first)
    {
      result.push_back(*left_iter); ++left_iter; ++right_iter;
    }
    else if(left_iter->first < right_iter->first)
    {
      result.push_back(*left_iter); ++left_iter;
    }
    else
    {
      result.push_back(*right_iter); ++right_iter;
    }
  }
  while(left_iter != left.end()) { result.push_back(*left_iter); ++left_iter; }
  while(right_iter != right.end()) { result.push_back(*right_iter); ++right_iter; }

  return result;
}

/*
  Recode the run from 'source' using 'outgoing' as the list of outgoing edges.
*/

run_type
recodeRun(run_type run, const DynamicRecord& source, const std::vector<edge_type>& outgoing)
{
  run.first = edgeTo(source.successor(run.first), outgoing);
  return run;
}

/*
  Merges 'right' into 'left' using the rank array. We need the node identifier for finding
  the correct range in the rank array and the number of sequences in 'left' for updating the
  samples from 'right'.
*/

void
mergeRecords(DynamicRecord& left, const DynamicRecord& right, ProducerBuffer<RankArray>& ra, node_type node, size_type left_sequences)
{
  // Rebuild the record using these structures.
  std::vector<edge_type> new_outgoing = mergeOutgoing(left.outgoing, right.outgoing);
  RunMerger new_body(new_outgoing.size());
  std::vector<sample_type> new_samples;

  // Rank array contains the number of elements from 'left' before each element from 'right'.
  auto left_iter = left.body.begin(); auto right_iter = right.body.begin();
  auto left_sample_iter = left.ids.begin(); auto right_sample_iter = right.ids.begin();
  run_type left_run(0, 0), right_run(0, 0);
  size_type insert_count = 0;
  while(!(ra.end()) && ra->first == node)
  {
    // Add runs from 'left'.
    while(new_body.size() < ra->second + insert_count)
    {
      if(left_run.second == 0) { left_run = recodeRun(*left_iter, left, new_outgoing); ++left_iter; }
      size_type insert_length = std::min(static_cast<size_type>(left_run.second), ra->second + insert_count - new_body.size());
      new_body.insert(run_type(left_run.first, insert_length));
      left_run.second -= insert_length;
    }
    // Add samples from 'left'.
    while(left_sample_iter != left.ids.end() && left_sample_iter->first < ra->second)
    {
      new_samples.emplace_back(left_sample_iter->first + insert_count, left_sample_iter->second);
      ++left_sample_iter;
    }
    // Add a single value and the possible sample from 'right'.
    if(right_run.second == 0) { right_run = recodeRun(*right_iter, right, new_outgoing); ++right_iter; }
    new_body.insert(right_run.first);
    if(right_sample_iter != right.ids.end() && right_sample_iter->first == insert_count)
    {
      new_samples.emplace_back(ra->second + insert_count, right_sample_iter->second + left_sequences);
      ++right_sample_iter;
    }
    right_run.second--; insert_count++;
    ++ra;
  }

  // Add the remaining runs from 'left'.
  if(left_run.second > 0) { new_body.insert(left_run); }
  while(left_iter != left.body.end())
  {
    new_body.insert(recodeRun(*left_iter, left, new_outgoing)); ++left_iter;
  }
  // Add the remaining samples from 'left'.
  while(left_sample_iter != left.ids.end())
  {
    new_samples.emplace_back(left_sample_iter->first + insert_count, left_sample_iter->second);
    ++left_sample_iter;
  }

  // Use the new data in 'left'.
  swapBody(left, new_body);
  left.outgoing.swap(new_outgoing);
  left.ids.swap(new_samples);
}

//------------------------------------------------------------------------------

void
DynamicGBWT::merge(const DynamicGBWT& source, const MergeParameters& parameters)
{
  double start = readTimer();

  if(source.empty())
  {
    if(Verbosity::level >= Verbosity::FULL)
    {
      std::cerr << "DynamicGBWT::merge(): The input GBWT is empty" << std::endl;
    }
    return;
  }
  bool index_was_empty = this->empty();

  // The merged index is bidirectional only if both indexes are bidirectional.
  if(!(source.bidirectional())) { this->header.unset(GBWTHeader::FLAG_BIDIRECTIONAL); }

  // Increase alphabet size and decrease offset if necessary.
  this->resize(source.header.offset, source.sigma());

  // Determine the node ranges for merge jobs.
  std::vector<range_type> node_ranges = Range::partition(range_type(0, this->effective() - 1), parameters.merge_jobs);
  for(range_type& range : node_ranges)
  {
    range.first = this->toNode(range.first);
    range.second = this->toNode(range.second);
  }

  // Build the rank array.
  double ra_start = readTimer();
  int threads = 1; //omp_get_max_threads();
  MergeBuffers mb(source.size(), threads, parameters, node_ranges);
  buildRA(*this, source, mb);
  if(Verbosity::level >= Verbosity::BASIC)
  {
    double seconds = readTimer() - ra_start;
    std::cerr << "DynamicGBWT::merge(): Rank array built in " << seconds << " seconds" << std::endl;
  }

  // Merge the records.
  double merge_start = readTimer();
  //#pragma omp parallel for schedule(static)
  for(size_type job = 0; job < node_ranges.size(); job++)
  {
    ProducerBuffer<RankArray> ra(*(mb.ra[job]));
    for(node_type node = node_ranges[job].first; node <= node_ranges[job].second; node++)
    {
      if(!(source.contains(node))) { continue; }
      mergeRecords(this->record(node), source.record(node), ra, node, this->sequences());
    }
  }
  if(Verbosity::level >= Verbosity::BASIC)
  {
    double seconds = readTimer() - merge_start;
    std::cerr << "DynamicGBWT::merge(): Records merged in " << seconds << " seconds" << std::endl;
  }

  // Merge the headers. Note that we need the original header for merging the records.
  this->header.size += source.size();
  this->header.sequences += source.sequences();

  // Rebuild the incoming edges from the record bodies and the outgoing edges. Then rebuild
  // the offsets in the outgoing edges from the incoming edges.
  if(Verbosity::level >= Verbosity::FULL)
  {
    std::cerr << "DynamicGBWT::merge(): Rebuilding the edges" << std::endl;
  }
  this->rebuildIncoming();
  this->rebuildOutgoing();

  // Merge the metadata.
  mergeMetadata(*this, source, index_was_empty);

  if(Verbosity::level >= Verbosity::BASIC)
  {
    double seconds = readTimer() - start;
    std::cerr << "DynamicGBWT::merge(): Inserted " << source.sequences() << " sequences of total length "
              << source.size() << " in " << seconds << " seconds" << std::endl;
  }
}

//------------------------------------------------------------------------------

size_type
DynamicGBWT::tryLocate(node_type node, size_type i) const
{
  const DynamicRecord& record = this->record(node);
  for(sample_type sample : record.ids)
  {
    if(sample.first == i) { return sample.second; }
    if(sample.first > i) { break; }
  }
  return invalid_sequence();
}

// FIXME This should really have a common implementation with GBWT::locate(state).
std::vector<size_type>
DynamicGBWT::locate(SearchState state) const
{
  std::vector<size_type> result;
  if(!(this->contains(state))) { return result; }

  // Initialize BWT positions for each offset in the range.
  std::vector<edge_type> positions(state.size());
  for(size_type i = state.range.first; i <= state.range.second; i++)
  {
    positions[i - state.range.first] = edge_type(state.node, i);
  }

  // Continue with LF() until samples have been found for all sequences.
  while(!(positions.empty()))
  {
    size_type tail = 0;
    node_type curr = invalid_node();
    const DynamicRecord* current = nullptr;
    std::vector<sample_type>::const_iterator sample;
    edge_type LF_result;
    range_type LF_range;

    for(size_type i = 0; i < positions.size(); i++)
    {
      if(positions[i].first != curr)  // Node changed.
      {
        curr = positions[i].first; current = &(this->record(curr));
        sample = current->nextSample(positions[i].second);
        LF_range.first = positions[i].second;
        LF_result = current->runLF(positions[i].second, LF_range.second);
      }
      while(sample != current->ids.end() && sample->first < positions[i].second)  // Went past the sample.
      {
        ++sample;
      }
      if(sample == current->ids.end() || sample->first > positions[i].second) // Not sampled.
      {
        if(positions[i].second > LF_range.second) // Went past the existing LF() result.
        {
          LF_range.first = positions[i].second;
          LF_result = current->runLF(positions[i].second, LF_range.second);
        }
        positions[tail] = edge_type(LF_result.first, LF_result.second + positions[i].second - LF_range.first);
        tail++;
      }
      else  // Found a sample.
      {
        result.push_back(sample->second);
      }
    }
    positions.resize(tail);
    sequentialSort(positions.begin(), positions.end());
  }

  removeDuplicates(result, false);
  return result;
}

//------------------------------------------------------------------------------

void
printStatistics(const DynamicGBWT& gbwt, const std::string& name, std::ostream& out)
{
  printHeader(indexType(gbwt), out) << name;
  if(gbwt.bidirectional()) { out << " (bidirectional)"; }
  out << std::endl;
  printHeader("Total length", out) << gbwt.size() << std::endl;
  printHeader("Sequences", out) << gbwt.sequences() << std::endl;
  printHeader("Alphabet size", out) << gbwt.sigma() << std::endl;
  printHeader("Effective", out) << gbwt.effective() << std::endl;
  std::pair<size_type, size_type> runs = gbwt.runs();
  printHeader("Runs", out) << runs.first << " concrete / " << runs.second << " logical" << std::endl;
  printHeader("DA samples", out) << gbwt.samples() << std::endl;
  if(gbwt.hasMetadata())
  {
    printHeader("Metadata", out) << gbwt.metadata << std::endl;
  }
  out << std::endl;
}

std::string
indexType(const DynamicGBWT&)
{
  return "Dynamic GBWT";
}

//------------------------------------------------------------------------------

GBWTBuilder::GBWTBuilder(size_type node_width, size_type buffer_size, size_type sample_interval) :
  input_buffer(buffer_size, 0, node_width), construction_buffer(buffer_size, 0, node_width),
  id_sample_interval(sample_interval),
  input_tail(0), construction_tail(0),
  inserted_sequences(0), batch_sequences(0),
  has_both_orientations(true)
{
}

GBWTBuilder::~GBWTBuilder()
{
  // Wait for the construction thread to finish.
  if(this->builder.joinable()) { this->builder.join(); }
}

void
GBWTBuilder::swapIndex(DynamicGBWT& another_index)
{
  this->index.swap(another_index);
}

void
GBWTBuilder::insert(const vector_type& sequence, bool both_orientations)
{
  size_type space_required = sequence.size() + 1;
  if(both_orientations) { space_required *= 2; }
  this->has_both_orientations &= both_orientations;

  // Flush the buffer if necessary.
  if(this->input_tail + space_required > this->input_buffer.size())
  {
    this->flush();
    if(space_required > this->input_buffer.size()) {
      if(Verbosity::level >= Verbosity::EXTENDED) {
        std::cerr << "GBWTBuilder::insert(): Sequence is too long; increasing buffer size to " << space_required << std::endl;
      }
      this->input_buffer.resize(space_required);
    }
  }

  // Forward orientation.
  for(auto node : sequence) { this->input_buffer[this->input_tail] = node; this->input_tail++; }
  this->input_buffer[this->input_tail] = ENDMARKER; this->input_tail++;
  this->batch_sequences++;

  // Reverse orientation.
  if(both_orientations)
  {
    reversePath(sequence, this->input_buffer, this->input_tail);
    this->input_buffer[this->input_tail] = ENDMARKER; this->input_tail++;
    this->batch_sequences++;
  }
}

void
GBWTBuilder::finish()
{
  // Flush the buffer if necessary.
  this->flush();

  // Wait for the construction thread to finish.
  if(this->builder.joinable()) { this->builder.join(); }

  // Finally recode the index to make it serializable.
  this->index.recode();
}

void
GBWTBuilder::flush()
{
  // Wait for the construction thread to finish.
  if(this->builder.joinable()) { this->builder.join(); }

  // Swap the input buffer and the construction buffer.
  this->input_buffer.swap(this->construction_buffer);
  this->construction_tail = this->input_tail;
  this->input_tail = 0;

  // Launch a new construction thread if necessary.
  if(this->construction_tail > 0)
  {
    this->builder =
      std::thread(gbwt::insertBatch<text_type>, std::ref(this->index), std::cref(this->construction_buffer),
                  this->construction_tail, this->inserted_sequences, this->has_both_orientations, this->id_sample_interval);
    this->inserted_sequences += this->batch_sequences;
    this->batch_sequences = 0;
  }
}

//------------------------------------------------------------------------------

} // namespace gbwt
