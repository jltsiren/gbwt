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

#include <gbwt/dynamic_gbwt.h>
#include <gbwt/internal.h>

namespace gbwt
{

//------------------------------------------------------------------------------

const std::string DynamicGBWT::EXTENSION = ".gbwt";

DynamicGBWT::DynamicGBWT()
{
}

DynamicGBWT::DynamicGBWT(const DynamicGBWT& source)
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
    this->bwt.swap(another.bwt);
  }
}

DynamicGBWT&
DynamicGBWT::operator=(const DynamicGBWT& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

DynamicGBWT&
DynamicGBWT::operator=(DynamicGBWT&& source)
{
  if(this != &source)
  {
    this->header = std::move(source.header);
    this->bwt = std::move(source.bwt);
  }
  return *this;
}

size_type
DynamicGBWT::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->header.serialize(out, child, "header");

  {
    RecordArray array(this->bwt);
    written_bytes += array.serialize(out, child, "bwt");
  }

  {
    DASamples compressed_samples(this->bwt);
    written_bytes += compressed_samples.serialize(out, child, "da_samples");
  }

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
DynamicGBWT::load(std::istream& in)
{
  // Read the header.
  this->header.load(in);
  if(!(this->header.check()))
  {
    std::cerr << "DynamicGBWT::load(): Invalid header: " << this->header << std::endl;
  }
  this->bwt.resize(this->effective());

  // Read and decompress the BWT.
  {
    RecordArray array;
    array.load(in);
    size_type offset = 0;
    for(comp_type comp = 0; comp < this->effective(); comp++)
    {
      size_type limit = array.limit(comp);
      DynamicRecord& current = this->bwt[comp];
      current.clear();

      // Decompress the outgoing edges.
      current.outgoing.resize(ByteCode::read(array.data, offset));
      node_type prev = 0;
      for(edge_type& outedge : current.outgoing)
      {
        outedge.first = ByteCode::read(array.data, offset) + prev;
        prev = outedge.first;
        outedge.second = ByteCode::read(array.data, offset);
      }

      // Decompress the body.
      if(current.outdegree() > 0)
      {
        Run decoder(current.outdegree());
        while(offset < limit)
        {
          run_type run = decoder.read(array.data, offset);
          current.body.push_back(run);
          current.body_size += run.second;
        }
      }
    }
  }

  // Read and decompress the samples.
  {
    DASamples samples;
    samples.load(in);
    sdsl::sd_vector<>::select_1_type offset_select(&(samples.sampled_offsets));
    size_type record_rank = 0, max_rank = samples.record_rank(samples.sampled_records.size());
    size_type record_start = 0;
    size_type sample_rank = 0;
    for(comp_type comp = 0; comp < this->effective(); comp++)
    {
      if(samples.sampled_records[comp] == 0) { continue; }
      DynamicRecord& current = this->bwt[comp];
      size_type limit = (record_rank + 1 < max_rank ? samples.bwt_select(record_rank + 2) : samples.bwt_ranges.size());
      while(sample_rank < samples.size())
      {
        size_type sample_offset = offset_select(sample_rank + 1);
        if(sample_offset >= limit) { break; }
        current.ids.push_back(sample_type(sample_offset - record_start, samples.array[sample_rank]));
        sample_rank++;
      }
      record_rank++; record_start = limit;
    }
  }

  // Rebuild the incoming edges.
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
        successor.addIncoming(edge_type(comp, counts[outrank]));
      }
    }
  }
}

void
DynamicGBWT::copy(const DynamicGBWT& source)
{
  this->header = source.header;
  this->bwt = source.bwt;
}

//------------------------------------------------------------------------------

size_type
DynamicGBWT::runs() const
{
  size_type total = 0;
  for(const DynamicRecord& node : this->bwt) { total += node.runs(); }
  return total;
}

size_type
DynamicGBWT::samples() const
{
  size_type total = 0;
  for(const DynamicRecord& node : this->bwt) { total += node.samples(); }
  return total;
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

  if(new_offset != this->header.offset || new_sigma != this->sigma())
  {
    if(Verbosity::level >= Verbosity::FULL)
    {
      if(new_offset != this->header.offset)
      {
        std::cerr << "DynamicGBWT::resize(): Changing alphabet offset to " << new_offset << std::endl;
      }
      if(new_sigma != this->sigma())
      {
        std::cerr << "DynamicGBWT::resize(): Increasing alphabet size to " << new_sigma << std::endl;
      }
    }

    std::vector<DynamicRecord> new_bwt(new_sigma - new_offset);
    if(this->effective() > 0) { new_bwt[0].swap(this->bwt[0]); }
    for(comp_type comp = 1; comp < this->effective(); comp++)
    {
      new_bwt[comp + this->header.offset - new_offset].swap(this->bwt[comp]);
    }
    this->bwt.swap(new_bwt);
    this->header.offset = new_offset; this->header.alphabet_size = new_sigma;
  }
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
  - Add sample (offset, id) if iteration % SAMPLE_INTERVAL == 0 or next == ENDMARKER.
  - Insert the 'next' node into position 'offset' in the body.
  - Set 'offset' to rank(next) within the record.
  - Update the predecessor count of 'curr' in the incoming edges of 'next'.

  We do not maintain incoming edges to the endmarker, because it can be expensive
  and because searching with the endmarker does not work in a multi-string BWT.
*/

void
updateRecords(DynamicGBWT& gbwt, std::vector<Sequence>& seqs, size_type iteration)
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
      rank_type outrank = current.edgeTo(seqs[i].next);
      if(outrank >= current.outdegree())  // Add edge (curr, next) if it does not exist.
      {
        current.outgoing.push_back(edge_type(seqs[i].next, 0));
        new_body.addEdge();
      }
      while(new_body.size() < seqs[i].offset)  // Add old runs until 'offset'.
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
      if(iteration % DynamicGBWT::SAMPLE_INTERVAL == 0 || seqs[i].next == ENDMARKER)  // Sample sequence id.
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
    current.ids = new_samples;
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
nextPosition(std::vector<Sequence>& seqs, const std::vector<node_type>&)
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
rebuildOffsets(DynamicGBWT& gbwt, std::vector<Sequence>& seqs)
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
      predecessor.offset(predecessor.edgeTo(next)) = offset;
      offset += inedge.second;
    }
  }

  for(Sequence& seq : seqs)
  {
    const DynamicRecord& current = gbwt.record(seq.curr);
    seq.offset += current.offset(current.edgeTo(seq.next));
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
advancePosition(std::vector<Sequence>& seqs, const std::vector<node_type>& text)
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
*/

template<class Source>
size_type
insert(DynamicGBWT& gbwt, std::vector<Sequence>& seqs, const Source& source)
{
  for(size_type iterations = 1; ; iterations++)
  {
    updateRecords(gbwt, seqs, iterations);  // Insert the next nodes into the GBWT.
    nextPosition(seqs, source); // Determine the next position for each sequence.
    sortSequences(seqs);  // Sort for the next iteration and remove the ones that have finished.
    if(seqs.empty()) { return iterations; }
    rebuildOffsets(gbwt, seqs); // Rebuild offsets in outgoing edges and sequences.
    advancePosition(seqs, source);  // Move the sequences to the next position.
  }
}

//------------------------------------------------------------------------------

/*
  Insert a batch of sequences with ids (in the current insertion) starting from 'start_id'.
  The template parameter should be an integer vector.
*/

template<class IntegerVector>
void
insertBatch(DynamicGBWT& index, const IntegerVector& text, size_type start_id)
{
  double start = readTimer();

  if(text.empty()) { return; }
  if(text[text.size() - 1] != ENDMARKER)
  {
    std::cerr << "insertBatch(): The text must end with an endmarker" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  /*
    Find the start of each sequence and initialize the sequence objects at the endmarker node.
    Increase alphabet size and decrease offset if necessary.
  */
  bool seq_start = true;
  node_type min_node = (index.empty() ? ~(node_type)0 : index.header.offset + 1);
  node_type max_node = (index.empty() ? 0 : index.sigma() - 1);
  std::vector<Sequence> seqs;
  for(size_type i = 0; i < text.size(); i++)
  {
    if(seq_start)
    {
      seqs.push_back(Sequence(text, i, index.sequences()));
      seq_start = false; index.header.sequences++;
    }
    if(text[i] == ENDMARKER) { seq_start = true; }
    else { min_node = std::min(text[i], min_node); }
    max_node = std::max(text[i], max_node);
  }
  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    std::cerr << "insertBatch(): Inserting sequences " << start_id << " to " << (start_id + seqs.size() - 1) << std::endl;
  }
  if(max_node == 0) { min_node = 1; } // No real nodes, setting offset to 0.
  index.resize(min_node - 1, max_node + 1);

  // Insert the sequences and sort the outgoing edges.
  size_type iterations = gbwt::insert(index, seqs, text);
  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    double seconds = readTimer() - start;
    std::cerr << "insertBatch(): " << iterations << " iterations in " << seconds << " seconds" << std::endl;
  }
}

//------------------------------------------------------------------------------

void
DynamicGBWT::insert(const text_type& text)
{
  if(text.empty())
  {
    if(Verbosity::level >= Verbosity::FULL)
    {
      std::cerr << "DynamicGBWT::insert(): The input text is empty" << std::endl;
    }
    return;
  }
  gbwt::insertBatch(*this, text, 0);
  this->recode();
}

void
DynamicGBWT::insert(const std::vector<node_type>& text)
{
  if(text.empty())
  {
    if(Verbosity::level >= Verbosity::FULL)
    {
      std::cerr << "DynamicGBWT::insert(): The input text is empty" << std::endl;
    }
    return;
  }
  gbwt::insertBatch(*this, text, 0);
  this->recode();
}

void
DynamicGBWT::insert(text_buffer_type& text, size_type batch_size)
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

  // Find the last endmarker in the batch, read the batch into memory, and insert the sequences.
  size_type old_sequences = this->sequences();
  for(size_type start_offset = 0; start_offset < text.size(); )
  {
    size_type limit = std::min(text.size(), start_offset + batch_size);
    while(limit > start_offset)
    {
      if(text[limit - 1] == ENDMARKER) { break; }
      limit--;
    }
    if(limit <= start_offset)
    {
      std::cerr << "DynamicGBWT::insert(): Cannot find an endmarker in the batch starting from offset " << start_offset << std::endl;
      std::exit(EXIT_FAILURE);
    }
    text_type batch(limit - start_offset, 0, text.width());
    for(size_type i = start_offset; i < limit; i++) { batch[i - start_offset] = text[i]; }
    gbwt::insertBatch(*this, batch, this->sequences() - old_sequences);
    start_offset = limit;
  }

  // Finally sort the outgoing edges.
  this->recode();

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
DynamicGBWT::merge(const GBWT& source, size_type batch_size)
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

  // Increase alphabet size and decrease offset if necessary.
  if(batch_size == 0) { batch_size = source.sequences(); }
  this->resize(source.header.offset, source.sigma());

  // Insert the sequences in batches.
  const CompressedRecord endmarker = source.record(ENDMARKER);
  CompressedRecordIterator iter(endmarker);
  size_type source_id = 0, run_offset = 0;
  while(source_id < source.sequences())
  {
    double batch_start = readTimer();
    size_type limit = std::min(source_id + batch_size, source.sequences());
    std::vector<Sequence> seqs; seqs.reserve(limit - source_id);
    while(source_id < limit)  // Create the new sequence iterators.
    {
      if(run_offset >= iter->second) { ++iter; run_offset = 0; }
      else
      {
        seqs.push_back(Sequence(endmarker.successor(iter->first), this->sequences(), source_id));
        this->header.sequences++; source_id++; run_offset++;
      }
    }
    if(Verbosity::level >= Verbosity::EXTENDED)
    {
      std::cerr << "DynamicGBWT::merge(): Inserting sequences " << (source_id - seqs.size())
                << " to " << (source_id - 1) << std::endl;
    }
    size_type iterations = gbwt::insert(*this, seqs, source);
    if(Verbosity::level >= Verbosity::EXTENDED)
    {
      double seconds = readTimer() - batch_start;
      std::cerr << "DynamicGBWT::merge(): " << iterations << " iterations in " << seconds << " seconds" << std::endl;
    }
  }

  // Finally sort the outgoing edges.
  this->recode();

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
    DynamicRecord current;
    std::vector<sample_type>::const_iterator sample;
    edge_type LF_result;
    range_type LF_range;

    for(size_type i = 0; i < positions.size(); i++)
    {
      if(positions[i].first != curr)  // Node changed.
      {
        curr = positions[i].first; current = this->record(curr);
        sample = current.nextSample(positions[i].second);
        LF_range.first = positions[i].second;
        LF_result = current.runLF(positions[i].second, LF_range.second);
      }
      while(sample != current.ids.end() && sample->first < positions[i].second)  // Went past the sample.
      {
        ++sample;
      }
      if(sample == current.ids.end() || sample->first > positions[i].second) // Not sampled.
      {
        if(positions[i].second > LF_range.second) // Went past the existing LF() result.
        {
          LF_range.first = positions[i].second;
          LF_result = current.runLF(positions[i].second, LF_range.second);
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
printStatistics(const DynamicGBWT& gbwt, const std::string& name)
{
  printHeader("Dynamic GBWT"); std::cout << name << std::endl;
  printHeader("Total length"); std::cout << gbwt.size() << std::endl;
  printHeader("Sequences"); std::cout << gbwt.sequences() << std::endl;
  printHeader("Alphabet size"); std::cout << gbwt.sigma() << std::endl;
  printHeader("Effective"); std::cout << gbwt.effective() << std::endl;
  printHeader("Runs"); std::cout << gbwt.runs() << std::endl;
  printHeader("Samples"); std::cout << gbwt.samples() << std::endl;
  std::cout << std::endl;
}

//------------------------------------------------------------------------------

} // namespace gbwt
