/*
  Copyright (c) 2020, 2021, 2022 Jouni Siren

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

#include <gbwt/fast_locate.h>
#include <gbwt/internal.h>

namespace gbwt
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr std::uint32_t FastLocate::Header::TAG;
constexpr std::uint32_t FastLocate::Header::VERSION;

constexpr size_type FastLocate::NO_POSITION;

//------------------------------------------------------------------------------

// Other class variables.

const std::string FastLocate::EXTENSION = ".ri";

//------------------------------------------------------------------------------

FastLocate::Header::Header() :
  tag(TAG), version(VERSION),
  max_length(1),
  flags(0)
{
}

size_type
FastLocate::Header::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;
  written_bytes += sdsl::write_member(this->tag, out, child, "tag");
  written_bytes += sdsl::write_member(this->version, out, child, "version");
  written_bytes += sdsl::write_member(this->max_length, out, child, "max_length");
  written_bytes += sdsl::write_member(this->flags, out, child, "flags");
  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
FastLocate::Header::load(std::istream& in)
{
  sdsl::read_member(this->tag, in);
  sdsl::read_member(this->version, in);
  sdsl::read_member(this->max_length, in);
  sdsl::read_member(this->flags, in);
}

void
FastLocate::Header::check() const
{
  if(this->tag != TAG)
  {
    throw sdsl::simple_sds::InvalidData("FastLocate: Invalid tag");
  }

  if(this->version != VERSION)
  {
    std::string msg = "FastLocate: Expected v" + std::to_string(VERSION) + ", got v" + std::to_string(this->version);
    throw sdsl::simple_sds::InvalidData(msg);
  }

  std::uint64_t mask = 0;
  switch(this->version)
  {
  case VERSION:
    mask = 0; break;
  }
  if((this->flags & mask) != this->flags)
  {
    throw sdsl::simple_sds::InvalidData("FastLocate: Invalid flags");
  }
}

//------------------------------------------------------------------------------

FastLocate::FastLocate() :
  index(nullptr)
{
}

FastLocate::FastLocate(const FastLocate& source)
{
  this->copy(source);
}

FastLocate::FastLocate(FastLocate&& source)
{
  *this = std::move(source);
}

FastLocate::~FastLocate()
{
}

void
FastLocate::swap(FastLocate& another)
{
  if(this != &another)
  {
    std::swap(this->index, another.index);
    std::swap(this->header, another.header);
    this->samples.swap(another.samples);
    this->last.swap(another.last);
    this->last_to_run.swap(another.last_to_run);
    this->comp_to_run.swap(another.comp_to_run);
  }
}

FastLocate&
FastLocate::operator=(const FastLocate& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

FastLocate&
FastLocate::operator=(FastLocate&& source)
{
  if(this != &source)
  {
    this->index = source.index;
    this->header = std::move(source.header);
    this->samples = std::move(source.samples);
    this->last = std::move(source.last);
    this->last_to_run = std::move(source.last_to_run);
    this->comp_to_run = std::move(source.comp_to_run);
  }
  return *this;
}

size_type
FastLocate::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->header.serialize(out, child, "header");
  written_bytes += this->samples.serialize(out, child, "samples");
  written_bytes += this->last.serialize(out, child, "last");
  written_bytes += this->last_to_run.serialize(out, child, "last_to_run");
  written_bytes += this->comp_to_run.serialize(out, child, "comp_to_run");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
FastLocate::load(std::istream& in)
{
  this->header.load(in);
  this->header.check();
  this->header.setVersion(); // Update to the current version.

  this->samples.load(in);
  this->last.load(in);
  this->last_to_run.load(in);
  this->comp_to_run.load(in);
}

void
FastLocate::copy(const FastLocate& source)
{
  this->header = source.header;
  this->samples = source.samples;
  this->last = source.last;
  this->last_to_run = source.last_to_run;
  this->comp_to_run = source.comp_to_run;
}

//------------------------------------------------------------------------------

FastLocate::FastLocate(const GBWT& source) :
  index(&source)
{
  double start = readTimer();

  if(this->index->empty())
  {
    if(Verbosity::level >= Verbosity::FULL)
    {
      std::cerr << "FastLocate::FastLocate(): The input GBWT is empty" << std::endl;
    }
    return;
  }

  // Determine the number of logical runs before each record.
  size_type total_runs = 0;
  this->comp_to_run.resize(this->index->effective());
  this->index->bwt.forEach([&](size_type comp, const CompressedRecord& record)
  {
    this->comp_to_run[comp] = total_runs; total_runs += record.runs().second;
  });
  sdsl::util::bit_compress(this->comp_to_run);
  if(Verbosity::level >= Verbosity::FULL)
  {
    std::cerr << "FastLocate::FastLocate(): " << total_runs << " logical runs in the GBWT" << std::endl;
  }

  // Global sample buffers.
  struct sample_record
  {
    size_type seq_id, seq_offset, run_id;

    // Sort by text position.
    bool operator<(const sample_record& another) const
    {
      return (this->seq_id < another.seq_id || (this->seq_id == another.seq_id && this->seq_offset < another.seq_offset));
    }
  };
  std::vector<sample_record> head_samples, tail_samples;
  head_samples.reserve(total_runs);
  tail_samples.reserve(total_runs);

  // Run identifier for each offset in the endmarker. We cannot get this
  // information efficiently with random access.
  if(Verbosity::level >= Verbosity::FULL)
  {
    std::cerr << "FastLocate::FastLocate(): Processing the endmarker record" << std::endl;
  }
  std::vector<size_type> endmarker_runs(this->index->sequences(), 0);
  {
    size_type run_id = 0;
    edge_type prev = this->index->start(0);
    for(size_type i = 1; i < this->index->sequences(); i++)
    {
      edge_type curr = this->index->start(i);
      if(curr.first == ENDMARKER || curr.first != prev.first) { run_id++; prev = curr; }
      endmarker_runs[i] = run_id;
    }
  }

  // Extract the samples from each sequence.
  double extract_start = readTimer();
  if(Verbosity::level >= Verbosity::FULL)
  {
    std::cerr << "FastLocate::FastLocate(): Extracting head/tail samples" << std::endl;
  }
  #pragma omp parallel for schedule(dynamic, 1)
  for(size_type i = 0; i < this->index->sequences(); i++)
  {
    std::vector<sample_record> head_buffer, tail_buffer;
    size_type seq_offset = 0, run_id = endmarker_runs[i];
    if(i == 0 || run_id != endmarker_runs[i - 1])
    {
      head_buffer.push_back({ i, seq_offset, this->globalRunId(ENDMARKER, run_id) });
    }
    if(i + 1 >= this->index->sequences() || run_id != endmarker_runs[i + 1])
    {
      tail_buffer.push_back({ i, seq_offset, this->globalRunId(ENDMARKER, run_id) });
    }
    edge_type curr = this->index->start(i); seq_offset++;
    range_type run(0, 0);
    while(curr.first != ENDMARKER)
    {
      edge_type next = this->index->record(curr.first).LF(curr.second, run, run_id);
      if(curr.second == run.first)
      {
        head_buffer.push_back({ i, seq_offset, this->globalRunId(curr.first, run_id) });
      }
      if(curr.second == run.second)
      {
        tail_buffer.push_back({ i, seq_offset, this->globalRunId(curr.first, run_id) });
      }
      curr = next; seq_offset++;
    }
    // GBWT is an FM-index of the reverse paths. The sequence offset r-index needs
    // is the distance to the BWT position with the endmarker (to the end of the
    // path, to the start of the string).
    for(sample_record& record : head_buffer) { record.seq_offset = seq_offset - 1 - record.seq_offset; }
    for(sample_record& record : tail_buffer) { record.seq_offset = seq_offset - 1 - record.seq_offset; }
    #pragma omp critical
    {
      this->header.max_length = std::max(this->header.max_length, seq_offset);
      head_samples.insert(head_samples.end(), head_buffer.begin(), head_buffer.end());
      tail_samples.insert(tail_samples.end(), tail_buffer.begin(), tail_buffer.end());
    }
  }
  sdsl::util::clear(endmarker_runs);
  if(Verbosity::level >= Verbosity::BASIC)
  {
    double seconds = readTimer() - extract_start;
    std::cerr << "FastLocate::FastLocate(): Extracted " << head_samples.size() << " / " << tail_samples.size() << " head/tail samples in " << seconds << " seconds" << std::endl;
  }

  // Store the head samples.
  if(Verbosity::level >= Verbosity::FULL)
  {
    std::cerr << "FastLocate::FastLocate(): Storing the head samples" << std::endl;
  }
  parallelQuickSort(head_samples.begin(), head_samples.end(), [](const sample_record& a, const sample_record& b)
  {
    return (a.run_id < b.run_id);
  });
  this->samples.width(sdsl::bits::length(this->pack(this->index->sequences() - 1, this->header.max_length - 1)));
  this->samples.resize(total_runs);
  for(size_type i = 0; i < total_runs; i++)
  {
    this->samples[i] = this->pack(head_samples[i].seq_id, head_samples[i].seq_offset);
  }
  sdsl::util::clear(head_samples);

  // Store the tail samples.
  if(Verbosity::level >= Verbosity::FULL)
  {
    std::cerr << "FastLocate::FastLocate(): Storing the tail samples" << std::endl;
  }
  parallelQuickSort(tail_samples.begin(), tail_samples.end());
  sdsl::sd_vector_builder builder(this->index->sequences() * this->header.max_length, total_runs);
  this->last_to_run.width(sdsl::bits::length(total_runs - 1));
  this->last_to_run.resize(total_runs);
  for(size_type i = 0; i < total_runs; i++)
  {
    builder.set_unsafe(this->pack(tail_samples[i].seq_id, tail_samples[i].seq_offset));
    this->last_to_run[i] = tail_samples[i].run_id;
  }
  sdsl::util::clear(tail_samples);
  this->last = sdsl::sd_vector<>(builder);

  if(Verbosity::level >= Verbosity::BASIC)
  {
    double seconds = readTimer() - start;
    std::cerr << "FastLocate::FastLocate(): Processed " << this->index->sequences() << " sequences of total length " << this->index->size() << " in " << seconds << " seconds" << std::endl;
  }
}

//------------------------------------------------------------------------------

SearchState
FastLocate::find(node_type node, size_type& first) const
{
  if(!(this->index->contains(node))) { return SearchState(); }

  CompressedRecord record = this->index->record(node);
  if(!(record.empty()))
  {
    first = this->getSample(node, 0);
  }
  return SearchState(node, 0, record.size() - 1);
}

SearchState
FastLocate::extend(SearchState state, node_type node, size_type& first) const
{
  if(state.empty() || !(this->index->contains(node))) { return SearchState(); }

  CompressedRecord record = this->index->record(state.node);
  bool starts_with_node = false;
  size_type run_id = invalid_offset();
  state.range = record.LF(state.range, node, starts_with_node, run_id);
  if(!(state.empty()))
  {
    // The position at the start of the resulting range is the successor of the
    // first occurrence of the query node in the query range. We decrement the
    // offset instead of incrementing it, because sequence offset is the distance
    // to the end of the sequence.
    if(starts_with_node) { first--; }
    else
    {
      first = this->getSample(state.node, run_id) - 1;
    }
  }
  state.node = node;

  return state;
}

std::vector<size_type>
FastLocate::locate(SearchState state, size_type first) const
{
  std::vector<size_type> result;
  if(!(this->index->contains(state))) { return result; }
  result.reserve(state.size());

  // Find the nearest run start and use the sample there as the first hit,
  // if the caller did not provide it.
  size_type offset_of_first = state.range.first;
  if(first == NO_POSITION)
  {
    CompressedRecord record = this->index->record(state.node);
    CompressedRecordIterator iter(record);
    size_type run_id = 0;
    while(!(iter.end()) && iter.offset() <= state.range.first)
    {
      ++iter; run_id++;
    }
    first = this->getSample(state.node, run_id);
    offset_of_first = iter.offset() - iter->second;
  }

  // Iterate until the start of the range and locate the first occurrence.
  while(offset_of_first < state.range.first)
  {
    first = this->locateNext(first);
    offset_of_first++;
  }
  result.push_back(this->seqId(first));

  // Locate the remaining occurrences.
  for(size_type i = state.range.first + 1; i <= state.range.second; i++)
  {
    first = this->locateNext(first);
    result.push_back(this->seqId(first));
  }

  removeDuplicates(result, false);
  return result;
}

std::vector<size_type>
FastLocate::decompressSA(node_type node) const
{
  std::vector<size_type> result;
  SearchState state = this->index->find(node);
  if(state.empty()) { return result; }

  result.reserve(state.size());
  result.push_back(this->locateFirst(node));
  for(size_type i = state.range.first; i < state.range.second; i++)
  {
    result.push_back(this->locateNext(result.back()));
  }

  return result;
}

std::vector<size_type>
FastLocate::decompressDA(node_type node) const
{
  std::vector<size_type> result = this->decompressSA(node);
  for(size_type i = 0; i < result.size(); i++)
  {
    result[i] = this->seqId(result[i]);
  }
  return result;
}

//------------------------------------------------------------------------------

size_type
FastLocate::locateNext(size_type prev) const
{
  auto iter = this->last.predecessor(prev);
  return this->samples[this->last_to_run[iter->first] + 1] + (prev - iter->second);
}

//------------------------------------------------------------------------------

void
printStatistics(const FastLocate& index, const std::string& name)
{
  printHeader(indexType(index)); std::cout << name << std::endl;
  printHeader("Runs"); std::cout << index.size() << std::endl;
  printHeader("Size"); std::cout << inMegabytes(sdsl::size_in_bytes(index)) << " MB" << std::endl;
  std::cout << std::endl;
}

std::string
indexType(const FastLocate&)
{
  return "R-index";
}

//------------------------------------------------------------------------------

} // namespace gbwt
