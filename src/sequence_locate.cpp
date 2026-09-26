#include <gbwt/sequence_locate.h>
#include <gbwt/internal.h>

#include <algorithm>

namespace gbwt
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr std::uint32_t SequenceLocate::Header::TAG;
constexpr std::uint32_t SequenceLocate::Header::VERSION;

constexpr size_type SequenceLocate::NO_POSITION;

//------------------------------------------------------------------------------

// Other class variables.

const std::string SequenceLocate::EXTENSION = ".sri";

//------------------------------------------------------------------------------

SequenceLocate::Header::Header() :
  tag(TAG), version(VERSION),
  max_length(1),
  flags(0)
{
}

size_type
SequenceLocate::Header::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
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
SequenceLocate::Header::load(std::istream& in)
{
  sdsl::read_member(this->tag, in);
  sdsl::read_member(this->version, in);
  sdsl::read_member(this->max_length, in);
  sdsl::read_member(this->flags, in);
}

void
SequenceLocate::Header::check() const
{
  if(this->tag != TAG)
  {
    throw sdsl::simple_sds::InvalidData("SequenceLocate: Invalid tag");
  }

  if(this->version != VERSION)
  {
    std::string msg = "SequenceLocate: Expected version " + std::to_string(VERSION) + ", got version " + std::to_string(this->version);
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
    throw sdsl::simple_sds::InvalidData("SequenceLocate: Invalid flags");
  }
}

//------------------------------------------------------------------------------

SequenceLocate::SequenceLocate() :
  index(nullptr)
{
}

SequenceLocate::SequenceLocate(const SequenceLocate& source)
{
  this->copy(source);
}

SequenceLocate::SequenceLocate(SequenceLocate&& source) noexcept
{
  *this = std::move(source);
}

SequenceLocate::~SequenceLocate()
{
}

void
SequenceLocate::swap(SequenceLocate& another) noexcept
{
  if(this != &another)
  {
    std::swap(this->index, another.index);
    this->length.swap(another.length);
    std::swap(this->header, another.header);
    this->samples.swap(another.samples);
    this->last.swap(another.last);
    this->last_to_run.swap(another.last_to_run);
    this->comp_to_run.swap(another.comp_to_run);
    this->sequence_length.swap(another.sequence_length);
  }
}

SequenceLocate&
SequenceLocate::operator=(const SequenceLocate& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

SequenceLocate&
SequenceLocate::operator=(SequenceLocate&& source) noexcept
{
  if(this != &source)
  {
    this->index = source.index;
    this->length = std::move(source.length);
    this->header = std::move(source.header);
    this->samples = std::move(source.samples);
    this->last = std::move(source.last);
    this->last_to_run = std::move(source.last_to_run);
    this->comp_to_run = std::move(source.comp_to_run);
    this->sequence_length = std::move(source.sequence_length);
  }
  return *this;
}

size_type
SequenceLocate::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->header.serialize(out, child, "header");
  written_bytes += this->samples.serialize(out, child, "samples");
  written_bytes += this->last.serialize(out, child, "last");
  written_bytes += this->last_to_run.serialize(out, child, "last_to_run");
  written_bytes += this->comp_to_run.serialize(out, child, "comp_to_run");
  written_bytes += this->sequence_length.serialize(out, child, "sequence_length");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
SequenceLocate::load(std::istream& in)
{
  this->header.load(in);
  this->header.check();
  this->header.setVersion(); // Update to the current version.

  this->samples.load(in);
  this->last.load(in);
  this->last_to_run.load(in);
  this->comp_to_run.load(in);
  this->sequence_length.load(in);
}

void
SequenceLocate::copy(const SequenceLocate& source)
{
  this->index = source.index;
  this->length = source.length;
  this->header = source.header;
  this->samples = source.samples;
  this->last = source.last;
  this->last_to_run = source.last_to_run;
  this->comp_to_run = source.comp_to_run;
  this->sequence_length = source.sequence_length;
}

//------------------------------------------------------------------------------

/*
  Same as the FastLocate construction, but the sequence offset grows by the node
  length instead of by one. We sample the offset at the end of each node and flip
  it at the end, so that the offset is the bp distance from the end of the node to
  the end of the sequence.
*/

SequenceLocate::SequenceLocate(const GBWT& source, const length_function& node_length) :
  index(&source), length(node_length)
{
  double start = readTimer();

  if(this->index->empty())
  {
    if(Verbosity::level >= Verbosity::FULL)
    {
      std::cerr << "SequenceLocate::SequenceLocate(): The input GBWT is empty" << std::endl;
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
    std::cerr << "SequenceLocate::SequenceLocate(): " << total_runs << " logical runs in the GBWT" << std::endl;
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
    std::cerr << "SequenceLocate::SequenceLocate(): Processing the endmarker record" << std::endl;
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
    std::cerr << "SequenceLocate::SequenceLocate(): Extracting head/tail samples" << std::endl;
  }
  std::vector<size_type> lengths(this->index->sequences(), 0);
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
    edge_type curr = this->index->start(i);
    if(curr.first != ENDMARKER) { seq_offset += this->length(curr.first); }
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
      curr = next;
      if(curr.first != ENDMARKER) { seq_offset += this->length(curr.first); }
    }
    // Flip the offsets to make them relative to the end of the sequence.
    for(sample_record& record : head_buffer) { record.seq_offset = seq_offset - record.seq_offset; }
    for(sample_record& record : tail_buffer) { record.seq_offset = seq_offset - record.seq_offset; }
    lengths[i] = seq_offset;
    #pragma omp critical
    {
      this->header.max_length = std::max(this->header.max_length, seq_offset + 1);
      head_samples.insert(head_samples.end(), head_buffer.begin(), head_buffer.end());
      tail_samples.insert(tail_samples.end(), tail_buffer.begin(), tail_buffer.end());
    }
  }
  sdsl::util::clear(endmarker_runs);
  this->sequence_length.width(sdsl::bits::length(this->header.max_length - 1));
  this->sequence_length.resize(lengths.size());
  for(size_type i = 0; i < lengths.size(); i++) { this->sequence_length[i] = lengths[i]; }
  sdsl::util::clear(lengths);
  if(Verbosity::level >= Verbosity::BASIC)
  {
    double seconds = readTimer() - extract_start;
    std::cerr << "SequenceLocate::SequenceLocate(): Extracted " << head_samples.size() << " / " << tail_samples.size() << " head/tail samples in " << seconds << " seconds" << std::endl;
  }

  // Store the head samples.
  if(Verbosity::level >= Verbosity::FULL)
  {
    std::cerr << "SequenceLocate::SequenceLocate(): Storing the head samples" << std::endl;
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
    std::cerr << "SequenceLocate::SequenceLocate(): Storing the tail samples" << std::endl;
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
    std::cerr << "SequenceLocate::SequenceLocate(): Processed " << this->index->sequences() << " sequences of total length " << this->index->size() << " in " << seconds << " seconds" << std::endl;
  }
}

//------------------------------------------------------------------------------

SearchState
SequenceLocate::find(node_type node, size_type& first) const
{
  if(!(this->index->contains(node))) { return SearchState(); }

  CompressedRecord record = this->index->record(node);
  if(!(record.empty()))
  {
    first = this->getSample(node, 0);
  }
  return SearchState(node, 0, record.size() - 1);
}

/*
  The first occurrence in the new range is the successor of the first occurrence
  of `node` in the old range. The offset decreases by the length of `node`, because
  the offset is the distance from the end of the node to the end of the sequence.
*/

SearchState
SequenceLocate::extend(SearchState state, node_type node, size_type& first) const
{
  if(state.empty() || !(this->index->contains(node))) { return SearchState(); }

  CompressedRecord record = this->index->record(state.node);
  bool starts_with_node = false;
  size_type run_id = invalid_offset();
  state.range = record.LF(state.range, node, starts_with_node, run_id);
  if(!(state.empty()))
  {
    if(starts_with_node) { first -= this->length(node); }
    else
    {
      first = this->getSample(state.node, run_id) - this->length(node);
    }
  }
  state.node = node;

  return state;
}

/*
  If the caller did not provide the first occurrence, we start from the sample at
  the head of the run containing the start of the range and iterate with locateNext().
*/

size_type
SequenceLocate::firstPosition(SearchState state, size_type first) const
{
  size_type offset_of_first = state.range.first;
  if(first == NO_POSITION)
  {
    CompressedRecord record = this->index->record(state.node);
    CompressedRecordIterator iter(record);
    while(!(iter.end()) && iter.offset() <= state.range.first) { ++iter; }
    // Each endmarker is a separate logical run, and runId() is that of the last one.
    size_type run_id = iter.runId();
    if(record.successor(iter->first) == ENDMARKER) { run_id -= iter->second - 1; }
    first = this->getSample(state.node, run_id);
    offset_of_first = iter.offset() - iter->second;
  }

  while(offset_of_first < state.range.first)
  {
    first = this->locateNext(first);
    offset_of_first++;
  }
  return first;
}

std::vector<size_type>
SequenceLocate::locate(SearchState state, size_type first) const
{
  std::vector<size_type> result;
  if(!(this->index->contains(state))) { return result; }
  result.reserve(state.size());

  first = this->firstPosition(state, first);
  result.push_back(this->seqId(first));
  for(size_type i = state.range.first + 1; i <= state.range.second; i++)
  {
    first = this->locateNext(first);
    result.push_back(this->seqId(first));
  }

  removeDuplicates(result, false);
  return result;
}

// Returns (sequence id, sequence offset) for each occurrence in the range.
std::vector<std::pair<size_type, size_type>>
SequenceLocate::locatePositions(SearchState state, size_type first) const
{
  std::vector<std::pair<size_type, size_type>> result;
  if(!(this->index->contains(state))) { return result; }
  result.reserve(state.size());

  first = this->firstPosition(state, first);
  result.push_back(this->unpack(first));
  for(size_type i = state.range.first + 1; i <= state.range.second; i++)
  {
    first = this->locateNext(first);
    result.push_back(this->unpack(first));
  }

  return result;
}

// Converts the positions from locatePositions() to forward starts of the node.
std::vector<std::pair<size_type, size_type>>
SequenceLocate::locateForward(SearchState state, size_type first) const
{
  std::vector<std::pair<size_type, size_type>> result = this->locatePositions(state, first);
  if(state.node == ENDMARKER) { return result; }
  for(auto& pos : result) { pos.second = this->forwardStart(pos.first, pos.second, state.node); }
  return result;
}

//------------------------------------------------------------------------------

/*
  The target offset is `len - forward_bp`, and the node containing it is the first
  visit with offset < target. Tail samples at or before the target have offset >=
  target, so we take the successor in `last`. From there, each LF step moves to the
  next node and reduces the offset by its length.
*/

edge_type
SequenceLocate::locateSequence(size_type seq_id, size_type forward_bp) const
{
  edge_type gbwt_pos;
  return this->locateSequence(seq_id, forward_bp, gbwt_pos);
}

edge_type
SequenceLocate::locateSequence(size_type seq_id, size_type forward_bp, edge_type& gbwt_pos) const
{
  gbwt_pos = invalid_edge();
  if(seq_id >= this->index->sequences() || forward_bp >= this->sequenceLength(seq_id)) { return invalid_edge(); }

  size_type total = this->sequenceLength(seq_id);
  size_type target = total - forward_bp;

  // Default to the endmarker if there is no usable tail sample.
  edge_type pos(ENDMARKER, seq_id);
  size_type offset_to_end = total;
  auto iter = this->last.successor(this->pack(seq_id, target));
  if(iter != this->last.one_end() && this->seqId(iter->second) == seq_id)
  {
    edge_type tail = this->runTail(this->last_to_run[iter->first]);
    if(tail.first != ENDMARKER)
    {
      pos = tail;
      offset_to_end = this->seqOffset(iter->second);
    }
  }

  while(offset_to_end >= target)
  {
    pos = this->index->LF(pos);
    offset_to_end -= this->length(pos.first);
  }

  gbwt_pos = pos;
  return edge_type(pos.first, forward_bp - (total - offset_to_end - this->length(pos.first)));
}

// Returns the GBWT position at the end of the given global run.
edge_type
SequenceLocate::runTail(size_type run) const
{
  // The last record whose first run is at or before `run`.
  size_type comp = std::upper_bound(this->comp_to_run.begin(), this->comp_to_run.end(), run) - this->comp_to_run.begin() - 1;
  node_type node = this->index->toNode(comp);
  if(node == ENDMARKER) { return edge_type(ENDMARKER, 0); }

  size_type local_run = run - this->comp_to_run[comp];
  CompressedRecord record = this->index->record(node);
  CompressedRecordIterator iter(record);
  while(!(iter.end()) && iter.runId() < local_run) { ++iter; }
  // In a concrete run of endmarkers, runId() is that of the last occurrence.
  return edge_type(node, iter.offset() - 1 - (iter.runId() - local_run));
}

//------------------------------------------------------------------------------

std::vector<size_type>
SequenceLocate::decompressSA(node_type node) const
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
SequenceLocate::decompressDA(node_type node) const
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
SequenceLocate::locateNext(size_type prev) const
{
  auto iter = this->last.predecessor(prev);
  return this->samples[this->last_to_run[iter->first] + 1] + (prev - iter->second);
}

//------------------------------------------------------------------------------

void
printStatistics(const SequenceLocate& index, const std::string& name)
{
  printHeader(indexType(index)); std::cout << name << std::endl;
  printHeader("Runs"); std::cout << index.size() << std::endl;
  printHeader("Size"); std::cout << inMegabytes(sdsl::size_in_bytes(index)) << " MB" << std::endl;
  std::cout << std::endl;
}

std::string
indexType(const SequenceLocate&)
{
  return "Sequence r-index";
}

//------------------------------------------------------------------------------

} // namespace gbwt
