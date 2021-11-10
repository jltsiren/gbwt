/*
  Copyright (c) 2017, 2019, 2020, 2021 Jouni Siren
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

#include <gbwt/gbwt.h>

#include <gbwt/dynamic_gbwt.h>
#include <gbwt/internal.h>

namespace gbwt
{

//------------------------------------------------------------------------------

// Other class variables.

const std::string GBWT::EXTENSION = ".gbwt";

//------------------------------------------------------------------------------

GBWT::GBWT()
{
  this->addSource();
}

GBWT::GBWT(const GBWT& source)
{
  this->copy(source);
}

GBWT::GBWT(const DynamicGBWT& source) :
  header(source.header),
  tags(source.tags),
  bwt(source.bwt), da_samples(source.bwt),
  metadata(source.metadata)
{
  this->cacheEndmarker();
}

GBWT::GBWT(GBWT&& source)
{
  *this = std::move(source);
}

GBWT::~GBWT()
{
}

void
GBWT::swap(GBWT& another)
{
  if(this != &another)
  {
    this->header.swap(another.header);
    this->tags.swap(another.tags);
    this->bwt.swap(another.bwt);
    this->da_samples.swap(another.da_samples);
    this->metadata.swap(another.metadata);
    this->endmarker_record.swap(another.endmarker_record);
  }
}

GBWT&
GBWT::operator=(const GBWT& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

GBWT&
GBWT::operator=(const DynamicGBWT& source)
{
  // Clear the data to save memory.
  this->bwt = RecordArray();
  this->da_samples = DASamples();

  this->header = source.header;
  this->tags = source.tags;
  this->bwt = RecordArray(source.bwt);
  this->da_samples = DASamples(source.bwt);
  this->metadata = source.metadata;
  this->cacheEndmarker();

  return *this;
}

GBWT&
GBWT::operator=(GBWT&& source)
{
  if(this != &source)
  {
    this->header = std::move(source.header);
    this->tags = std::move(source.tags);
    this->bwt = std::move(source.bwt);
    this->da_samples = std::move(source.da_samples);
    this->metadata = std::move(source.metadata);
    this->endmarker_record = std::move(source.endmarker_record);
  }
  return *this;
}

void
GBWT::resample(size_type sample_interval)
{
  this->da_samples = DASamples(); // Delete samples to save memory.
  std::vector<std::pair<node_type, sample_type>> samples = gbwt::resample(*this, sample_interval);
  this->da_samples = DASamples(this->bwt, samples);
}

size_type
GBWT::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->header.serialize(out, child, "header");
  written_bytes += this->tags.serialize(out, child, "tags");
  written_bytes += this->bwt.serialize(out, child, "bwt");
  written_bytes += this->da_samples.serialize(out, child, "da_samples");

  if(this->hasMetadata())
  {
    written_bytes += this->metadata.serialize(out, child, "metadata");
  }

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
GBWT::load(std::istream& in)
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

  // Read the BWT.
  if(simple_sds) { this->bwt.simple_sds_load(in); }
  else { this->bwt.load(in); }
  if(this->bwt.size() != this->effective())
  {
    throw sdsl::simple_sds::InvalidData("GBWT: BWT record count / alphabet size mismatch");
  }

  // Cache the endmarker before the potential resampling.
  this->cacheEndmarker();

  // Read the DA samples.
  if(simple_sds)
  {
    // The samples may be absent or from an unknown source.
    bool found_samples = false;
    if(source_is_self) { found_samples = sdsl::simple_sds::load_option(this->da_samples, in); }
    else { sdsl::simple_sds::skip_option(in); }
    if(!found_samples) { this->resample(DynamicGBWT::SAMPLE_INTERVAL); }
  }
  else { this->da_samples.load(in); }
  if(this->da_samples.records() != this->effective())
  {
    throw sdsl::simple_sds::InvalidData("GBWT: Sample record count / alphabet size mismatch");
  }

  // Read the metadata.
  if(simple_sds)
  {
    bool loaded_metadata = sdsl::simple_sds::load_option(this->metadata, in);
    if(loaded_metadata != this->hasMetadata())
    {
      throw sdsl::simple_sds::InvalidData("GBWT: Invalid metadata flag in the header");
    }
  }
  else if(this->hasMetadata()) { this->metadata.load(in); }
  if(this->hasMetadata() && this->metadata.hasPathNames())
  {
    size_type expected_paths = (this->bidirectional() ? this->sequences() / 2 : this->sequences());
    if(this->metadata.paths() != expected_paths)
    {
      throw sdsl::simple_sds::InvalidData("GBWT: Path name / sequence count mismatch");
    }
  }
}

void
GBWT::simple_sds_serialize(std::ostream& out) const
{
  GBWTHeader h = this->header;
  h.set(GBWTHeader::FLAG_SIMPLE_SDS); // We only set this flag in the serialized header.
  sdsl::simple_sds::serialize_value(h, out);

  this->tags.simple_sds_serialize(out);
  this->bwt.simple_sds_serialize(out);
  sdsl::simple_sds::serialize_option(this->da_samples, out);
  if(this->hasMetadata()) { sdsl::simple_sds::serialize_option(this->metadata, out); }
  else { sdsl::simple_sds::empty_option(out); }
}

void
GBWT::simple_sds_load(std::istream& in)
{
  // The same load() function can handle the SDSL and simple-sds formats.
  this->load(in);
}

size_t
GBWT::simple_sds_size() const
{
  size_t result = sdsl::simple_sds::value_size(this->header);
  result += this->tags.simple_sds_size();
  result += this->bwt.simple_sds_size();
  result += sdsl::simple_sds::option_size(this->da_samples);
  if(this->hasMetadata()) { result += sdsl::simple_sds::option_size(this->metadata); }
  else { result += sdsl::simple_sds::empty_option_size(); }
  return result;
}

void
GBWT::copy(const GBWT& source)
{
  this->header = source.header;
  this->tags = source.tags;
  this->bwt = source.bwt;
  this->da_samples = source.da_samples;
  this->metadata = source.metadata;
  this->endmarker_record = source.endmarker_record;
}

void
GBWT::resetTags()
{
  this->tags.clear();
  this->addSource();
}

void
GBWT::addSource()
{
  this->tags.set(Version::SOURCE_KEY, Version::SOURCE_VALUE);
}

//------------------------------------------------------------------------------

GBWT::GBWT(const std::vector<GBWT>& sources)
{
  if(sources.empty()) { return; }

  // We cannot know what to do with the tags in the sources.
  this->addSource();

  // Merge the headers.
  size_type valid_sources = 0;
  bool is_bidirectional = true;
  for(const GBWT& source : sources)
  {
    if(source.empty()) { continue; }
    this->header.sequences += source.header.sequences;
    this->header.size += source.header.size;
    if(valid_sources == 0)
    {
      this->header.offset = source.header.offset;
      this->header.alphabet_size = source.header.alphabet_size;
    }
    else
    {
      this->header.offset = std::min(this->header.offset, source.header.offset);
      this->header.alphabet_size = std::max(this->header.alphabet_size, source.header.alphabet_size);
    }
    if(!(source.bidirectional())) { is_bidirectional = false; }
    valid_sources++;
  }
  if(valid_sources == 0) { return; }
  if(is_bidirectional) { this->header.set(GBWTHeader::FLAG_BIDIRECTIONAL); }
  else { this->header.unset(GBWTHeader::FLAG_BIDIRECTIONAL); }

  // Determine the mapping between source comp values and merged comp values.
  std::vector<size_type> record_offsets(sources.size());
  for(size_type i = 0; i < sources.size(); i++)
  {
    record_offsets[i] = sources[i].header.offset - this->header.offset;
  }

  // Determine the origin of each record.
  sdsl::int_vector<0> origins(this->effective(), sources.size(), sdsl::bits::length(sources.size()));
  for(size_type source_id = 0; source_id < sources.size(); source_id++)
  {
    const GBWT& source = sources[source_id];
    source.bwt.forEach([&](size_type source_comp, const CompressedRecord& record)
    {
      if(source_comp == 0 || record.empty()) { return; }
      comp_type merged_comp = source_comp + record_offsets[source_id];
      if(origins[merged_comp] != sources.size())
      {
        std::cerr << "GBWT::GBWT(): Sources " << origins[merged_comp] << " and " << source_id << " both have node " << this->toNode(merged_comp) << std::endl;
        std::exit(EXIT_FAILURE);
      }
      origins[merged_comp] = source_id;
    });
  }

  // Interleave the BWTs.
  {
    std::vector<RecordArray const*> bwt_sources(sources.size());
    for(size_type i = 0; i < sources.size(); i++)
    {
      bwt_sources[i] = &(sources[i].bwt);
    }
    this->bwt = RecordArray(bwt_sources, origins);
  }

  // Interleave the samples.
  {
    std::vector<DASamples const*> sample_sources(sources.size());
    std::vector<size_type> sequence_counts(sources.size());
    for(size_type i = 0; i < sources.size(); i++)
    {
      sample_sources[i] = &(sources[i].da_samples);
      sequence_counts[i] = sources[i].sequences();
    }
    this->da_samples = DASamples(sample_sources, origins, record_offsets, sequence_counts);
  }

  // Merge the metadata from non-empty sources.
  {
    bool has_metadata = false, all_metadata = true;
    for(size_type i = 0; i < sources.size(); i++)
    {
      if(sources[i].empty()) { continue; }
      has_metadata |= sources[i].hasMetadata();
      all_metadata &= sources[i].hasMetadata();
    }
    if(has_metadata)
    {
      if(all_metadata)
      {
        std::vector<const Metadata*> source_metadata;
        source_metadata.reserve(sources.size());
        for(size_type i = 0; i < sources.size(); i++)
        {
          if(sources[i].empty()) { continue; }
          source_metadata.push_back(&(sources[i].metadata));
        }
        this->metadata = Metadata(source_metadata, true, false); // Same samples, different contigs.
        this->addMetadata();
      }
      else if(Verbosity::level >= Verbosity::BASIC)
      {
        std::cerr << "GBWT::merge(): All inputs do not have metadata" << std::endl;
      }
    }
  }

  this->cacheEndmarker();
}

//------------------------------------------------------------------------------

std::pair<size_type, size_type>
GBWT::runs() const
{
  std::pair<size_type, size_type> result(0, 0);
  this->bwt.forEach([&result](size_type, const CompressedRecord& record)
  {
    std::pair<size_type, size_type> temp = record.runs();
    result.first += temp.first; result.second += temp.second;
  });
  return result;
}

//------------------------------------------------------------------------------

std::vector<size_type>
GBWT::locate(SearchState state) const
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
    CompressedRecord current;
    sample_type sample;
    edge_type LF_result;
    range_type LF_range;

    for(size_type i = 0; i < positions.size(); i++)
    {
      if(positions[i].first != curr)              // Node changed.
      {
        curr = positions[i].first; current = this->record(curr);
        sample = this->da_samples.nextSample(this->toComp(curr), positions[i].second);
        LF_range.first = positions[i].second;
        LF_result = current.runLF(positions[i].second, LF_range.second);
      }
      if(sample.first < positions[i].second)      // Went past the sample.
      {
        sample = this->da_samples.nextSample(this->toComp(curr), positions[i].second);
      }
      if(sample.first > positions[i].second)      // Not sampled, also valid for invalid_sample().
      {
        if(positions[i].second > LF_range.second) // Went past the existing LF() result.
        {
          LF_range.first = positions[i].second;
          LF_result = current.runLF(positions[i].second, LF_range.second);
        }
        positions[tail] = edge_type(LF_result.first, LF_result.second + positions[i].second - LF_range.first);
        tail++;
      }
      else                                        // Found a sample.
      {
        result.push_back(sample.second);
      }
    }
    positions.resize(tail);
    sequentialSort(positions.begin(), positions.end());
  }

  removeDuplicates(result, false);
  return result;
}

//------------------------------------------------------------------------------

edge_type
GBWT::inverseLF(node_type from, size_type i) const
{
  if(!(this->bidirectional()) || from == ENDMARKER) { return invalid_edge(); }

  // Find the predecessor node id.
  CompressedRecord reverse_record = this->record(Node::reverse(from));
  node_type predecessor = reverse_record.predecessorAt(i);
  if(predecessor == invalid_node()) { return invalid_edge(); }

  // Determine the offset.
  CompressedRecord pred_record = this->record(predecessor);
  size_type offset = pred_record.offsetTo(from, i);
  if(offset == invalid_offset()) { return invalid_edge(); }

  return edge_type(predecessor, offset);
}

//------------------------------------------------------------------------------

CompressedRecord
GBWT::record(node_type node) const
{
  comp_type comp = this->toComp(node);
  std::pair<size_type, size_type> range = this->bwt.getRange(comp);
  return CompressedRecord(this->bwt.data, range.first, range.second);
}

void
GBWT::cacheEndmarker()
{
  if(this->empty()) { return; }
  this->endmarker_record = DecompressedRecord(this->record(ENDMARKER));
}

//------------------------------------------------------------------------------

void
printStatistics(const GBWT& gbwt, const std::string& name, std::ostream& out)
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
  printHeader("BWT", out) << inMegabytes(sdsl::size_in_bytes(gbwt.bwt)) << " MB" << std::endl;
  printHeader("DA samples", out) << inMegabytes(sdsl::size_in_bytes(gbwt.da_samples)) << " MB" << std::endl;
  printHeader("Total", out) << inMegabytes(sdsl::size_in_bytes(gbwt)) << " MB" << std::endl;
  if(gbwt.hasMetadata())
  {
    printHeader("Metadata", out) << gbwt.metadata << std::endl;
  }
  out << std::endl;
}

std::string
indexType(const GBWT&)
{
  return "Compressed GBWT";
}

//------------------------------------------------------------------------------

} // namespace gbwt
