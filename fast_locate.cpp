/*
  Copyright (c) 2020 Jouni Siren

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

bool
FastLocate::Header::check() const
{
  return (this->tag == TAG && this->version == VERSION && this->flags == 0);
}

//------------------------------------------------------------------------------

FastLocate::FastLocate()
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
  if(!(this->header.check()))
  {
    std::cerr << "FastLocate::load(): Invalid or old header" << std::endl;
    std::cerr << "FastLocate::load(): Header version is " << this->header.version << "; expected " << Header::VERSION << std::endl;
    return;
  }
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

SearchState
FastLocate::find(node_type node, size_type& first) const
{
  if(!(this->index.contains(node))) { return SearchState(); }

  CompressedRecord record = this->index.record(node);
  if(!(record.empty()))
  {
    first = this->getSample(node, 0);
  }
  return SearchState(node, 0, record.size() - 1);
}

template<class Iterator>
SearchState
FastLocate::find(Iterator begin, Iterator end, size_type& first) const
{
  if(begin == end) { return SearchState(); }

  SearchState state = this->find(*begin, first);
  ++begin;

  return this->extend(state, begin, end, first);
}

SearchState
FastLocate::extend(SearchState state, node_type node, size_type& first) const
{
  if(state.empty() || !(this->index.contains(node))) { return SearchState(); }

  CompressedRecord record = this->index.record(node);
  bool starts_with_node = false;
  size_type run_id = invalid_offset();
  size_type first_occ = invalid_offset();
  state.range = record.LF(state, node, starts_with_node, run_id);
  state.node = node;
  if(!(state.empty()))
  {
    // The position at the start of the resulting range is the successor of the
    // first occurrence of the query node in the query range. We decrement the
    // offset instead of incrementing it, because sequence offset is the distance
    // to the end of the sequence.
    if(starts_with_node) { first--; }
    else
    {
      first = this->getSample(node, run_id) - 1;
    }
  }

  return state;
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

std::vector<size_type>
FastLocate::locate(SearchState state, size_type first) const
{
  std::vector<size_type> result;
  if(!(this->contains(state))) { return result; }
  result.reserve(state.size());

  // Find the nearest run start and use the sample there as the first hit,
  // if the caller did not provide it.
  size_type offset_of_first = state.range.first;
  if(first == NO_POSITION)
  {
    CompressedRecord record = this->index.record(state.node);
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

//------------------------------------------------------------------------------

size_type
FastLocate::locateNext(size_type prev) const
{
  SDIterator iter(this->last, prev, SDIterator::query_predecessor);
  return this->samples[this->last_to_run[iter.rank()] + 1] + (prev - *iter);
}

//------------------------------------------------------------------------------

} // namespace gbwt
