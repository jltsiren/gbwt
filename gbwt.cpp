/*
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
#include <gbwt/internal.h>

namespace gbwt
{

//------------------------------------------------------------------------------

const std::string GBWT::EXTENSION = ".gbwt";

GBWT::GBWT()
{
}

GBWT::GBWT(const GBWT& source)
{
  this->copy(source);
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
    this->bwt.swap(another.bwt);
    this->da_samples.swap(another.da_samples);
  }
}

GBWT&
GBWT::operator=(const GBWT& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

GBWT&
GBWT::operator=(GBWT&& source)
{
  if(this != &source)
  {
    this->header = std::move(source.header);
    this->bwt = std::move(source.bwt);
    this->da_samples = std::move(source.da_samples);
  }
  return *this;
}

size_type
GBWT::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->header.serialize(out, child, "header");
  written_bytes += this->bwt.serialize(out, child, "bwt");
  written_bytes += this->da_samples.serialize(out, child, "da_samples");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
GBWT::load(std::istream& in)
{
  this->header.load(in);
  if(!(this->header.check()))
  {
    std::cerr << "GBWT::load(): Invalid header: " << this->header << std::endl;
  }

  this->bwt.load(in);
  this->da_samples.load(in);
}

void
GBWT::copy(const GBWT& source)
{
  this->header = source.header;
  this->bwt = source.bwt;
  this->da_samples = source.da_samples;
}

//------------------------------------------------------------------------------

size_type
GBWT::runs() const
{
  size_type start = 0, result = 0;
  for(comp_type comp = 0; comp < this->effective(); comp++)
  {
    size_type limit = this->bwt.limit(comp);
    CompressedRecord record(this->bwt.data, start, limit);
    result += record.runs();
    start = limit;
  }
  return result;
}

//------------------------------------------------------------------------------

size_type
GBWT::locate(node_type node, size_type i) const
{
  if(!(this->contains(node))) { return invalid_sequence(); }

  while(true)
  {
    size_type result = this->tryLocate(node, i);
    if(result != invalid_sequence()) { return result; }
    std::tie(node, i) = this->LF(node, i);
  }
}

std::vector<size_type>
GBWT::locate(node_type node, range_type range) const
{
  std::vector<size_type> result;
  if(!(this->contains(node)) || Range::empty(range)) { return result; }
  if(range.second >= this->count(node)) { return result; }

  // Initialize BWT positions for each offset in the range.
  std::vector<edge_type> positions(Range::length(range));
  for(size_type i = range.first; i <= range.second; i++)
  {
    positions[i - range.first] = edge_type(node, i);
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

CompressedRecord
GBWT::record(node_type node) const
{
  comp_type comp = this->toComp(node);
  size_type start = this->bwt.start(comp), limit = this->bwt.limit(comp);
  return CompressedRecord(this->bwt.data, start, limit);
}

//------------------------------------------------------------------------------

void
printStatistics(const GBWT& gbwt, const std::string& name)
{
  printHeader("Compressed GBWT"); std::cout << name << std::endl;
  printHeader("Total length"); std::cout << gbwt.size() << std::endl;
  printHeader("Sequences"); std::cout << gbwt.sequences() << std::endl;
  printHeader("Alphabet size"); std::cout << gbwt.sigma() << std::endl;
  printHeader("Effective"); std::cout << gbwt.effective() << std::endl;
  printHeader("Runs"); std::cout << gbwt.runs() << std::endl;
  printHeader("Samples"); std::cout << gbwt.samples() << std::endl;
  printHeader("BWT"); std::cout << inMegabytes(sdsl::size_in_bytes(gbwt.bwt)) << " MB" << std::endl;
  printHeader("Samples"); std::cout << inMegabytes(sdsl::size_in_bytes(gbwt.da_samples)) << " MB" << std::endl;
  printHeader("Total"); std::cout << inMegabytes(sdsl::size_in_bytes(gbwt)) << " MB" << std::endl;
  std::cout << std::endl;
}

//------------------------------------------------------------------------------

} // namespace gbwt
