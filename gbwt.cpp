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

#include "gbwt.h"
#include "internal.h"

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
    this->record_index.swap(another.record_index);
    sdsl::util::swap_support(this->record_select, another.record_select, &(this->record_index), &(another.record_index));
    this->bwt.swap(another.bwt);
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
    this->record_index = std::move(source.record_index);
    this->record_select = std::move(source.record_select);
    this->bwt = std::move(source.bwt);
  }
  return *this;
}

size_type
GBWT::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->header.serialize(out, child, "header");
  written_bytes += this->record_index.serialize(out, child, "record_index");
  written_bytes += this->record_select.serialize(out, child, "record_select");

  // Serialize the BWT.
  size_type bwt_bytes = this->bwt.size() * sizeof(byte_type);
  sdsl::structure_tree_node* bwt_node =
    sdsl::structure_tree::add_child(child, "bwt", "std::vector<gbwt::byte_type>");
  out.write((const char*)(this->bwt.data()), bwt_bytes);
  sdsl::structure_tree::add_size(bwt_node, bwt_bytes);
  written_bytes += bwt_bytes;

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
GBWT::load(std::istream& in)
{
  // Read the header.
  this->header.load(in);
  if(!(this->header.check()))
  {
    std::cerr << "GBWT::load(): Invalid header: " << this->header << std::endl;
  }

  // Read the record index.
  this->record_index.load(in);
  this->record_select.load(in, &(this->record_index));

  // Read the BWT.
  this->bwt.resize(this->record_index.size());
  in.read((char*)(this->bwt.data()), this->bwt.size() * sizeof(byte_type));
}

void
GBWT::copy(const GBWT& source)
{
  this->header = source.header;
  this->record_index = source.record_index;
  this->record_select = source.record_select;
  this->record_select.set_vector(&(this->record_index));
  this->bwt = source.bwt;
}

//------------------------------------------------------------------------------

size_type
GBWT::runs() const
{
  size_type start = 0, result = 0;
  for(comp_type comp = 0; comp < this->effective(); comp++)
  {
    size_type limit = (comp + 1 < this->effective() ? this->record_select(comp + 2) : this->record_index.size());
    CompressedRecord record(this->bwt, start, limit);
    result += record.runs();
    start = limit;
  }
  return result;
}

//------------------------------------------------------------------------------

edge_type
GBWT::LF(node_type from, size_type i) const
{
  if(from >= this->sigma()) { return invalid_edge(); }
  return this->record(from).LF(i);
}

size_type
GBWT::LF(node_type from, size_type i, node_type to) const
{
  if(to >= this->sigma()) { return invalid_offset(); }
  if(from >= this->sigma()) { return this->count(to); }
  return this->record(from).LF(i, to);
}

range_type
GBWT::LF(node_type from, range_type range, node_type to) const
{
  if(to >= this->sigma()) { return Range::empty_range(); }
  if(from >= this->sigma()) { range.first = range.second = this->count(to); }
  return this->record(from).LF(range, to);
}

//------------------------------------------------------------------------------

CompressedRecord
GBWT::record(node_type node) const
{
  comp_type comp = (node == 0 ? node : node - this->header.offset);
  size_type start = this->record_select(comp + 1);
  size_type limit = (comp + 1 < this->effective() ? this->record_select(comp + 2) : this->record_index.size());
  return CompressedRecord(this->bwt, start, limit);
}

//------------------------------------------------------------------------------

} // namespace gbwt
