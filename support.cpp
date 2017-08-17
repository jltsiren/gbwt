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

#include "internal.h"

namespace gbwt
{

//------------------------------------------------------------------------------

void
DynamicRecord::clear()
{
  DynamicRecord temp;
  this->swap(temp);
}

void
DynamicRecord::swap(DynamicRecord& another)
{
  if(this != &another)
  {
    std::swap(this->body_size, another.body_size);
    this->incoming.swap(another.incoming);
    this->outgoing.swap(another.outgoing);
    this->body.swap(another.body);
  }
}

//------------------------------------------------------------------------------

void
DynamicRecord::recode()
{
  if(this->empty()) { return; }

  bool sorted = true;
  for(rank_type outrank = 1; outrank < this->outdegree(); outrank++)
  {
    if(this->successor(outrank) < this->successor(outrank - 1)) { sorted = false; break; }
  }
  if(sorted) { return; }

  for(run_type& run : this->body) { run.first = this->successor(run.first); }
  sequentialSort(this->outgoing.begin(), this->outgoing.end());
  for(run_type& run : this->body) { run.first = this->edgeTo(run.first); }
}

//------------------------------------------------------------------------------

edge_type
DynamicRecord::LF(size_type i) const
{
  if(i >= this->size()) { return invalid_edge(); }

  std::vector<edge_type> result(this->outgoing);
  rank_type last_edge = 0;
  size_type offset = 0;
  for(run_type run : this->body)
  {
    last_edge = run.first;
    result[run.first].second += run.second;
    offset += run.second;
    if(offset > i) { break; }
  }

  result[last_edge].second -= (offset - i);
  return result[last_edge];
}

size_type
DynamicRecord::LF(size_type i, node_type to) const
{
  size_type outrank = this->edgeTo(to);
  if(outrank >= this->outdegree()) { return invalid_offset(); }

  size_type result = this->offset(outrank), offset = 0;
  for(run_type run : this->body)
  {
    if(run.first == outrank) { result += run.second; }
    offset += run.second;
    if(offset >= i)
    {
      if(run.first == outrank) { result -= offset - i; }
      break;
    }
  }
  return result;
}

range_type
DynamicRecord::LF(range_type range, node_type to) const
{
  if(Range::empty(range)) { return Range::empty_range(); }

  size_type outrank = this->edgeTo(to);
  if(outrank >= this->outdegree()) { return Range::empty_range(); }

  std::vector<run_type>::const_iterator iter = this->body.begin();
  run_type run = *iter;
  size_type result = this->offset(outrank) + (run.first == outrank ? run.second : 0), offset = run.second;

  while(iter != body.end() && offset < range.first)
  {
    ++iter;
    if(iter == body.end()) { break; }
    run = *iter;
    if(run.first == outrank) { result += run.second; }
    offset += run.second;
  }
  range.first = result - (run.first == outrank ? offset - range.first : 0);

  while(iter != body.end() && offset < range.second)
  {
    ++iter;
    if(iter == body.end()) { break; }
    run = *iter;
    if(run.first == outrank) { result += run.second; }
    offset += run.second;
  }
  range.second = result - (run.first == outrank ? offset - range.second : 0);

  return range;
}

node_type
DynamicRecord::operator[](size_type i) const
{
  if(i >= this->size()) { return ENDMARKER; }

  size_type offset = 0;
  for(run_type run : this->body)
  {
    offset += run.second;
    if(offset > i) { return this->successor(run.first); }
  }

  return ENDMARKER;
}

//------------------------------------------------------------------------------

rank_type
DynamicRecord::edgeTo(node_type to) const
{
  for(rank_type outrank = 0; outrank < this->outdegree(); outrank++)
  {
    if(this->successor(outrank) == to) { return outrank; }
  }
  return this->outdegree();
}

//------------------------------------------------------------------------------

rank_type
DynamicRecord::findFirst(node_type from) const
{
  for(size_type i = 0; i < this->indegree(); i++)
  {
    if(this->incoming[i].first >= from) { return i; }
  }
  return this->indegree();
}

void
DynamicRecord::increment(node_type from)
{
  for(rank_type inrank = 0; inrank < this->indegree(); inrank++)
  {
    if(this->predecessor(inrank) == from) { this->count(inrank)++; return; }
  }
  this->addIncoming(edge_type(from, 1));
}

void
DynamicRecord::addIncoming(edge_type inedge)
{
  this->incoming.push_back(inedge);
  sequentialSort(this->incoming.begin(), this->incoming.end());
}

//------------------------------------------------------------------------------

std::ostream&
operator<<(std::ostream& out, const DynamicRecord& record)
{
  out << "(size " << record.size() << ", "
      << record.runs() << " runs, "
      << "indegree " << record.indegree()
      << ", outdegree " << record.outdegree()
      << ", incoming = " << record.incoming
      << ", outgoing = " << record.outgoing
      << ", body = " << record.body << ")";

  return out;
}

//------------------------------------------------------------------------------

CompressedRecord::CompressedRecord(const std::vector<byte_type>& source, size_type start, size_type limit)
{
  this->outgoing.resize(ByteCode::read(source, start));
  node_type prev = 0;
  for(edge_type& outedge : this->outgoing)
  {
    outedge.first = ByteCode::read(source, start) + prev;
    prev = outedge.first;
    outedge.second = ByteCode::read(source, start);
  }

  this->body = source.data() + start;
  this->data_size = limit - start;
}

size_type
CompressedRecord::size() const
{
  size_type result = 0;
  if(this->outdegree() > 0)
  {
    for(CompressedRecordIterator iter(*this); !(iter.end()); ++iter) { result += iter->second; }
  }
  return result;
}

size_type
CompressedRecord::runs() const
{
  size_type result = 0;
  if(this->outdegree() > 0)
  {
    for(CompressedRecordIterator iter(*this); !(iter.end()); ++iter) { result++; }
  }
  return result;
}

edge_type
CompressedRecord::LF(size_type i) const
{
  if(this->outdegree() == 0) { return invalid_edge(); }

  for(CompressedRecordFullIterator iter(*this); !(iter.end()); ++iter)
  {
    if(iter.offset() > i)
    {
      edge_type result = iter.edge(); result.second -= (iter.offset() - i);
      return result;
    }
  }
  return invalid_edge();
}

size_type
CompressedRecord::LF(size_type i, node_type to) const
{
  size_type outrank = this->edgeTo(to);
  if(outrank >= this->outdegree()) { return invalid_offset(); }
  CompressedRecordRankIterator iter(*this, outrank);

  while(!(iter.end()) && iter.offset() < i) { ++iter; }
  return iter.rankAt(i);
}

range_type
CompressedRecord::LF(range_type range, node_type to) const
{
  if(Range::empty(range)) { return Range::empty_range(); }

  size_type outrank = this->edgeTo(to);
  if(outrank >= this->outdegree()) { return Range::empty_range(); }
  CompressedRecordRankIterator iter(*this, outrank);

  while(!(iter.end()) && iter.offset() < range.first) { ++iter; }
  range.first = iter.rankAt(range.first);
  while(!(iter.end()) && iter.offset() < range.second) { ++iter; }
  range.second = iter.rankAt(range.second);

  return range;
}

node_type
CompressedRecord::operator[](size_type i) const
{
  if(this->outdegree() == 0) { return ENDMARKER; }

  for(CompressedRecordIterator iter(*this); !(iter.end()); ++iter)
  {
    if(iter.offset() > i) { return this->successor(iter->first); }
  }
  return ENDMARKER;
}

rank_type
CompressedRecord::edgeTo(node_type to) const
{
  for(rank_type outrank = 0; outrank < this->outdegree(); outrank++)
  {
    if(this->successor(outrank) == to) { return outrank; }
  }
  return this->outdegree();
}

//------------------------------------------------------------------------------

} // namespace gbwt
