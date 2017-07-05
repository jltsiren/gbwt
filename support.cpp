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

#include "support.h"

namespace gbwt
{

//------------------------------------------------------------------------------

bool
DynamicRecord::operator==(const DynamicRecord& another) const
{
  if(this->size() != another.size() || this->runs() != another.runs() ||
     this->indegree() != another.indegree() || this->outdegree() != another.outdegree())
  {
    return false;
  }

  for(rank_type inrank = 0; inrank < this->indegree(); inrank++)
  {
    if(this->incoming[inrank] != another.incoming[inrank]) { return false; }
  }
  for(rank_type outrank = 0; outrank < this->outdegree(); outrank++)
  {
    if(this->outgoing[outrank] != another.outgoing[outrank]) { return false; }
  }
  for(size_type i = 0; i < this->runs(); i++)
  {
    if(this->body[i] != another.body[i]) { return false; }
  }

  return true;
}

//------------------------------------------------------------------------------

void
DynamicRecord::recode()
{
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

size_type
DynamicRecord::LF(size_type i, rank_type outrank) const
{
  size_type res = this->offset(outrank);
  if(i == 0) { return res; }

  size_type j = 0;
  for(run_type run : this->body)
  {
    if(run.first == outrank) { res += run.second; }
    j += run.second;
    if(j + 1 >= i)
    {
      if(run.first == outrank) { res -= j + 1 - i; }
      break;
    }
  }
  return res;
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

} // namespace gbwt
