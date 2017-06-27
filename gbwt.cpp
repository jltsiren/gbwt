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

namespace gbwt
{

//------------------------------------------------------------------------------

Sequence::Sequence() :
  id(0), curr(ENDMARKER), next(ENDMARKER), offset(0), pos(0)
{
}

Sequence::Sequence(const text_type& text, size_type i, size_type seq_id) :
  id(seq_id), curr(ENDMARKER), next(text[i]), offset(seq_id), pos(i)
{
}

//------------------------------------------------------------------------------

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

/*
  A support structure for run-length encoding outrank sequences.
*/

struct RunMerger
{
  typedef Run::value_type value_type;
  typedef Run::run_type   run_type;

  size_type              total_size;
  run_type               accumulator;
  std::vector<run_type>  runs;
  std::vector<size_type> counts;

  RunMerger(size_type sigma) : total_size(0), accumulator(0, 0), counts(sigma) {}

  inline size_type size() const { return this->total_size; }

  inline void insert(run_type run)
  {
    this->total_size += run.second; counts[run.first] += run.second;
    if(run.first == accumulator.first) { accumulator.second += run.second; }
    else { this->flush(); this->accumulator = run; }
  }

  inline void insert(value_type value) { this->insert(run_type(value, 1)); }

  inline void flush()
  {
    if(this->accumulator.second > 0)
    {
      this->runs.push_back(this->accumulator);
      this->accumulator.second = 0;
    }
  }

  inline void addEdge() { this->counts.push_back(0); }

  // Flush the merger and swap body and size with the record. Counts are invalid after this.
  void swap(DynamicRecord& record);
};

void
RunMerger::swap(DynamicRecord& record)
{
  this->flush();
  this->runs.swap(record.body);
  std::swap(this->total_size, record.body_size);
}

//------------------------------------------------------------------------------

size_type
DynamicGBWT::LF(node_type from, size_type i, node_type to) const
{
  if(to >= this->sigma()) { return this->size(); }
  if(from >= this->sigma()) { return this->bwt[to].size(); }

  const DynamicRecord& from_node = this->bwt[from];
  rank_type outrank = from_node.edgeTo(to);
  if(outrank < from_node.outdegree())
  {
    return from_node.LF(i, outrank);
  }

  /*
    Edge (from, to) has not been observed. We find the first edge from a node >= 'from' to 'to'.
    If 'inrank' is equal to indegree, all incoming edges are from nodes < 'from'.
    Otherwise the result is the stored offset in the node we found.
  */
  const DynamicRecord& to_node = this->bwt[to];
  rank_type inrank = to_node.findFirst(from);
  if(inrank >= to_node.indegree()) { return to_node.size(); }
  const DynamicRecord& next_from = this->bwt[to_node.predecessor(inrank)];
  return next_from.offset(next_from.edgeTo(to));
}

//------------------------------------------------------------------------------

void
DynamicGBWT::insert(const text_type& text)
{
  if(text.empty()) { return; }
  if(text[text.size() - 1] != ENDMARKER)
  {
    std::cerr << "DynamicGBWT::insert(): The text must end an endmarker" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  for(node_type node : text)
  {
    if(node >= this->sigma())
    {
      std::cerr << "DynamicGBWT::insert(): Cannot insert " << node << " with alphabet size " << this->sigma() << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  std::vector<Sequence> sequences;
  bool insert = true;
  for(size_type i = 0; i < text.size(); i++)
  {
    if(insert)
    {
      Sequence temp(text, i, this->count(text[0]) + sequences.size());
      sequences.push_back(temp); insert = false;
    }
    if(text[i] == ENDMARKER) { insert = true; }
  }

  // Invariant: Sequences are sorted by (curr, offset).
  while(true)
  {
    /*
      Process ranges of sequences sharing the same 'curr' node.
      - Add the outgoing edge (curr, next) if necessary.
      - Insert the 'next' node into position 'offset' in the body.
      - Set 'offset' to rank(next) within the record.
      - Update the predecessor count of 'curr' in the incoming edges of 'next'.
    */
    size_type i = 0;
    while(i < sequences.size())
    {
      node_type curr = sequences[i].curr;
      DynamicRecord& current = this->bwt[curr];
      RunMerger new_body(current.outdegree());
      std::vector<run_type>::iterator iter = current.body.begin();
      while(i < sequences.size() && sequences[i].curr == curr)
      {
        rank_type outrank = current.edgeTo(sequences[i].next);
        if(outrank >= current.outdegree())  // Add edge (curr, next) if it does not exist.
        {
          current.outgoing.push_back(edge_type(sequences[i].next, 0));
          new_body.addEdge();
        }
        while(new_body.size() < sequences[i].offset)  // Add old runs until 'offset'.
        {
          if(iter->second <= sequences[i].offset - new_body.size()) { new_body.insert(*iter); ++iter; }
          else
          {
            run_type temp(iter->first, sequences[i].offset - new_body.size());
            new_body.insert(temp);
            iter->second -= temp.second;
          }
        }
        sequences[i].offset = new_body.counts[outrank]; // rank(next) within the record.
        new_body.insert(outrank);
        this->bwt[sequences[i].next].increment(curr);
        i++;
      }
      new_body.swap(current);
    }

    /*
      Sort the sequences for the next iteration and remove the ones that have reached the endmarker.
      Note that sorting by (next, curr, offset) now is equivalent to sorting by (curr, offset) in the
      next interation.
    */
    std::sort(sequences.begin(), sequences.end());
    size_type head = 0;
    while(head < sequences.size() && sequences[head].next == ENDMARKER) { head++; }
    if(head > 0)
    {
      for(size_type j = 0; head + j < sequences.size(); j++) { sequences[j] = sequences[head + j]; }
      sequences.resize(sequences.size() - head);
    }
    if(sequences.empty()) { break; }

    /*
      Rebuild the edge offsets in the outgoing edges to each 'next' node. The offsets will be
      valid after the insertions in the next iteration.
    */
    node_type next = this->sigma();
    for(Sequence& seq : sequences)
    {
      if(seq.next == next) { continue; }
      next = seq.next;
      size_type offset = 0;
      for(edge_type inedge : this->bwt[next].incoming)
      {
        DynamicRecord& predecessor = this->bwt[inedge.first];
        predecessor.offset(predecessor.edgeTo(next)) = offset;
        offset += inedge.second;
      }
    }

    /*
      Until now sequence offsets have been rank(next) within the record. We add edge offsets
      to them to get valid offsets in the next record and then advance the text position.
    */
    for(Sequence& seq : sequences)
    {
      const DynamicRecord& current = this->bwt[seq.curr];
      seq.offset += current.offset(current.edgeTo(seq.next));
      seq.advance(text);
    }
  }
}

//------------------------------------------------------------------------------

} // namespace gbwt
