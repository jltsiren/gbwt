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

DynamicGBWT::size_type
DynamicGBWT::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->header.serialize(out);

  std::vector<byte_type> compressed_bwt;
  std::vector<size_type> bwt_offsets(this->sigma());
  for(node_type node = 0; node < this->sigma(); node++)
  {
    bwt_offsets[node] = compressed_bwt.size();
    const DynamicRecord& current = this->bwt[node];

    // Write the outgoing edges.
    ByteCode::write(compressed_bwt, current.outdegree());
    for(edge_type outedge : current.outgoing)
    {
      ByteCode::write(compressed_bwt, outedge.first);
      ByteCode::write(compressed_bwt, outedge.second);
    }

    // Write the body.
    if(current.outdegree() > 0)
    {
      Run encoder(current.outdegree());
      for(run_type run : current.body) { encoder.write(compressed_bwt, run); }
    }
  }

  // Build and serialize index.
  sdsl::sd_vector_builder builder(compressed_bwt.size(), bwt_offsets.size());
  for(size_type offset : bwt_offsets) { builder.set(offset); }
  sdsl::sd_vector<> node_index(builder);
  sdsl::sd_vector<>::select_1_type node_select(&node_index);
  written_bytes += node_index.serialize(out);
  written_bytes += node_select.serialize(out);

  // Serialize BWT.
  out.write((char*)(compressed_bwt.data()), compressed_bwt.size());
  written_bytes += compressed_bwt.size();

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
  this->bwt.resize(this->sigma());

  // Read the node index.
  sdsl::sd_vector<> node_index; node_index.load(in);
  sdsl::sd_vector<>::select_1_type node_select; node_select.load(in, &node_index);

  for(node_type node = 0; node < this->sigma(); node++)
  {
    DynamicRecord& current = this->bwt[node];
    current.incoming.clear(); // We rebuild the incoming edges later.

    // Read the current node.
    size_type start = node_select(node + 1);
    size_type stop = (node + 1 < this->sigma() ? node_select(node + 2) : node_index.size());
    std::vector<byte_type> node_encoding(stop - start);
    in.read((char*)(node_encoding.data()), node_encoding.size());
    size_type offset = 0;

    // Decompress the outgoing edges.
    current.outgoing.resize(ByteCode::read(node_encoding, offset));
    for(edge_type& outedge : current.outgoing)
    {
      outedge.first = ByteCode::read(node_encoding, offset);
      outedge.second = ByteCode::read(node_encoding, offset);
    }

    // Decompress the body.
    current.body.clear();
    Run decoder(current.outdegree());
    while(offset < node_encoding.size()) { current.body.push_back(decoder.read(node_encoding, offset)); }
  }

  // Rebuild the incoming edges.
  for(node_type node = 0; node < this->sigma(); node++)
  {
    DynamicRecord& current = this->bwt[node];
    std::vector<size_type> counts(current.outdegree());
    for(run_type run : current.body) { counts[run.first] += run.second; }
    for(rank_type outrank = 0; outrank < current.outdegree(); outrank++)
    {
      if(current.successor(outrank) != ENDMARKER)
      {
        this->bwt[current.successor(outrank)].addIncoming(edge_type(node, counts[outrank]));
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
  for(const DynamicRecord& record : this->bwt) { total += record.runs(); }
  return total;
}

bool
DynamicGBWT::compare(const DynamicGBWT& another, std::ostream& out) const
{
  out << "Comparing dynamic GBWTs" << std::endl;
  out << std::endl;

  if(this->header != another.header)
  {
    out << "This:    " << this->header << std::endl;
    out << "Another: " << another.header << std::endl;
    out << std::endl;
    return false;
  }

  for(node_type node = 0; node < this->sigma(); node++)
  {
    if(this->bwt[node] != another.bwt[node])
    {
      out << "This[" << node << "]:    " << this->bwt[node] << std::endl;
      out << "Another[" << node << "]: " << another.bwt[node] << std::endl;
      out << std::endl;
      return false;
    }
  }

  out << "The GBWTs are identical" << std::endl;
  out << std::endl;
  return true;
}

//------------------------------------------------------------------------------

void
swapBody(DynamicRecord& record, RunMerger& merger)
{
  merger.flush();
  merger.runs.swap(record.body);
  std::swap(merger.total_size, record.body_size);
}

void
DynamicGBWT::insert(const text_type& text)
{
  double start = readTimer();

  if(text.empty()) { return; }
  if(text[text.size() - 1] != ENDMARKER)
  {
    std::cerr << "DynamicGBWT::insert(): The text must end an endmarker" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  /*
    Find the start of each sequence and initialize the sequence objects at the endmarker node.
    Increase alphabet size if necessary.
  */
  bool seq_start = true;
  size_type max_node = 0;
  std::vector<Sequence> seqs;
  for(size_type i = 0; i < text.size(); i++)
  {
    max_node = std::max(text[i], max_node);
    if(seq_start)
    {
      Sequence temp(text, i, this->sequences());
      seqs.push_back(temp); seq_start = false;
      this->header.sequences++;
    }
    if(text[i] == ENDMARKER) { seq_start = true; }
  }
  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    std::cerr << "DynamicGBWT::insert(): Inserting " << seqs.size() << " sequences of total length " << text.size() << std::endl;
  }
  if(max_node >= this->sigma())
  {
    if(Verbosity::level >= Verbosity::FULL)
    {
      std::cerr << "DynamicGBWT::insert(): Increasing alphabet size to " << (max_node + 1) << std::endl;
    }
    this->header.alphabet_size = max_node + 1;
    this->bwt.resize(this->sigma());
  }

  // Invariant: Sequences are sorted by (curr, offset).
  size_type iteration = 0;
  while(true)
  {
    iteration++;  // We use 1-based iterations.

    /*
      Process ranges of sequences sharing the same 'curr' node.
      - Add the outgoing edge (curr, next) if necessary.
      - Insert the 'next' node into position 'offset' in the body.
      - Set 'offset' to rank(next) within the record.
      - Update the predecessor count of 'curr' in the incoming edges of 'next'.

      We do not maintain incoming edges to the endmarker, because it can be expensive
      and because searching with the endmarker does not work in a multi-string BWT.
    */
    size_type i = 0;
    while(i < seqs.size())
    {
      node_type curr = seqs[i].curr;
      DynamicRecord& current = this->bwt[curr];
      RunMerger new_body(current.outdegree());
      std::vector<run_type>::iterator iter = current.body.begin();
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
        seqs[i].offset = new_body.counts[outrank]; // rank(next) within the record.
        new_body.insert(outrank);
        if(seqs[i].next != ENDMARKER)  // The endmarker does not have incoming edges.
        {
          this->bwt[seqs[i].next].increment(curr);
        }
        i++;
      }
      while(iter != current.body.end()) // Add the rest of the old body.
      {
        new_body.insert(*iter); ++iter;
      }
      swapBody(current, new_body);
    }
    this->header.size += seqs.size();

    /*
      Sort the sequences for the next iteration and remove the ones that have reached the endmarker.
      Note that sorting by (next, curr, offset) now is equivalent to sorting by (curr, offset) in the
      next interation.
    */
    parallelMergeSort(seqs.begin(), seqs.end());
    size_type head = 0;
    while(head < seqs.size() && seqs[head].next == ENDMARKER) { head++; }
    if(head > 0)
    {
      for(size_type j = 0; head + j < seqs.size(); j++) { seqs[j] = seqs[head + j]; }
      seqs.resize(seqs.size() - head);
    }
    if(seqs.empty()) { break; }

    /*
      Rebuild the edge offsets in the outgoing edges to each 'next' node. The offsets will be
      valid after the insertions in the next iteration.
    */
    node_type next = this->sigma();
    for(Sequence& seq : seqs)
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
    for(Sequence& seq : seqs)
    {
      const DynamicRecord& current = this->bwt[seq.curr];
      seq.offset += current.offset(current.edgeTo(seq.next));
      seq.advance(text);
    }
  }

  // Update the effective alphabet size and sort the outgoing edges.
  std::atomic<size_type> effective_alphabet(0);
  #pragma omp parallel for schedule(static)
  for(node_type node = 0; node < this->sigma(); node++)
  {
    if(this->count(node) > 0)
    {
      effective_alphabet++;
      this->bwt[node].recode();
    }
  }
  this->header.nodes = effective_alphabet;
  if(Verbosity::level >= Verbosity::FULL)
  {
    std::cerr << "DynamicGBWT::insert(): Effective alphabet size " << this->effective() << std::endl;
  }

  if(Verbosity::level >= Verbosity::EXTENDED)
  {
    double seconds = readTimer() - start;
    std::cerr << "DynamicGBWT::insert(): " << iteration << " iterations in " << seconds << " seconds" << std::endl;
  }
}

//------------------------------------------------------------------------------

size_type
DynamicGBWT::LF(node_type from, size_type i, node_type to) const
{
  if(to >= this->sigma()) { return invalid_offset(); }
  if(from >= this->sigma()) { return this->count(to); }

  size_type result = this->bwt[from].LF(i, to);
  if(result != invalid_offset()) { return result; }

  /*
    Edge (from, to) has not been observed. We find the first edge from a node >= 'from' to 'to'.
    If 'inrank' is equal to indegree, all incoming edges are from nodes < 'from'.
    Otherwise the result is the stored offset in the node we found.
  */
  const DynamicRecord& to_node = this->bwt[to];
  rank_type inrank = to_node.findFirst(from);
  if(inrank >= to_node.indegree()) { return this->count(to); }
  const DynamicRecord& next_from = this->bwt[to_node.predecessor(inrank)];
  return next_from.offset(next_from.edgeTo(to));
}

edge_type
DynamicGBWT::LF(node_type from, size_type i) const
{
  if(from >= this->sigma()) { return invalid_edge(); }
  return this->bwt[from].LF(i);
}

//------------------------------------------------------------------------------

} // namespace gbwt
