/*
  Copyright (c) 2019, 2021 Jouni Siren

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

#include <gbwt/cached_gbwt.h>

namespace gbwt
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_type CachedGBWT::INITIAL_CAPACITY;
constexpr size_type CachedGBWT::SINGLE_CAPACITY;
constexpr double CachedGBWT::MAX_LOAD_FACTOR;

//------------------------------------------------------------------------------

CachedGBWT::CachedGBWT()
{
}

CachedGBWT::CachedGBWT(const CachedGBWT& source)
{
  this->copy(source);
}

CachedGBWT::CachedGBWT(CachedGBWT&& source)
{
  *this = std::move(source);
}

void
CachedGBWT::swap(CachedGBWT& another)
{
  if(this != &another)
  {
    std::swap(this->index, another.index);
    this->cache_index.swap(another.cache_index);
    this->cached_records.swap(another.cached_records);
  }
}

CachedGBWT&
CachedGBWT::operator=(const CachedGBWT& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

CachedGBWT&
CachedGBWT::operator=(CachedGBWT&& source)
{
  if(this != &source)
  {
    this->index = std::move(source.index);
    this->cache_index = std::move(source.cache_index);
    this->cached_records = std::move(source.cached_records);
  }
  return *this;
}

void
CachedGBWT::copy(const CachedGBWT& source)
{
  this->index = source.index;
  this->cache_index = source.cache_index;
  this->cached_records = source.cached_records;
}

CachedGBWT::~CachedGBWT()
{
}

//------------------------------------------------------------------------------

CachedGBWT::CachedGBWT(const GBWT& gbwt_index, bool single_record) :
  index(&gbwt_index), cache_index((single_record ? SINGLE_CAPACITY : INITIAL_CAPACITY), invalid_edge())
{
  this->cached_records.reserve((single_record ? SINGLE_CAPACITY : INITIAL_CAPACITY));
}

//------------------------------------------------------------------------------

void
CachedGBWT::clearCache()
{
  for (edge_type& cell : this->cache_index)
  {
    cell = invalid_edge();
  }
  this->cached_records.clear();
}

size_type
CachedGBWT::findRecord(node_type node) const
{
  size_type index_offset = this->indexOffset(node);
  if(this->cache_index[index_offset].first == node) { return this->cache_index[index_offset].second; }

  // Insert the new record into the cache. Rehash if needed.
  this->cache_index[index_offset] = edge_type(node, this->cacheSize());
  this->cached_records.emplace_back(this->index->record(node));
  if(this->cacheSize() > MAX_LOAD_FACTOR * this->cacheCapacity()) { this->rehash(); }

  return this->cacheSize() - 1;
}

SearchState
CachedGBWT::cachedExtend(SearchState state, size_type cache_offset, size_type i) const
{
  if(state.empty()) { return SearchState(); }
  node_type node = this->successor(cache_offset, i);
  state.range = this->cached_records[cache_offset].LF(state.range, node);
  state.node = node;
  return state;
}

BidirectionalState
CachedGBWT::cachedExtendForward(BidirectionalState state, size_type cache_offset, size_type i) const
{
  if(state.empty()) { return BidirectionalState(); }
  size_type reverse_offset = 0;
  node_type node = this->successor(cache_offset, i);
  state.forward.range = this->cached_records[cache_offset].bdLF(state.forward.range, node, reverse_offset);
  state.forward.node = node;
  state.backward.range.first += reverse_offset;
  state.backward.range.second = state.backward.range.first + state.forward.size() - 1;
  return state;
}

BidirectionalState
CachedGBWT::cachedExtendBackward(BidirectionalState state, size_type cache_offset, size_type i) const
{
  state.flip();
  state = this->cachedExtendForward(state, cache_offset, i);
  state.flip();
  return state;
}

//------------------------------------------------------------------------------

// TODO This should have a common implementation with GBWT::locate().
std::vector<size_type>
CachedGBWT::locate(SearchState state) const
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
    size_type cache_offset = invalid_offset();
    sample_type sample;
    edge_type LF_result;
    range_type LF_range;

    for(size_type i = 0; i < positions.size(); i++)
    {
      if(positions[i].first != curr)              // Node changed.
      {
        curr = positions[i].first; cache_offset = this->findRecord(curr);
        sample = this->index->da_samples.nextSample(this->toComp(curr), positions[i].second);
        LF_range.first = positions[i].second;
        LF_result = this->cached_records[cache_offset].runLF(positions[i].second, LF_range.second);
      }
      if(sample.first < positions[i].second)      // Went past the sample.
      {
        sample = this->index->da_samples.nextSample(this->toComp(curr), positions[i].second);
      }
      if(sample.first > positions[i].second)      // Not sampled, also valid for invalid_sample().
      {
        if(positions[i].second > LF_range.second) // Went past the existing LF() result.
        {
          LF_range.first = positions[i].second;
          LF_result = this->cached_records[cache_offset].runLF(positions[i].second, LF_range.second);
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

// TODO This should have a common implementation with GBWT::inverseLF().
edge_type
CachedGBWT::inverseLF(node_type from, size_type i) const
{
  if(!(this->bidirectional()) || from == ENDMARKER) { return invalid_edge(); }

  // Find the predecessor node id.
  const CompressedRecord& reverse_record = this->record(Node::reverse(from));
  node_type predecessor = reverse_record.predecessorAt(i);
  if(predecessor == invalid_node()) { return invalid_edge(); }

  // Determine the offset.
  const CompressedRecord& pred_record = this->record(predecessor);
  size_type offset = pred_record.offsetTo(from, i);
  if(offset == invalid_offset()) { return invalid_edge(); }

  return edge_type(predecessor, offset);
}

//------------------------------------------------------------------------------

size_type
CachedGBWT::indexOffset(node_type node) const
{
  size_type offset = wang_hash_64(node) & (this->cacheCapacity() - 1);
  for(size_type attempt = 0; attempt < this->cacheCapacity(); attempt++)
  {
    if(this->cache_index[offset].first == node || this->cache_index[offset].second >= this->cacheSize()) { return offset; }

    // Quadratic probing with triangular numbers.
    offset = (offset + attempt + 1) & (this->cacheCapacity() - 1);
  }

  // This should not happen.
  std::cerr << "CachedGBWT::indexOffset(): Cannot find index offset for node " << node << std::endl;
  return invalid_offset();
}

void
CachedGBWT::rehash() const
{
  std::vector<edge_type> old_cache_index(2 * this->cacheCapacity(), invalid_edge());
  this->cache_index.swap(old_cache_index);
  for(size_type i = 0; i < old_cache_index.size(); i++)
  {
    if(old_cache_index[i].second >= this->cacheSize()) { continue; }
    size_type offset = this->indexOffset(old_cache_index[i].first);
    this->cache_index[offset] = old_cache_index[i];
  }
}

//------------------------------------------------------------------------------

} // namespace gbwt
