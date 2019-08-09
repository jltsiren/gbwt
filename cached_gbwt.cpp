/*
  Copyright (c) 2019 Jouni Siren

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
constexpr double CachedGBWT::MAX_LOAD_FACTOR;

//------------------------------------------------------------------------------

CachedGBWT::CachedGBWT(const GBWT& gbwt_index) :
  index(gbwt_index), cache_index(INITIAL_CAPACITY, invalid_edge())
{
  this->cached_records.reserve(INITIAL_CAPACITY);
}

CachedGBWT::~CachedGBWT()
{
}

//------------------------------------------------------------------------------

void CachedGBWT::clearCache() const
{
  for (edge_type& cell : this->cache_index)
  {
    cell = invalid_edge();
  }
  this->cached_records.clear();
}

size_type CachedGBWT::findRecord(node_type node) const
{
  size_type index_offset = this->findOffset(node);
  if(this->cache_index[index_offset].first == node) { return this->cache_index[index_offset].second; }

  // Insert the new record into the cache. Rehash if needed.
  this->cache_index[index_offset] = edge_type(node, this->cacheSize());
  this->cached_records.emplace_back(this->index.record(node));
  if(this->cacheSize() > MAX_LOAD_FACTOR * this->cacheCapacity()) { this->rehash(); }

  return this->cacheSize() - 1;
}

//------------------------------------------------------------------------------

size_type
CachedGBWT::findOffset(node_type node) const
{
  size_type offset = wang_hash_64(node) & (this->cacheCapacity() - 1);
  for(size_type attempt = 0; attempt < this->cacheCapacity(); attempt++)
  {
    if(this->cache_index[offset].first >= this->cacheSize() || this->cache_index[offset].first == node) { return offset; }

    // Quadratic probing with triangular numbers.
    offset = (offset + attempt + 1) & (this->cacheCapacity() - 1);
  }

  // This should not happen.
  std::cerr << "CachedGBWT::findOffset(): Cannot find index offset for node " << node << std::endl;
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
    size_type offset = this->findOffset(old_cache_index[i].first);
    this->cache_index[offset] = old_cache_index[i];
  }
}

const CompressedRecord&
CachedGBWT::record(node_type node) const
{
  size_type cache_offset = this->findRecord(node);
  return this->cached_records[cache_offset];
}

//------------------------------------------------------------------------------

} // namespace gbwt
