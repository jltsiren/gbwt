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

#include <gtest/gtest.h>

#include <gbwt/cached_gbwt.h>
#include <gbwt/dynamic_gbwt.h>

using namespace gbwt;

namespace
{

//------------------------------------------------------------------------------

/*
  TODO Also test DynamicGBWT vs GBWT and GBWT vs correct answers.
  TODO Test using multiple indexes.
*/

//------------------------------------------------------------------------------

/*
  The test paths were copied from the VG gbwt_helper.cpp unit tests.
*/

vector_type alt_path
{
  static_cast<vector_type::value_type>(Node::encode(1, false)),
  static_cast<vector_type::value_type>(Node::encode(2, false)),
  static_cast<vector_type::value_type>(Node::encode(4, false)),
  static_cast<vector_type::value_type>(Node::encode(5, false)),
  static_cast<vector_type::value_type>(Node::encode(6, false)),
  static_cast<vector_type::value_type>(Node::encode(8, false)),
  static_cast<vector_type::value_type>(Node::encode(9, false))
};

vector_type short_path
{
  static_cast<vector_type::value_type>(Node::encode(1, false)),
  static_cast<vector_type::value_type>(Node::encode(4, false)),
  static_cast<vector_type::value_type>(Node::encode(5, false)),
  static_cast<vector_type::value_type>(Node::encode(6, false)),
  static_cast<vector_type::value_type>(Node::encode(7, false)),
  static_cast<vector_type::value_type>(Node::encode(9, false))
};

// Build a bidirectional GBWT of the paths.
GBWT
buildGBWT(const std::vector<vector_type>& paths)
{
  size_type node_width = 1, total_length = 0;
  for(auto& path : paths)
  {
    for(auto node : path) { node_width = std::max(node_width, size_type(sdsl::bits::length(Node::encode(node, true)))); }
    total_length += 2 * (path.size() + 1);
  }

  Verbosity::set(Verbosity::SILENT);
  GBWTBuilder builder(node_width, total_length);
  for(auto& path : paths) { builder.insert(path, true); }
  builder.finish();

  return GBWT(builder.index);
}

// Build a bidirectional GBWT with three paths including a duplicate.
GBWT
getGBWT()
{
  std::vector<vector_type> paths
  {
    short_path, alt_path, short_path
  };

  return buildGBWT(paths);
}

std::vector<vector_type>
getPatterns()
{
  std::vector<vector_type> result;
  result.emplace_back(alt_path.begin(), alt_path.begin() + 2);
  result.emplace_back(short_path.begin(), short_path.begin() + 2);
  result.emplace_back(alt_path.begin() + 2, alt_path.begin() + 5); // Shared.
  result.emplace_back(alt_path.begin() + 2, alt_path.begin() + 6);
  result.emplace_back(short_path.begin() + 2, short_path.begin() + 5);
  return result;
}

//------------------------------------------------------------------------------

TEST(CachedGBWTTest, Statistics)
{
  GBWT index = getGBWT();
  CachedGBWT cached(index);

  EXPECT_EQ(cached.size(), index.size()) << "Wrong index size";
  EXPECT_EQ(cached.empty(), index.empty()) << "Wrong empty flag";
  EXPECT_EQ(cached.sequences(), index.sequences()) << "Wrong number of sequences";
  EXPECT_EQ(cached.sigma(), index.sigma()) << "Wrong alphabet size";
  EXPECT_EQ(cached.effective(), index.effective()) << "Wrong effective alphabet size";

  EXPECT_EQ(cached.runs(), index.runs()) << "Wrong number of runs";
  EXPECT_EQ(cached.samples(), index.samples()) << "Wrong number of samples";

  EXPECT_EQ(cached.bidirectional(), index.bidirectional()) << "Wrong bidirectional flag";
}

TEST(CachedGBWTTest, HighLevel)
{
  GBWT index = getGBWT();
  CachedGBWT cached(index);
  std::vector<vector_type> patterns = getPatterns();

  // find(node)
  for(node_type node = index.firstNode(); node < index.sigma(); node++)
  {
    EXPECT_EQ(cached.find(node), index.find(node)) << "Wrong find() result for node " << node;
  }

  // find(pattern)
  std::vector<SearchState> query_states;
  for(size_type i = 0; i < patterns.size(); i++)
  {
    const vector_type& pattern = patterns[i];
    SearchState cached_state = cached.find(pattern.begin(), pattern.end());
    SearchState index_state = index.find(pattern.begin(), pattern.end());
    query_states.push_back(index_state);
    EXPECT_EQ(cached_state, index_state) << "Wrong find() result for pattern " << i;
  }

  // locate(position), locate(state)
  for(SearchState state : query_states)
  {
    for(size_type offset = state.range.first; offset <= state.range.second; offset++)
    {
      EXPECT_EQ(cached.locate(state.node, offset), index.locate(state.node, offset)) << "Wrong locate() result for node " << state.node << ", offset " << offset;
    }
    EXPECT_EQ(cached.locate(state), index.locate(state)) << "Wrong locate() result for state " << state;
  }

  // extract(sequence)
  for(size_type i = 0; i < index.sequences(); i++)
  {
    EXPECT_EQ(cached.extract(i), index.extract(i)) << "Wrong extract() result for sequence " << i;
  }
}

TEST(CachedGBWTTest, Bidirectional)
{
  GBWT index = getGBWT();
  CachedGBWT cached(index);
  std::vector<vector_type> patterns = getPatterns();

  for(size_type i = 0; i < patterns.size(); i++)
  {
    const vector_type& pattern = patterns[i];
    size_type midpoint = pattern.size() / 2;
    BidirectionalState cached_state = cached.bdFind(pattern[midpoint]);
    BidirectionalState index_state = index.bdFind(pattern[midpoint]);
    for(size_type j = midpoint + 1; j < pattern.size(); j++)
    {
      cached_state = cached.bdExtendForward(cached_state, pattern[j]);
      index_state = index.bdExtendForward(index_state, pattern[j]);
    }
    for(size_type j = midpoint; j > 0; j--)
    {
      cached_state = cached.bdExtendBackward(cached_state, pattern[j - 1]);
      index_state = index.bdExtendBackward(index_state, pattern[j - 1]);
    }
    EXPECT_EQ(cached_state, index_state) << "Wrong bidirectional search result for pattern " << i;
  }
}

TEST(CachedGBWTTest, Nodes)
{
  GBWT index = getGBWT();
  CachedGBWT cached(index);

  // contains(node)
  for(node_type node = 0; node < index.sigma() + 2; node++)
  {
    ASSERT_EQ(cached.contains(node), index.contains(node)) << "Wrong existence flag for node " << node;
  }

  // hasEdge(from, to)
  for(node_type from = 0; from < index.sigma(); from++)
  {
    for(node_type to = 0; to < index.sigma(); to++)
    {
      EXPECT_EQ(cached.hasEdge(from, to), index.hasEdge(from, to)) << "From existence flag for the edge from " << from << " to " << to;
    }
  }

  // edges(from), nodeSize(node), empty(node)
  for(node_type node = 0; node < index.sigma(); node++)
  {
    if(!(index.contains(node))) { continue; }
    EXPECT_EQ(cached.edges(node), index.edges(node)) << "Wrong edges from node " << node;
    EXPECT_EQ(cached.nodeSize(node), index.nodeSize(node)) << "Wrong size for node " << node;
    EXPECT_EQ(cached.empty(node), index.empty(node)) << "Wrong empty flag for node " << node;
  }

  // firstNode(), toComp(node), toNode(comp)
  ASSERT_EQ(cached.firstNode(), index.firstNode()) << "Wrong identifier for the first real node";
  EXPECT_EQ(cached.toComp(ENDMARKER), index.toComp(ENDMARKER)) << "Wrong comp value for the endmarker";
  for(node_type node = index.firstNode(); node < index.sigma(); node++)
  {
    EXPECT_EQ(cached.toComp(node), index.toComp(node)) << "Wrong comp value for node " << node;
  }
  for(comp_type comp = 0; comp < index.effective(); comp++)
  {
    EXPECT_EQ(cached.toNode(comp), index.toNode(comp)) << "Wrong node for comp value " << comp;
  }
}

TEST(CachedGBWTTest, Navigation)
{
  GBWT index = getGBWT();
  CachedGBWT cached(index);

  // LF(node, i), LF(node, i, to)
  for(node_type node = 0; node < index.sigma(); node++)
  {
    if(!(index.contains(node))) { continue; }
    size_type node_size = index.nodeSize(node);
    for(size_type i = 0; i < node_size; i++)
    {
      EXPECT_EQ(cached.LF(node, i), index.LF(node, i)) << "Wrong LF() result from node " << node << ", offset " << i;
      for(node_type to = 0; to < index.sigma(); to++)
      {
        if(!(index.contains(to))) { continue; }
        EXPECT_EQ(cached.LF(node, i, to), index.LF(node, i, to)) << "Wrong LF() result from node " << ", offset " << i << " to node " << to;
      }
    }
  }

  // LF(state, to), bdLF(state, to, reverse_offset)
  for(node_type node = 0; node < index.sigma(); node++)
  {
    if(!(index.contains(node))) { continue; }
    size_type node_size = index.nodeSize(node);
    for(size_type i = 0; i < node_size; i++)
    {
      range_type range(i, std::min(i + 1, node_size - 1));
      SearchState state(node, range);
      for(node_type to = 0; to < index.sigma(); to++)
      {
        if(!(index.contains(to))) { continue; }
        EXPECT_EQ(cached.LF(state, to), index.LF(state, to)) << "Wrong LF() result from state " << state << " to node " << to;
        size_type index_offset = 0, cached_offset = 0;
        EXPECT_EQ(cached.bdLF(state, to, cached_offset), index.bdLF(state, to, index_offset)) << "Wrong bdLF() result from state " << state << " to node " << to;
        EXPECT_EQ(cached_offset, index_offset) << "Wrong reverse offset after bdLF() from state " << state << " to node " << node;
      }
    }
  }
}

TEST(CachedGBWTTest, Sequences)
{
  GBWT index = getGBWT();
  CachedGBWT cached(index);

  // start(sequence)
  for(size_type i = 0; i < index.sequences(); i++)
  {
    EXPECT_EQ(cached.start(i), index.start(i)) << "Wrong starting position for sequence " << i;
  }

  // tryLocate(node, i)
  for(node_type node = 0; node < index.sigma(); node++)
  {
    if(!(index.contains(node))) { continue; }
    size_type node_size = index.nodeSize(node);
    for(size_type i = 0; i < node_size; i++)
    {
      EXPECT_EQ(cached.tryLocate(node, i), index.tryLocate(node, i)) << "Wrong tryLocate() result for node " << node << ", offset " << i;
    }
  }
}

TEST(CachedGBWTTest, Cache)
{
  GBWT index = getGBWT();
  CachedGBWT cached(index);

  // Empty cache.
  ASSERT_EQ(cached.cacheSize(), static_cast<size_type>(0)) << "The cache does not start empty";
  ASSERT_GT(cached.cacheCapacity(), static_cast<size_type>(0)) << "The cache has no capacity";

  // Small cache.
  {
    CachedGBWT small(index, true);
    ASSERT_EQ(small.cacheSize(), static_cast<size_type>(0)) << "Small cache does not start empty";
    ASSERT_GT(small.cacheCapacity(), static_cast<size_type>(0)) << "Small cache has no capacity";
  }

  // Insert all records twice and then clear the cache.
  for(node_type node = cached.firstNode(); node < cached.sigma(); node++)
  {
    cached.findRecord(node);
    cached.findRecord(node);
  }
  ASSERT_EQ(cached.cacheSize(), cached.effective() - 1) << "Wrong number of records in the cache";
  cached.clearCache();
  ASSERT_EQ(cached.cacheSize(), static_cast<size_type>(0)) << "Cleared cache is not empty";
  ASSERT_GT(cached.cacheCapacity(), static_cast<size_type>(0)) << "Cleared cache has no capacity";

  // outdegree(), successor()
  for(node_type node = cached.firstNode(); node < cached.sigma(); node++)
  {
    size_type offset = cached.findRecord(node);
    ASSERT_EQ(cached.outdegree(offset), index.record(node).outdegree()) << "Wrong outdegree for node " << node;
    std::vector<edge_type> edges = cached.edges(node);
    for(size_type i = 0; i < cached.outdegree(offset); i++)
    {
      EXPECT_EQ(cached.successor(offset, i), edges[i].first) << "Wrong successor " << i << " for node " << node;
    }
  }

  // cachedExtend(), cachedExtendForward(), cachedExtendBackward()
  for(node_type node = cached.firstNode(); node < cached.sigma(); node++)
  {
    BidirectionalState state = cached.bdFind(node);
    {
      size_type offset = cached.findRecord(state.forward.node);
      std::vector<edge_type> edges = cached.edges(state.forward.node);
      for(size_type i = 0; i < cached.outdegree(offset); i++)
      {
        if(cached.successor(offset, i) == ENDMARKER) { continue; }
        BidirectionalState correct = cached.bdExtendForward(state, edges[i].first);
        EXPECT_EQ(cached.cachedExtend(state.forward, offset, i), correct.forward) << "Wrong cachedExtend() result from state " << state.forward << " to successor " << i << " that is " << cached.successor(offset, i);
        EXPECT_EQ(cached.cachedExtendForward(state, offset, i), correct) << "Wrong cachedExtendForward() result from state " << state << " to successor " << i << " that is " << cached.successor(offset, i);
      }
    }
    {
      size_type offset = cached.findRecord(state.backward.node);
      std::vector<edge_type> edges = cached.edges(state.backward.node);
      for(size_type i = 0; i < cached.outdegree(offset); i++)
      {
        if(cached.successor(offset, i) == ENDMARKER) { continue; }
        // bdExtendBackward reverses the node, so we need to reverse it here as well, because the orientation
        // is already correct.
        BidirectionalState correct = cached.bdExtendBackward(state, Node::reverse(edges[i].first));
        EXPECT_EQ(cached.cachedExtendBackward(state, offset, i), correct) << "Wrong cachedExtendBackward() result from state " << state << " to successor " << i << " that is " << cached.successor(offset, i);
      }
    }
  }
}

//------------------------------------------------------------------------------

} // namespace
