#include <gtest/gtest.h>

#include <gbwt/dynamic_gbwt.h>
#include <gbwt/fast_locate.h>
#include <gbwt/sequence_locate.h>

#include <map>
#include <random>
#include <sstream>

using namespace gbwt;

namespace
{

//------------------------------------------------------------------------------

// Builds a GBWT of the paths.
GBWT
buildGBWT(const std::vector<vector_type>& paths, bool bidirectional)
{
  size_type node_width = 1, total_length = 0;
  for(auto& path : paths)
  {
    for(auto node : path) { node_width = std::max(node_width, size_type(sdsl::bits::length(Node::encode(node, true)))); }
    total_length += (bidirectional ? 2 : 1) * (path.size() + 1);
  }

  Verbosity::set(Verbosity::SILENT);
  GBWTBuilder builder(node_width, total_length);
  for(auto& path : paths) { builder.insert(path, bidirectional); }
  builder.finish();

  return GBWT(builder.index);
}

// Correct packed positions for each GBWT position, computed by walking the sequences.
std::map<edge_type, size_type>
truePositions(const GBWT& index, const SequenceLocate& r_index)
{
  std::map<edge_type, size_type> result;
  for(size_type i = 0; i < index.sequences(); i++)
  {
    std::vector<edge_type> visits;
    size_type total = 0;
    for(edge_type pos = index.start(i); pos.first != ENDMARKER; pos = index.LF(pos))
    {
      visits.push_back(pos); total += r_index.length(pos.first);
    }
    result[edge_type(ENDMARKER, i)] = r_index.pack(i, total);
    size_type offset = total;
    for(edge_type pos : visits)
    {
      offset -= r_index.length(pos.first);
      result[pos] = r_index.pack(i, offset);
    }
  }
  return result;
}

// Forward positions of every base as (node, offset in node).
std::vector<edge_type>
trueBases(const GBWT& index, const SequenceLocate& r_index, size_type seq_id)
{
  std::vector<edge_type> result;
  for(edge_type pos = index.start(seq_id); pos.first != ENDMARKER; pos = index.LF(pos))
  {
    for(size_type j = 0; j < r_index.length(pos.first); j++) { result.emplace_back(pos.first, j); }
  }
  return result;
}

// Copies an SDSL vector into a std::vector.
template<class T>
std::vector<size_type>
toVector(const T& source)
{
  return std::vector<size_type>(source.begin(), source.end());
}

// Returns the positions of the set bits in `last`.
std::vector<size_type>
lastBits(const SequenceLocate& r_index)
{
  std::vector<size_type> result;
  for(auto iter = r_index.last.one_begin(); iter != r_index.last.one_end(); ++iter) { result.push_back(iter->second); }
  return result;
}

//------------------------------------------------------------------------------

/*
  Nodes a..f = 1..6 with lengths 3, 1, 2, 4, 1, 5.
  P0 = a b d e f (14 bp), P1 = a c e f (11 bp).
*/

class ExampleTest : public ::testing::Test
{
public:
  GBWT index;
  SequenceLocate r_index;

  constexpr static node_type A = 1, B = 2, C = 3, D = 4, E = 5, F = 6;

  void SetUp() override
  {
    std::vector<vector_type> paths { { A, B, D, E, F }, { A, C, E, F } };
    this->index = buildGBWT(paths, false);
    std::vector<size_type> lengths { 0, 3, 1, 2, 4, 1, 5 };
    this->r_index = SequenceLocate(this->index, [lengths](node_type node) -> size_type { return lengths[node]; });
  }
};

TEST_F(ExampleTest, Arrays)
{
  ASSERT_EQ(this->index.effective(), size_type(7)) << "Unexpected alphabet";
  EXPECT_EQ(this->r_index.header.max_length, size_type(15)) << "Invalid max_length";
  EXPECT_EQ(toVector(this->r_index.comp_to_run), std::vector<size_type>({ 0, 1, 3, 4, 5, 6, 7 })) << "Invalid comp_to_run";
  EXPECT_EQ(toVector(this->r_index.samples), std::vector<size_type>({ 14, 11, 23, 10, 21, 6, 20, 15, 0 })) << "Invalid samples";
  EXPECT_EQ(lastBits(this->r_index), std::vector<size_type>({ 0, 5, 6, 10, 11, 15, 21, 23, 26 })) << "Invalid last";
  EXPECT_EQ(toVector(this->r_index.last_to_run), std::vector<size_type>({ 8, 6, 5, 3, 1, 7, 4, 2, 0 })) << "Invalid last_to_run";
  EXPECT_EQ(toVector(this->r_index.sequence_length), std::vector<size_type>({ 14, 11 })) << "Invalid sequence lengths";
}

TEST_F(ExampleTest, NodeOffsetsUnchanged)
{
  // Sanity check: the node-based r-index uses the same runs.
  FastLocate node_index(this->index);
  EXPECT_EQ(node_index.header.max_length, size_type(6)) << "Invalid node max_length";
  EXPECT_EQ(toVector(node_index.samples), std::vector<size_type>({ 5, 4, 9, 3, 8, 2, 7, 6, 0 })) << "Invalid node samples";
  EXPECT_EQ(toVector(node_index.last_to_run), toVector(this->r_index.last_to_run)) << "Different last_to_run";
}

TEST_F(ExampleTest, LocateNext)
{
  // P0 at e is the last entry of e, so the next one is P1 at f.
  EXPECT_EQ(this->r_index.locateNext(this->r_index.pack(0, 5)), this->r_index.pack(1, 0)) << "Invalid locateNext()";
}

TEST_F(ExampleTest, LocatePositions)
{
  SearchState state = this->index.find(E);
  std::vector<std::pair<size_type, size_type>> raw { { 1, 5 }, { 0, 5 } };
  std::vector<std::pair<size_type, size_type>> forward { { 1, 5 }, { 0, 8 } };
  EXPECT_EQ(this->r_index.locatePositions(state), raw) << "Invalid locatePositions()";
  EXPECT_EQ(this->r_index.locateForward(state), forward) << "Invalid locateForward()";
  EXPECT_EQ(this->r_index.locate(state), std::vector<size_type>({ 0, 1 })) << "Invalid locate()";
}

TEST_F(ExampleTest, FindExtend)
{
  vector_type pattern { A, B, D };
  size_type first = SequenceLocate::NO_POSITION;
  SearchState state = this->r_index.find(pattern.begin(), pattern.end(), first);
  EXPECT_EQ(state, SearchState(D, 0, 0)) << "Invalid state for a b d";
  EXPECT_EQ(first, this->r_index.pack(0, 6)) << "Invalid first for a b d";

  pattern = { E, F };
  state = this->r_index.find(pattern.begin(), pattern.end(), first);
  EXPECT_EQ(state, SearchState(F, 0, 1)) << "Invalid state for e f";
  EXPECT_EQ(first, this->r_index.pack(1, 0)) << "Invalid first for e f";
}

TEST_F(ExampleTest, LocateSequence)
{
  EXPECT_EQ(this->r_index.locateSequence(0, 6), edge_type(D, 2)) << "Invalid position for P0:6";
  edge_type gbwt_pos;
  this->r_index.locateSequence(0, 6, gbwt_pos);
  EXPECT_EQ(gbwt_pos, edge_type(D, 0)) << "Invalid GBWT position for P0:6";
  EXPECT_EQ(this->index.LF(gbwt_pos), edge_type(E, 1)) << "Invalid LF from P0:6";
  this->r_index.locateSequence(0, 14, gbwt_pos);
  EXPECT_EQ(gbwt_pos, invalid_edge()) << "GBWT position past the end";
  for(size_type seq_id = 0; seq_id < this->index.sequences(); seq_id++)
  {
    std::vector<edge_type> bases = trueBases(this->index, this->r_index, seq_id);
    ASSERT_EQ(bases.size(), this->r_index.sequenceLength(seq_id)) << "Invalid length for sequence " << seq_id;
    for(size_type i = 0; i < bases.size(); i++)
    {
      EXPECT_EQ(this->r_index.locateSequence(seq_id, i), bases[i]) << "Invalid position for " << seq_id << ":" << i;
    }
    EXPECT_EQ(this->r_index.locateSequence(seq_id, bases.size()), invalid_edge()) << "Position past the end for " << seq_id;
  }
}

TEST_F(ExampleTest, Serialization)
{
  std::stringstream buffer;
  this->r_index.serialize(buffer);
  SequenceLocate loaded;
  loaded.load(buffer);
  loaded.setGBWT(this->index, this->r_index.length);

  EXPECT_EQ(loaded.header.max_length, this->r_index.header.max_length) << "Invalid max_length";
  EXPECT_EQ(toVector(loaded.samples), toVector(this->r_index.samples)) << "Invalid samples";
  EXPECT_EQ(lastBits(loaded), lastBits(this->r_index)) << "Invalid last";
  EXPECT_EQ(toVector(loaded.last_to_run), toVector(this->r_index.last_to_run)) << "Invalid last_to_run";
  EXPECT_EQ(toVector(loaded.comp_to_run), toVector(this->r_index.comp_to_run)) << "Invalid comp_to_run";
  EXPECT_EQ(toVector(loaded.sequence_length), toVector(this->r_index.sequence_length)) << "Invalid sequence lengths";
  EXPECT_EQ(loaded.locateSequence(0, 6), edge_type(D, 2)) << "Invalid query after loading";
}

//------------------------------------------------------------------------------

/*
  Random bidirectional GBWTs with shared prefixes, duplicates, repeated nodes, empty
  paths, and paths ending at the same node. Every query is checked against the
  positions obtained by walking the sequences.
*/

// Generates a random set of paths.
std::vector<vector_type>
randomPaths(std::mt19937_64& rng)
{
  std::vector<vector_type> paths;
  size_type n = 2 + rng() % 8;
  for(size_type i = 0; i < n; i++)
  {
    size_type choice = rng() % 6;
    if(choice == 0 && !paths.empty()) { paths.push_back(paths[rng() % paths.size()]); continue; }
    vector_type path;
    if(choice == 1 && !paths.empty())
    {
      const vector_type& source = paths[rng() % paths.size()];
      path.assign(source.begin(), source.begin() + rng() % (source.size() + 1));
    }
    if(choice == 2) { paths.push_back(path); continue; }
    size_type curr = (path.empty() ? 1 + rng() % 3 : Node::id(path.back()));
    size_type steps = rng() % 12;
    for(size_type j = 0; j < steps; j++)
    {
      if(rng() % 5 == 0 && curr > 2) { curr -= 1 + rng() % 2; } // Revisit a node.
      else { curr += 1 + rng() % 3; }
      path.push_back(Node::encode(curr, rng() % 7 == 0));
    }
    paths.push_back(path);
  }
  return paths;
}

TEST(SequenceLocateTest, EmptyPaths)
{
  std::vector<vector_type> paths { {}, { Node::encode(1, false), Node::encode(2, false) }, {} };
  GBWT index = buildGBWT(paths, true);
  SequenceLocate r_index(index, [](node_type) -> size_type { return 3; });
  std::vector<size_type> all { 0, 1, 2, 3, 4, 5 };
  EXPECT_EQ(r_index.locate(ENDMARKER, range_type(0, index.sequences() - 1)), all) << "Invalid locate() at the endmarker";
  for(size_type i = 0; i < index.sequences(); i++)
  {
    size_type expected = (Path::id(i) == 1 ? 6 : 0);
    EXPECT_EQ(r_index.sequenceLength(i), expected) << "Invalid length for sequence " << i;
    EXPECT_EQ(r_index.locateSequence(i, expected), invalid_edge()) << "Position past the end for " << i;
  }
  EXPECT_EQ(r_index.locateSequence(2, 4), edge_type(Node::encode(2, false), 1)) << "Invalid position in the non-empty path";
  EXPECT_EQ(r_index.locateSequence(3, 4), edge_type(Node::encode(1, true), 1)) << "Invalid position in the reverse path";
}

TEST(SequenceLocateTest, OnlyEmptyPaths)
{
  std::vector<vector_type> paths { {}, {}, {} };
  GBWT index = buildGBWT(paths, true);
  SequenceLocate r_index(index, [](node_type) -> size_type { return 1; });
  std::vector<size_type> all { 0, 1, 2, 3, 4, 5 };
  EXPECT_EQ(r_index.locate(ENDMARKER, range_type(0, index.sequences() - 1)), all) << "Invalid locate() at the endmarker";
  EXPECT_EQ(r_index.header.max_length, size_type(1)) << "Invalid max_length";
  for(size_type i = 0; i < index.sequences(); i++)
  {
    EXPECT_EQ(r_index.sequenceLength(i), size_type(0)) << "Invalid length for sequence " << i;
    EXPECT_EQ(r_index.locateSequence(i, 0), invalid_edge()) << "Position in an empty sequence " << i;
  }
}

TEST(SequenceLocateTest, Random)
{
  std::mt19937_64 rng(0xDEADBEEF);
  for(size_type round = 0; round < 300; round++)
  {
    std::vector<vector_type> paths = randomPaths(rng);
    GBWT index = buildGBWT(paths, true);
    std::uint64_t mult = 1 + rng() % 97;
    SequenceLocate r_index(index, [mult](node_type node) -> size_type { return 1 + (Node::id(node) * mult) % 11; });
    std::map<edge_type, size_type> truth = truePositions(index, r_index);

    for(node_type node = (index.empty() ? 0 : index.firstNode()); node < index.sigma(); node++)
    {
      if(!(index.contains(node)) || index.nodeSize(node) == 0) { continue; }
      size_type n = index.nodeSize(node);
      std::vector<size_type> sa = r_index.decompressSA(node);
      ASSERT_EQ(sa.size(), n) << "Round " << round << ": invalid SA size for node " << node;
      for(size_type i = 0; i < n; i++)
      {
        ASSERT_EQ(sa[i], truth[edge_type(node, i)]) << "Round " << round << ": invalid SA[" << i << "] for node " << node;
      }
      for(size_type i = 0; i < n; i++)
      {
        std::vector<std::pair<size_type, size_type>> result = r_index.locatePositions(SearchState(node, i, n - 1));
        ASSERT_EQ(result.size(), n - i) << "Round " << round << ": invalid result size for node " << node;
        for(size_type j = 0; j < result.size(); j++)
        {
          ASSERT_EQ(r_index.pack(result[j].first, result[j].second), sa[i + j]) << "Round " << round << ": invalid position " << (i + j) << " for node " << node;
        }
      }
    }

    // find() / extend() with every subpath.
    for(size_type i = 0; i < index.sequences(); i++)
    {
      vector_type path = index.extract(i);
      for(size_type start = 0; start < path.size(); start++)
      {
        size_type first = SequenceLocate::NO_POSITION;
        SearchState state = r_index.find(path[start], first);
        for(size_type end = start + 1; end <= path.size(); end++)
        {
          ASSERT_FALSE(state.empty()) << "Round " << round << ": empty state for sequence " << i;
          ASSERT_EQ(first, truth[edge_type(state.node, state.range.first)]) << "Round " << round << ": invalid first for sequence " << i;
          if(end < path.size()) { state = r_index.extend(state, path[end], first); }
        }
      }
    }

    // locateSequence() for every base, with the GBWT position of the visit.
    for(size_type i = 0; i < index.sequences(); i++)
    {
      std::vector<edge_type> bases = trueBases(index, r_index, i);
      std::vector<edge_type> visits;
      for(edge_type pos = index.start(i); pos.first != ENDMARKER; pos = index.LF(pos))
      {
        for(size_type j = 0; j < r_index.length(pos.first); j++) { visits.push_back(pos); }
      }
      ASSERT_EQ(bases.size(), r_index.sequenceLength(i)) << "Round " << round << ": invalid length for sequence " << i;
      for(size_type j = 0; j < bases.size(); j++)
      {
        edge_type gbwt_pos;
        ASSERT_EQ(r_index.locateSequence(i, j, gbwt_pos), bases[j]) << "Round " << round << ": invalid position " << i << ":" << j;
        ASSERT_EQ(gbwt_pos, visits[j]) << "Round " << round << ": invalid GBWT position " << i << ":" << j;
      }
    }
  }
}

//------------------------------------------------------------------------------

} // namespace
