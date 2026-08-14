#include <gtest/gtest.h>

#include <set>
#include <unordered_map>
#include <unordered_set>

#include <gbwt/dynamic_gbwt.h>

using namespace gbwt;

namespace
{

//------------------------------------------------------------------------------

// This is the components_ref test case from GBWTGraph, except we have added a
// gap (node 24) into the node id space of the second component.

vector_type
createPath(const std::vector<std::pair<node_type, bool>>& nodes)
{
  vector_type result;
  for(const std::pair<node_type, bool>& node_visit : nodes)
  {
    result.push_back(Node::encode(node_visit.first, node_visit.second));
  }
  return result;
}

std::vector<vector_type>
createPaths(bool first_component, bool second_component)
{
  std::vector<vector_type> result;
  if(first_component)
  {
    result.push_back(createPath({ { 11, false }, { 13, false }, { 14, false}, { 15, false }, { 17, false } }));
    result.push_back(createPath({ { 11, false }, { 12, false }, { 14, false}, { 15, false }, { 17, false } }));
    result.push_back(createPath({ { 11, false }, { 13, false }, { 14, false}, { 16, false }, { 17, false } }));
  }
  if(second_component)
  {
    result.push_back(createPath({ { 21, false }, { 23, false }, { 25, true }, { 22, true }, { 21, true } }));
    result.push_back(createPath({ { 21, false }, { 22, false }, { 25, false }, { 23, true }, { 21, true } }));
    result.push_back(createPath({ { 21, false }, { 22, false }, { 25, false }, { 26, false } }));
  }
  return result;
}

std::vector<FullPathName>
createPathNames(bool first_component, bool second_component)
{
  std::vector<FullPathName> result;
  if(first_component)
  {
    result.push_back(FullPathName { "ref", "A", 0, 0 });
    result.push_back(FullPathName { "sample", "A", 1, 0 });
    result.push_back(FullPathName { "sample", "A", 2, 0 });
  }
  if(second_component)
  {
    result.push_back(FullPathName { "ref", "B", 0, 0 });
    result.push_back(FullPathName { "sample", "B", 1, 0 });
    result.push_back(FullPathName { "sample", "B", 2, 0 });
  }
  return result;
}

template <typename IndexType>
IndexType buildGBWT(const std::vector<vector_type>& paths, const std::vector<FullPathName>& path_names, bool bidirectional)
{
  size_type total_length = 0;
  for(const vector_type& path : paths) { total_length += path.size(); }
  total_length += paths.size(); // Endmarkers.
  if(bidirectional) { total_length *= 2; }

  // 32 bits should be enough, because we always compile with GBWT_SAVE_MEMORY.
  // We sample every other position to be able to test the edge cases with merging/splitting DASamples.
  GBWTBuilder builder(32, total_length, 2);
  for(const vector_type& path : paths) { builder.insert(path, bidirectional); }
  builder.finish();
  IndexType index(builder.index);
  index.addMetadata();

  // TODO: We should integrate metadata construction in the builder.
  std::unordered_map<std::string, size_type> sample_mapping, contig_mapping;
  std::vector<std::string> sample_names, contig_names;
  std::set<std::pair<size_type, size_type>> haplotypes;
  for(const FullPathName& name : path_names)
  {
    auto sample_result = sample_mapping.try_emplace(name.sample_name, sample_mapping.size());
    if(sample_result.second) { sample_names.push_back(name.sample_name); }
    auto contig_result = contig_mapping.try_emplace(name.contig_name, contig_mapping.size());
    if(contig_result.second) { contig_names.push_back(name.contig_name); }
    haplotypes.insert(std::make_pair(name.haplotype, name.offset));
  }
  index.metadata.setSamples(sample_names);
  index.metadata.setContigs(contig_names);
  index.metadata.setHaplotypes(haplotypes.size());
  for(size_type i = 0; i < path_names.size(); i++)
  {
    const FullPathName& name = path_names[i];
    size_type sample_id = sample_mapping[name.sample_name];
    size_type contig_id = contig_mapping[name.contig_name];
    index.metadata.addPath(sample_id, contig_id, name.haplotype, name.offset);
  }

  return index;
}

//------------------------------------------------------------------------------

void compareGBWTs(const GBWT& index, const GBWT& truth, const std::string& test_case_name)
{
  std::string test_case = (test_case_name.empty() ? std::string() : std::string(" in ") + test_case_name);

  bool same_header = (index.header == truth.header);
  if(!same_header)
  {
    EXPECT_EQ(index.header.sequences, truth.header.sequences) << "Wrong number of sequences in GBWT header" << test_case;
    EXPECT_EQ(index.header.size, truth.header.size) << "Wrong size in GBWT header" << test_case;
    EXPECT_EQ(index.header.offset, truth.header.offset) << "Wrong offset in GBWT header" << test_case;
    EXPECT_EQ(index.header.alphabet_size, truth.header.alphabet_size) << "Wrong alphabet size in GBWT header" << test_case;
    EXPECT_EQ(index.header.flags, truth.header.flags) << "Wrong flags in GBWT header" << test_case;
  }
  ASSERT_TRUE(same_header) << "Wrong GBWT header" << test_case;

  bool same_tags = (index.tags == truth.tags);
  ASSERT_TRUE(same_tags) << "Wrong GBWT tags" << test_case;

  bool same_recordarray_index = (index.bwt.index == truth.bwt.index);
  ASSERT_TRUE(same_recordarray_index) << "Wrong index for GBWT records" << test_case;
  bool same_recordarray_data = (index.bwt.data == truth.bwt.data);
  ASSERT_TRUE(same_recordarray_data) << "Wrong data for GBWT records" << test_case;

  bool same_dasamples_sampled = (index.da_samples.sampled_records == truth.da_samples.sampled_records);
  ASSERT_TRUE(same_dasamples_sampled) << "Wrong sampled records in DA samples" << test_case;
  bool same_dasamples_ranges = (index.da_samples.bwt_ranges == truth.da_samples.bwt_ranges);
  ASSERT_TRUE(same_dasamples_ranges) << "Wrong BWT ranges in DA samples" << test_case;
  bool same_dasamples_offsets = (index.da_samples.sampled_offsets == truth.da_samples.sampled_offsets);
  ASSERT_TRUE(same_dasamples_offsets) << "Wrong sampled offsets in DA samples" << test_case;
  bool same_dasamples_samples = (index.da_samples.array == truth.da_samples.array);
  ASSERT_TRUE(same_dasamples_samples) << "Wrong samples in DA samples" << test_case;

  bool same_metadata = (index.metadata == truth.metadata);
  if(!same_metadata)
  {
    EXPECT_EQ(index.metadata.samples(), truth.metadata.samples()) << "Wrong number of samples in GBWT metadata" << test_case;
    EXPECT_EQ(index.metadata.haplotypes(), truth.metadata.haplotypes()) << "Wrong number of haplotypes in GBWT metadata" << test_case;
    EXPECT_EQ(index.metadata.contigs(), truth.metadata.contigs()) << "Wrong number of contigs in GBWT metadata" << test_case;
    EXPECT_EQ(index.metadata.paths(), truth.metadata.paths()) << "Wrong number of paths in GBWT metadata" << test_case;
    EXPECT_EQ(index.metadata.hasPathNames(), truth.metadata.hasPathNames()) << "Wrong path name flag in GBWT metadata" << test_case;
    EXPECT_EQ(index.metadata.hasSampleNames(), truth.metadata.hasSampleNames()) << "Wrong sample name flag in GBWT metadata" << test_case;
    EXPECT_EQ(index.metadata.hasContigNames(), truth.metadata.hasContigNames()) << "Wrong contig name flag in GBWT metadata" << test_case;
  }
  ASSERT_TRUE(same_metadata) << "Wrong GBWT metadata" << test_case;
}

//------------------------------------------------------------------------------

using IndexTypes = ::testing::Types<GBWT, DynamicGBWT>;

template <typename T>
class SerializeTest : public ::testing::Test
{
public:
  std::vector<vector_type> paths;
  std::vector<FullPathName> path_names;

  void SetUp() override
  {
    this->paths = createPaths(true, true);
    this->path_names = createPathNames(true, true);
  }
};

TYPED_TEST_SUITE(SerializeTest, IndexTypes);

TYPED_TEST(SerializeTest, SDSL)
{
  using index_type = TypeParam;

  index_type index = buildGBWT<index_type>(this->paths, this->path_names, true);
  std::string filename = TempFile::getName("gbwt");
  sdsl::store_to_file(index, filename);

  index_type loaded_index;
  sdsl::load_from_file(loaded_index, filename);
  TempFile::remove(filename);

  // Convert to GBWT, as we can only compare them.
  GBWT loaded(loaded_index);
  GBWT truth(index);
  compareGBWTs(loaded, truth, "SDSL");
}

TYPED_TEST(SerializeTest, SimpleSDS)
{
  using index_type = TypeParam;

  index_type index = buildGBWT<index_type>(this->paths, this->path_names, true);
  std::string filename = TempFile::getName("gbwt");
  sdsl::simple_sds::serialize_to(index, filename);

  index_type loaded_index;
  sdsl::simple_sds::load_from(loaded_index, filename);
  TempFile::remove(filename);

  // Convert to GBWT, as we can only compare them.
  GBWT loaded(loaded_index);
  GBWT truth(index);
  compareGBWTs(loaded, truth, "Simple-SDS");
}

TYPED_TEST(SerializeTest, SimpleSDSVersion)
{
  using index_type = TypeParam;

  index_type index = buildGBWT<index_type>(this->paths, this->path_names, true);
  GBWT truth(index);

  for(std::uint32_t version = GBWTHeader::MIN_SERIALIZE_VERSION; version <= GBWTHeader::VERSION; version++)
  {
    std::string filename = TempFile::getName("gbwt");
    std::ofstream out(filename, std::ios::binary);
    index.simple_sds_serialize_version(out, version);
    out.close();

    index_type loaded_index;
    sdsl::simple_sds::load_from(loaded_index, filename);
    TempFile::remove(filename);

    // Convert to GBWT, as we can only compare them.
    GBWT loaded(loaded_index);
    compareGBWTs(loaded, truth, "Simple-SDS version " + std::to_string(version));
  }
}

//------------------------------------------------------------------------------

class MergeSplitTest : public ::testing::Test
{
public:
  std::vector<vector_type> first_paths, second_paths;
  std::vector<FullPathName> first_path_names, second_path_names;
  std::unordered_set<size_type> first_ids, second_ids;

  void SetUp() override
  {
    this->first_paths = createPaths(true, false);
    this->second_paths = createPaths(false, true);

    this->first_path_names = createPathNames(true, false);
    this->second_path_names = createPathNames(false, true);

    this->first_ids = { 11, 12, 13, 14, 15, 16, 17 };
    this->second_ids = { 21, 22, 23, 25, 26 };
  }
};

TEST_F(MergeSplitTest, Merge)
{
  for(const auto& test_case : { std::make_pair(false, "unidirectional"), std::make_pair(true, "bidirectional") })
  {
    bool bidirectional = test_case.first;
    const std::string& test_case_name = test_case.second;

    GBWT first = buildGBWT<GBWT>(this->first_paths, this->first_path_names, bidirectional);
    GBWT second = buildGBWT<GBWT>(this->second_paths, this->second_path_names, bidirectional);
    GBWT merged({ first, second });

    std::vector<vector_type> merged_paths = this->first_paths;
    merged_paths.insert(merged_paths.end(), this->second_paths.begin(), this->second_paths.end());
    std::vector<FullPathName> merged_path_names = this->first_path_names;
    merged_path_names.insert(merged_path_names.end(), this->second_path_names.begin(), this->second_path_names.end());
    GBWT truth = buildGBWT<GBWT>(merged_paths, merged_path_names, bidirectional);

    compareGBWTs(merged, truth, test_case_name);
  }
}

TEST_F(MergeSplitTest, Split)
{
  for(const auto& test_case : { std::make_pair(false, "unidirectional"), std::make_pair(true, "bidirectional") })
  {
    bool bidirectional = test_case.first;
    const std::string& test_case_name = test_case.second;

    std::vector<vector_type> merged_paths = this->first_paths;
    merged_paths.insert(merged_paths.end(), this->second_paths.begin(), this->second_paths.end());
    std::vector<FullPathName> merged_path_names = this->first_path_names;
    merged_path_names.insert(merged_path_names.end(), this->second_path_names.begin(), this->second_path_names.end());
    GBWT merged = buildGBWT<GBWT>(merged_paths, merged_path_names, bidirectional);

    size_t subgraphs = 2;
    std::vector<GBWT> parts = merged.split(subgraphs, [&](node_type node) -> size_type
    {
      size_type id = Node::id(node);
      if(this->first_ids.find(id) != this->first_ids.end()) { return 0; }
      if(this->second_ids.find(id) != this->second_ids.end()) { return 1; }
      return subgraphs;
    });
    ASSERT_EQ(parts.size(), subgraphs) << "Wrong number of subgraph indexes in " << test_case_name;

    GBWT first_truth = buildGBWT<GBWT>(this->first_paths, this->first_path_names, bidirectional);
    GBWT second_truth = buildGBWT<GBWT>(this->second_paths, this->second_path_names, bidirectional);

    std::string first_name(test_case_name); first_name += " (first)";
    compareGBWTs(parts[0], first_truth, first_name);
    std::string second_name(test_case_name); second_name += " (second)";
    compareGBWTs(parts[1], second_truth, second_name);
  }
}

TEST_F(MergeSplitTest, SplitSingle)
{
  for(const auto& test_case : { std::make_pair(false, "unidirectional"), std::make_pair(true, "bidirectional") })
  {
    bool bidirectional = test_case.first;
    const std::string& test_case_name = test_case.second;
    size_t subgraphs = 1;
    auto node_to_subgraph = [&](node_type node) -> size_type
    {
      size_type id = Node::id(node);
      if(this->first_ids.find(id) != this->first_ids.end()) { return 0; }
      if(this->second_ids.find(id) != this->second_ids.end()) { return 1; }
      return 2;
    };

    std::vector<vector_type> merged_paths = this->first_paths;
    merged_paths.insert(merged_paths.end(), this->second_paths.begin(), this->second_paths.end());
    std::vector<FullPathName> merged_path_names = this->first_path_names;
    merged_path_names.insert(merged_path_names.end(), this->second_path_names.begin(), this->second_path_names.end());
    GBWT merged = buildGBWT<GBWT>(merged_paths, merged_path_names, bidirectional);

    std::vector<GBWT> truths;
    truths.push_back(buildGBWT<GBWT>(this->first_paths, this->first_path_names, bidirectional));
    truths.push_back(buildGBWT<GBWT>(this->second_paths, this->second_path_names, bidirectional));

    for(node_type part : { 0, 1 })
    {
      std::string name(test_case_name); name += " (part " + std::to_string(part) + ")";
      std::vector<GBWT> parts = merged.split(subgraphs, [&](node_type node) -> size_type
      {
        return (node_to_subgraph(node) == part ? 0 : 1);
      });
      ASSERT_EQ(parts.size(), subgraphs) << "Wrong number of subgraph indexes in " << name;
      compareGBWTs(parts[0], truths[part], name);
    }
  }
}

//------------------------------------------------------------------------------

} // namespace
