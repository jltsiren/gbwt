/*
  Copyright (c) 2019, 2021, 2024, 2025 Jouni Siren

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

#include <random>
#include <sstream>

#include <gbwt/metadata.h>

using namespace gbwt;

namespace
{

//------------------------------------------------------------------------------

class MetadataTest : public ::testing::Test
{
public:
  std::vector<std::string> keys, first_keys, second_keys, third_keys;
  std::vector<PathName> paths;
  size_type path_samples, path_haplotypes, path_contigs;

  MetadataTest()
  {
    Verbosity::set(Verbosity::SILENT);
  }

  void SetUp() override
  {
    this->keys =
    {
      "first", "second", "third", "fourth", "fifth"
    };
    this->first_keys =
    {
      "first", "second", "third"
    };
    this->second_keys =
    {
      "fourth", "fifth", "sixth"
    };
    this->third_keys =
    {
      "first", "fourth", "fifth"
    };
    this->paths =
    {
      PathName(0, 0, 0, 0),
      PathName(0, 0, 1, 0),
      PathName(1, 0, 0, 0),
      PathName(1, 0, 1, 0),
      PathName(1, 1, 0, 0),
      PathName(1, 1, 1, 0),
      PathName(2, 0, 0, 0),
      PathName(2, 0, 0, 1),
      PathName(2, 0, 1, 0),
      PathName(2, 0, 1, 1)
    };
    this->path_samples = 3;
    this->path_haplotypes = 6;
    this->path_contigs = 2;
  }

  void testMergeConstructor(const Metadata& first, const Metadata& second, const Metadata& correct_result, bool same_samples, bool same_contigs, const std::string& test_name) const
  {
    std::vector<const Metadata*> sources { &first, &second };
    Metadata constructed(sources, same_samples, same_contigs);
    EXPECT_EQ(constructed, correct_result) << "Merge constructor does not work correctly " << test_name;
  }

  // first[start, limit) should be equal to second.
  bool sameSamples(const Metadata& first, const Metadata& second, size_type start, size_type limit) const
  {
    for(size_type i = start; i < limit; i++)
    {
      if(first.sample(i) != second.sample(i - start)) { return false; }
    }
    return true;
  }

  void compareSamples(const Metadata& first, const Metadata& second, bool same_samples) const
  {
    std::stringstream stream;
    stream << "(" << (first.hasSampleNames() ? "names" : "nonames") << ", "
           << (second.hasSampleNames() ? "names" : "nonames") << ", "
           << (same_samples ? "same" : "not_same") << ")";
    std::string test_name = stream.str();
    bool merge_by_names = first.hasSampleNames() & second.hasSampleNames();
    size_type correct_count = ((same_samples & !merge_by_names) ? first.samples() : first.samples() + second.samples());

    Metadata merged = first;
    merged.merge(second, same_samples, false);
    ASSERT_EQ(merged.samples(), correct_count) << "Expected " << correct_count << " samples, got " << merged.samples() << " " << test_name;
    testMergeConstructor(first, second, merged, same_samples, false, test_name);

    if(merge_by_names)
    {
      ASSERT_TRUE(merged.hasSampleNames()) << "Merged metadata object does not have sample names " << test_name;
      bool correct_samples = true;
      correct_samples &= sameSamples(merged, first, 0, first.samples());
      correct_samples &= sameSamples(merged, second, first.samples(), first.samples() + second.samples());
      EXPECT_TRUE(correct_samples) << "Sample names were not merged " << test_name;
    }
    else if(same_samples)
    {
      bool check_sample_names = false, correct_samples = true;
      if(first.hasSampleNames())
      {
        correct_samples = sameSamples(merged, first, 0, merged.samples());
        check_sample_names = true;
      }
      else if(second.hasSampleNames())
      {
        correct_samples = sameSamples(merged, second, 0, merged.samples());
        check_sample_names = true;
      }
      if(check_sample_names)
      {
        ASSERT_TRUE(merged.hasSampleNames()) << "Merged metadata object does not have sample names " << test_name;
        EXPECT_TRUE(correct_samples) << "Sample names were not merged " << test_name;
      }
      else
      {
        ASSERT_FALSE(merged.hasSampleNames()) << "Merged metadata object has sample names " << test_name;
      }
    }
    else
    {
      ASSERT_FALSE(merged.hasSampleNames()) << "Merged metadata object has sample names " << test_name;
    }
  }

  // first[start, limit) should be equal to second.
  bool sameContigs(const Metadata& first, const Metadata& second, size_type start, size_type limit) const
  {
    for(size_type i = start; i < limit; i++)
    {
      if(first.contig(i) != second.contig(i - start)) { return false; }
    }
    return true;
  }

  void compareContigs(const Metadata& first, const Metadata& second, bool same_contigs) const
  {
    std::stringstream stream;
    stream << "(" << (first.hasContigNames() ? "names" : "nonames") << ", "
           << (second.hasContigNames() ? "names" : "nonames") << ", "
           << (same_contigs ? "same" : "not_same") << ")";
    std::string test_name = stream.str();
    bool merge_by_names = first.hasContigNames() & second.hasContigNames();
    size_type correct_count = ((same_contigs & !merge_by_names) ? first.contigs() : first.contigs() + second.contigs());

    Metadata merged = first;
    merged.merge(second, false, same_contigs);
    ASSERT_EQ(merged.contigs(), correct_count) << "Expected " << correct_count << " contigs, got " << merged.contigs() << " " << test_name;
    testMergeConstructor(first, second, merged, false, same_contigs, test_name);

    if(merge_by_names)
    {
      ASSERT_TRUE(merged.hasContigNames()) << "Merged metadata object does not have contig names " << test_name;
      bool correct_contigs = true;
      correct_contigs &= sameContigs(merged, first, 0, first.contigs());
      correct_contigs &= sameContigs(merged, second, first.contigs(), first.contigs() + second.contigs());
      EXPECT_TRUE(correct_contigs) << "Contig names were not merged " << test_name;
    }
    else if(same_contigs)
    {
      bool check_contig_names = false, correct_contigs = true;
      if(first.hasContigNames())
      {
        correct_contigs = sameContigs(merged, first, 0, merged.contigs());
        check_contig_names = true;
      }
      else if(second.hasContigNames())
      {
        correct_contigs = sameContigs(merged, second, 0, merged.contigs());
        check_contig_names = true;
      }
      if(check_contig_names)
      {
        ASSERT_TRUE(merged.hasContigNames()) << "Merged metadata object does not have contig names " << test_name;
        EXPECT_TRUE(correct_contigs) << "Contig names were not merged " << test_name;
      }
      else
      {
        ASSERT_FALSE(merged.hasContigNames()) << "Merged metadata object has contig names " << test_name;
      }
    }
    else
    {
      ASSERT_FALSE(merged.hasContigNames()) << "Merged metadata object has contig names " << test_name;
    }
  }

  void comparePaths(const Metadata& first, const Metadata& second, bool same) const
  {
    std::stringstream stream;
    stream << "(" << (first.hasSampleNames() ? "names" : "nonames") << ", "
           << (second.hasSampleNames() ? "names" : "nonames") << ", "
           << (same ? "same" : "not_same") << ")";
    std::string test_name = stream.str();
    bool merge_by_names = first.hasSampleNames() & second.hasSampleNames();
    size_type correct_count = first.paths() + second.paths();

    Metadata merged = first;
    merged.merge(second, same, same);
    ASSERT_TRUE(merged.hasPathNames()) << "Merged metadata object does not have path names" << test_name;
    ASSERT_EQ(merged.paths(), correct_count) << "Expected " << correct_count << " paths, got " << merged.paths() << " " << test_name;
    testMergeConstructor(first, second, merged, same, same, test_name);

    if(merge_by_names)
    {
      bool correct_paths = true;
      for(size_type i = 0; i < first.paths(); i++)
      {
        PathName path = first.path(i);
        path.sample = merged.sample(first.sample(path.sample));
        path.contig = merged.contig(first.contig(path.contig));
        correct_paths &= (merged.path(i) == path);
      }
      for(size_type i = first.paths(); i < first.paths() + second.paths(); i++)
      {
        PathName path = second.path(i - first.paths());
        path.sample = merged.sample(second.sample(path.sample));
        path.contig = merged.contig(second.contig(path.contig));
        correct_paths &= (merged.path(i) == path);
      }
      EXPECT_TRUE(correct_paths) << "Path names were not merged correctly " << test_name;
    }
    else if(same)
    {
      bool correct_paths = true;
      for(size_type i = 0; i < first.paths(); i++)
      {
        correct_paths &= (merged.path(i) == first.path(i));
      }
      for(size_type i = first.paths(); i < first.paths() + second.paths(); i++)
      {
        correct_paths &= (merged.path(i) == second.path(i - first.paths()));
      }
      EXPECT_TRUE(correct_paths) << "Path names were not merged correctly " << test_name;
    }
    else
    {
      bool correct_paths = true;
      for(size_type i = 0; i < first.paths(); i++)
      {
        for(size_type i = 0; i < first.paths(); i++)
        {
          correct_paths &= (merged.path(i) == first.path(i));
        }
      }
      for(size_type i = first.paths(); i < first.paths() + second.paths(); i++)
      {
        PathName path = second.path(i - first.paths());
        path.sample += first.samples();
        path.contig += first.contigs();
        correct_paths &= (merged.path(i) == path);
      }
      EXPECT_TRUE(correct_paths) << "Path names were not merged correctly " << test_name;
    }
  }
};

//------------------------------------------------------------------------------

TEST_F(MetadataTest, BasicTest)
{
  // Empty metadata.
  Metadata empty;
  EXPECT_EQ(empty.samples(), static_cast<size_type>(0)) << "Empty metadata object contains samples";
  EXPECT_EQ(empty.haplotypes(), static_cast<size_type>(0)) << "Empty metadata object contains haplotypes";
  EXPECT_EQ(empty.contigs(), static_cast<size_type>(0)) << "Empty metadata object contains contigs";
  EXPECT_FALSE(empty.hasSampleNames()) << "Empty metadata object contains sample names";
  EXPECT_FALSE(empty.hasContigNames()) << "Empty metadata object contains contig names";
  EXPECT_FALSE(empty.hasPathNames()) << "Empty metadata object contains path names";

  // Set counts.
  Metadata nonempty;
  size_type samples = 1, haplotypes = 2, contigs = 3;
  nonempty.setSamples(samples);
  nonempty.setHaplotypes(haplotypes);
  nonempty.setContigs(contigs);
  EXPECT_EQ(nonempty.samples(), samples) << "Expected " << samples << " samples, got " << nonempty.samples();
  EXPECT_EQ(nonempty.haplotypes(), haplotypes) << "Expected " << haplotypes << " haplotypes, got " << nonempty.haplotypes();
  EXPECT_EQ(nonempty.contigs(), contigs) << "Expected " << contigs << " contigs, got " << nonempty.contigs();

  // Comparisons and clear().
  EXPECT_NE(empty, nonempty) << "Empty and nonempty metadata objects are equal";
  nonempty.clear();
  EXPECT_EQ(empty, nonempty) << "Cleared metadata object is not empty";
}

//------------------------------------------------------------------------------

TEST_F(MetadataTest, Samples)
{
  Metadata metadata;
  metadata.setSamples(keys);
  EXPECT_TRUE(metadata.hasSampleNames()) << "Metadata object does not contain sample names";
  ASSERT_EQ(metadata.samples(), static_cast<size_type>(keys.size())) << "Sample count is incorrect";
  bool ok = true;
  for(size_type i = 0; i < keys.size(); i++)
  {
    ok &= (metadata.sample(i) == keys[i]);
    ok &= (metadata.sample(keys[i]) == i);
  }
  EXPECT_TRUE(ok) << "Sample names are incorrect";

  {
    std::vector<std::string> first_half(keys.begin(), keys.begin() + keys.size() / 2);
    std::vector<std::string> second_half(keys.begin() + keys.size() / 2, keys.end());
    Metadata partial;
    partial.setSamples(first_half);
    partial.addSamples(second_half);
    EXPECT_EQ(metadata, partial) << "addSamples() does not work correctly";
  }

  metadata.clearSampleNames();
  EXPECT_FALSE(metadata.hasSampleNames()) << "Metadata object contains sample names";
  EXPECT_EQ(metadata.samples(), static_cast<size_type>(keys.size())) << "Clearing sample names also cleared sample count";
}

TEST_F(MetadataTest, RemoveSample)
{
  for(size_type removed_key = 0; removed_key <= keys.size(); removed_key++)
  {
    Metadata metadata;
    metadata.setSamples(keys);

    std::vector<size_type> removed = metadata.removeSample(removed_key);
    size_type correct_size = (removed_key >= keys.size() ? keys.size() : keys.size() - 1);
    EXPECT_TRUE(removed.empty()) << "Path names were removed from an object without path names";
    ASSERT_EQ(metadata.samples(), correct_size) << "Expected " << correct_size << " samples after removing key " << removed_key << ", got " << metadata.samples();
    bool ok = true;
    for(size_type i = 0; i < keys.size(); i++)
    {
      if(i < removed_key)
      {
        ok &= (metadata.sample(i) == keys[i]);
        ok &= (metadata.sample(keys[i]) == i);
      }
      else if(i == removed_key)
      {
        ok &= (metadata.sample(keys[i]) == metadata.samples());
      }
      else
      {
        ok &= (metadata.sample(i - 1) == keys[i]);
        ok &= (metadata.sample(keys[i]) == i - 1);
      }
    }
    EXPECT_TRUE(ok) << "Metadata object does not contain the correct samples after removing sample " << removed_key;
  }
}

TEST_F(MetadataTest, SampleMerging)
{
  Metadata first_names, first_nonames, second_names, second_nonames;
  first_names.setSamples(first_keys);
  first_nonames.setSamples(first_keys.size());
  second_names.setSamples(second_keys);
  second_nonames.setSamples(second_keys.size());

  compareSamples(first_nonames, second_nonames, false);
  compareSamples(first_nonames, second_nonames, true);

  compareSamples(first_names, second_nonames, false);
  compareSamples(first_names, second_nonames, true);

  compareSamples(first_nonames, second_names, false);
  compareSamples(first_nonames, second_names, true);

  compareSamples(first_names, second_names, false);
  compareSamples(first_names, second_names, true);
}

//------------------------------------------------------------------------------

TEST_F(MetadataTest, Contigs)
{
  Metadata metadata;
  metadata.setContigs(keys);
  EXPECT_TRUE(metadata.hasContigNames()) << "Metadata object does not contain contig names";
  ASSERT_EQ(metadata.contigs(), static_cast<size_type>(keys.size())) << "Contig count is incorrect";
  bool ok = true;
  for(size_type i = 0; i < keys.size(); i++)
  {
    ok &= (metadata.contig(i) == keys[i]);
    ok &= (metadata.contig(keys[i]) == i);
  }
  EXPECT_TRUE(ok) << "Contig names are incorrect";

  {
    std::vector<std::string> first_half(keys.begin(), keys.begin() + keys.size() / 2);
    std::vector<std::string> second_half(keys.begin() + keys.size() / 2, keys.end());
    Metadata partial;
    partial.setContigs(first_half);
    partial.addContigs(second_half);
    EXPECT_EQ(metadata, partial) << "addContigs() does not work correctly";
  }

  metadata.clearContigNames();
  EXPECT_FALSE(metadata.hasContigNames()) << "Metadata object contains contig names";
  EXPECT_EQ(metadata.contigs(), static_cast<size_type>(keys.size())) << "Clearing contig names also cleared contig count";
}

TEST_F(MetadataTest, RemoveContig)
{
  for(size_type removed_key = 0; removed_key <= keys.size(); removed_key++)
  {
    Metadata metadata;
    metadata.setContigs(keys);

    std::vector<size_type> removed = metadata.removeContig(removed_key);
    size_type correct_size = (removed_key >= keys.size() ? keys.size() : keys.size() - 1);
    EXPECT_TRUE(removed.empty()) << "Path names were removed from an object without path names";
    ASSERT_EQ(metadata.contigs(), correct_size) << "Expected " << correct_size << " contigs after removing key " << removed_key << ", got " << metadata.contigs();

    bool ok = true;
    for(size_type i = 0; i < keys.size(); i++)
    {
      if(i < removed_key)
      {
        ok &= (metadata.contig(i) == keys[i]);
        ok &= (metadata.contig(keys[i]) == i);
      }
      else if(i == removed_key)
      {
        ok &= (metadata.contig(keys[i]) == metadata.contigs());
      }
      else
      {
        ok &= (metadata.contig(i - 1) == keys[i]);
        ok &= (metadata.contig(keys[i]) == i - 1);
      }
    }
    EXPECT_TRUE(ok) << "Metadata object does not contain the correct contigs after removing contig " << removed_key;
  }
}

TEST_F(MetadataTest, ContigMerging)
{
  Metadata first_names, first_nonames, second_names, second_nonames;
  first_names.setContigs(first_keys);
  first_nonames.setContigs(first_keys.size());
  second_names.setContigs(second_keys);
  second_nonames.setContigs(second_keys.size());

  compareContigs(first_nonames, second_nonames, false);
  compareContigs(first_nonames, second_nonames, true);

  compareContigs(first_names, second_nonames, false);
  compareContigs(first_names, second_nonames, true);

  compareContigs(first_nonames, second_names, false);
  compareContigs(first_nonames, second_names, true);

  compareContigs(first_names, second_names, false);
  compareContigs(first_names, second_names, true);
}

//------------------------------------------------------------------------------

TEST_F(MetadataTest, Paths)
{
  Metadata metadata;
  for(const PathName& path : paths) { metadata.addPath(path); }
  EXPECT_TRUE(metadata.hasPathNames()) << "Metadata object does not contain path names";
  ASSERT_EQ(metadata.paths(), static_cast<size_type>(paths.size())) << "Path count is incorrect";
  bool ok = true;
  for(size_type i = 0; i < keys.size(); i++)
  {
    ok &= (metadata.path(i) == paths[i]);
  }
  EXPECT_TRUE(ok) << "Path names are incorrect";

  // Find paths by sample and contig.
  {
    std::vector<size_type> correct_paths;
    for(size_type i = 0; i < paths.size(); i++)
    {
      if(paths[i].sample == 1 && paths[i].contig == 0) { correct_paths.push_back(i); }
    }
    bool ok = (correct_paths == metadata.findPaths(1, 0));
    EXPECT_TRUE(ok) << "Path selection by sample and contig failed";
  }

  // Find paths by sample.
  {
    std::vector<size_type> correct_paths;
    for(size_type i = 0; i < paths.size(); i++)
    {
      if(paths[i].sample == 1) { correct_paths.push_back(i); }
    }
    bool ok = (correct_paths == metadata.pathsForSample(1));
    EXPECT_TRUE(ok) << "Path selection by sample failed";
  }

  // Find paths by contig.
  {
    std::vector<size_type> correct_paths;
    for(size_type i = 0; i < paths.size(); i++)
    {
      if(paths[i].contig == 1) { correct_paths.push_back(i); }
    }
    bool ok = (correct_paths == metadata.pathsForContig(1));
    EXPECT_TRUE(ok) << "Path selection by contig failed";
  }

  metadata.clearPathNames();
  EXPECT_FALSE(metadata.hasPathNames()) << "Metadata object contains path names";
}

TEST_F(MetadataTest, RemovePaths)
{
  Metadata metadata;
  metadata.setSamples(path_samples);
  metadata.setHaplotypes(path_haplotypes);
  metadata.setContigs(path_contigs);
  for(const PathName& path : paths) { metadata.addPath(path); }

  for(size_type sample = 0; sample <= metadata.samples(); sample++)
  {
    Metadata curr = metadata;

    size_type correct_samples = (sample >= metadata.samples() ? path_samples : path_samples - 1);
    std::vector<size_type> removed = curr.removeSample(sample);
    ASSERT_EQ(curr.samples(), correct_samples) << "Expected " << correct_samples << " samples after removing sample " << sample << ", got " << curr.samples();

    bool ok = true;
    size_type path_tail = 0, removed_tail = 0;
    for(size_type i = 0; i < paths.size(); i++)
    {
      if(paths[i].sample == sample)
      {
        if(removed_tail >= removed.size() || removed[removed_tail] != i) { ok = false; break; }
        removed_tail++;
      }
      else
      {
        PathName expected = paths[i];
        if(expected.sample > sample) { expected.sample--; }
        if(path_tail >= curr.paths() || curr.path(path_tail) != expected) { ok = false; break; }
        path_tail++;
      }
    }
    ASSERT_TRUE(ok) << "Metadata object does not contain the correct contigs after removing sample " << sample;
    EXPECT_EQ(path_tail, curr.paths()) << "Expected " << path_tail << " paths after removing sample " << sample << ", got " << curr.paths();
    EXPECT_EQ(removed_tail, removed.size()) << "Expected " << removed_tail << " removed paths after removing sample " << sample << ", got " << removed.size();
  }

  for(size_type contig = 0; contig <= metadata.contigs(); contig++)
  {
    Metadata curr = metadata;

    size_type correct_contigs = (contig >= metadata.contigs() ? path_contigs : path_contigs - 1);
    std::vector<size_type> removed = curr.removeContig(contig);
    ASSERT_EQ(curr.contigs(), correct_contigs) << "Expected " << correct_contigs << " contigs after removing contig " << contig << ", got " << curr.contigs();

    bool ok = true;
    size_type path_tail = 0, removed_tail = 0;
    for(size_type i = 0; i < paths.size(); i++)
    {
      if(paths[i].contig == contig)
      {
        if(removed_tail >= removed.size() || removed[removed_tail] != i) { ok = false; break; }
        removed_tail++;
      }
      else
      {
        PathName expected = paths[i];
        if(expected.contig > contig) { expected.contig--; }
        if(path_tail >= curr.paths() || curr.path(path_tail) != expected) { ok = false; break; }
        path_tail++;
      }
    }
    ASSERT_TRUE(ok) << "Metadata object does not contain the correct contigs after removing contig " << contig;
    EXPECT_EQ(path_tail, curr.paths()) << "Expected " << path_tail << " paths after removing contig " << contig << ", got " << curr.paths();
    EXPECT_EQ(removed_tail, removed.size()) << "Expected " << removed_tail << " removed paths after removing contig " << contig << ", got " << removed.size();
  }
}

TEST_F(MetadataTest, PathMerging)
{
  Metadata first_names, first_nonames, third_names, third_nonames;
  first_names.setSamples(first_keys);
  first_nonames.setSamples(first_keys.size());
  first_names.setContigs(first_keys);
  first_nonames.setContigs(first_keys.size());
  third_names.setSamples(third_keys);
  third_nonames.setSamples(third_keys.size());
  third_names.setContigs(third_keys);
  third_nonames.setContigs(third_keys.size());
  for(const PathName& path : paths)
  {
    first_names.addPath(path);
    first_nonames.addPath(path);
    third_names.addPath(path);
    third_nonames.addPath(path);
  }

  comparePaths(first_names, third_names, false);
  comparePaths(first_nonames, third_nonames, true);
  comparePaths(first_nonames, third_nonames, false);
}

TEST_F(MetadataTest, Serialization)
{
  Metadata original;
  original.setSamples(std::vector<std::string>(keys.begin(), keys.begin() + path_samples));
  original.setHaplotypes(path_haplotypes);
  original.setContigs(std::vector<std::string>(keys.begin(), keys.begin() + path_contigs));
  for(const PathName& path : paths) { original.addPath(path); }
  size_t expected_size = original.simple_sds_size() * sizeof(sdsl::simple_sds::element_type);

  std::string sdsl_filename = TempFile::getName("Metadata");
  sdsl::store_to_file(original, sdsl_filename);
  Metadata sdsl_copy; sdsl::load_from_file(sdsl_copy, sdsl_filename);
  TempFile::remove(sdsl_filename);
  EXPECT_EQ(original, sdsl_copy) << "SDSL serialization failed";

  std::string simple_sds_filename = TempFile::getName("Metadata");
  sdsl::simple_sds::serialize_to(original, simple_sds_filename);
  std::ifstream in(simple_sds_filename, std::ios_base::binary);
  size_t bytes = fileSize(in);
  ASSERT_EQ(bytes, expected_size) << "Invalid Simple-SDS file size";
  Metadata simple_sds_copy; simple_sds_copy.simple_sds_load(in);
  in.close();
  TempFile::remove(simple_sds_filename);
  EXPECT_EQ(original, simple_sds_copy) << "Simple-SDS serialization failed";
}

//------------------------------------------------------------------------------

class PathNameTest : public ::testing::Test
{
public:
  Metadata metadata;

  PathNameTest()
  {
    Verbosity::set(Verbosity::SILENT);
  }

  void SetUp() override
  {
    std::vector<std::string> samples { "sample1", "sample2" };
    this->metadata.setSamples(samples);
    std::vector<std::string> contigs { "contig1", "contig2" };
    this->metadata.setContigs(contigs);
    this->metadata.setHaplotypes(4);

    // Sample 1: Not fragmented
    this->metadata.addPath(0, 0, 0, 0);
    this->metadata.addPath(0, 0, 1, 0);
    this->metadata.addPath(0, 1, 0, 0);
    this->metadata.addPath(0, 1, 1, 0);

    // Sample 2: Fragmented
    this->metadata.addPath(1, 0, 0, 0);
    this->metadata.addPath(1, 0, 0, 1000);
    this->metadata.addPath(1, 0, 0, 2000);
    this->metadata.addPath(1, 0, 1, 0);
    this->metadata.addPath(1, 0, 1, 800);
    this->metadata.addPath(1, 0, 1, 2100);
    this->metadata.addPath(1, 1, 0, 0);
    this->metadata.addPath(1, 1, 0, 1400);
    this->metadata.addPath(1, 1, 1, 0);
    this->metadata.addPath(1, 1, 1, 1500);
    this->metadata.addPath(1, 1, 1, 2500);
  }

  void findFragment(const FullPathName& name, size_type truth_id) const
  {
    size_type fragment_id = this->metadata.findFragment(name);
    ASSERT_EQ(fragment_id, truth_id) << "Wrong fragment for sample " << name.sample_name <<
      ", contig " << name.contig_name << ", haplotype " << name.haplotype << ", offset " << name.offset;
  }
};

TEST_F(PathNameTest, NameConversions)
{
  for(size_type path_id = 0; path_id < this->metadata.paths(); path_id++)
  {
    PathName path = this->metadata.path(path_id);
    FullPathName full_path = this->metadata.fullPath(path_id);
    PathName::path_name_type sample_id = this->metadata.sample(full_path.sample_name);
    ASSERT_EQ(sample_id, path.sample) << "Wrong sample id for FullPathName " << path_id;
    PathName::path_name_type contig_id = this->metadata.contig(full_path.contig_name);
    ASSERT_EQ(contig_id, path.contig) << "Wrong contig id for FullPathName " << path_id;

    PathName converted = this->metadata.path(full_path);
    ASSERT_EQ(converted, path) << "FullPathName to PathName conversion failed for path " << path_id;
  }

  {
    FullPathName wrong_sample { "sample3", "contig1", 1, 0 };
    PathName converted = this->metadata.path(wrong_sample);
    PathName truth(2, 0, 1, 0);
    ASSERT_EQ(converted, truth) << "FullPathName to PathName conversion failed with a non-existing sample";
  }

  {
    FullPathName wrong_contig { "sample1", "contig3", 0, 0 };
    PathName converted = this->metadata.path(wrong_contig);
    PathName truth(0, 2, 0, 0);
    ASSERT_EQ(converted, truth) << "FullPathName to PathName conversion failed with a non-existing contig";
  }

  {
    FullPathName wrong_haplotype { "sample1", "contig1", 4, 0 };
    PathName converted = this->metadata.path(wrong_haplotype);
    PathName truth(0, 0, 4, 0);
    ASSERT_EQ(converted, truth) << "FullPathName to PathName conversion failed with a non-existing haplotype";
  }

  {
    FullPathName wrong_offset { "sample2", "contig2", 0, 3000 };
    PathName converted = this->metadata.path(wrong_offset);
    PathName truth(1, 1, 0, 3000);
    ASSERT_EQ(converted, truth) << "FullPathName to PathName conversion failed with a non-existing offset";
  }
}

TEST_F(PathNameTest, FindFragment)
{
  // All actual paths with offset - 1, offset, offset + 1
  for(size_type path_id = 0; path_id < this->metadata.paths(); path_id++)
  {
    FullPathName name = this->metadata.fullPath(path_id);
    if(name.offset > 0)
    {
      // In our test data, path fragments appear in order.
      FullPathName prev = name; prev.offset--;
      this->findFragment(prev, path_id - 1);
    }
    this->findFragment(name, path_id);
    name.offset++;
    // In our test data, there are no fragments starting at successive offsets.
    this->findFragment(name, path_id);
  }

  {
    FullPathName wrong_sample { "sample3", "contig1", 1, 0 };
    ASSERT_EQ(this->metadata.findFragment(wrong_sample), this->metadata.paths()) << "Found a fragment for a non-existing sample";
  }

  {
    FullPathName wrong_contig { "sample1", "contig3", 0, 0 };
    ASSERT_EQ(this->metadata.findFragment(wrong_contig), this->metadata.paths()) << "Found a fragment for a non-existing contig";
  }

  {
    FullPathName wrong_haplotype { "sample1", "contig1", 4, 0 };
    ASSERT_EQ(this->metadata.findFragment(wrong_haplotype), this->metadata.paths()) << "Found a fragment for a non-existing haplotype";
  }

  size_type old_paths = this->metadata.paths();
  size_type sample_id = this->metadata.samples();
  this->metadata.addSamples({ "sample3" });
  this->metadata.setHaplotypes(this->metadata.haplotypes() + 1);
  this->metadata.addPath(sample_id, 0, 0, 1000);
  this->metadata.addPath(sample_id, 0, 0, 501);
  this->metadata.addPath(sample_id, 0, 0, 500);

  {
    // Fragments are not in order.
    FullPathName unordered { "sample3", "contig1", 0, 800 };
    this->findFragment(unordered, old_paths + 1);
  }

  {
    // Before the first fragment.
    FullPathName before { "sample3", "contig1", 0, 499 };
    this->findFragment(before, this->metadata.paths());
  }

  {
    // First fragment.
    FullPathName first { "sample3", "contig1", 0, 500 };
    this->findFragment(first, old_paths + 2);
  }

  {
    // Fragment starting immediately after the first fragment.
    FullPathName after { "sample3", "contig1", 0, 501 };
    this->findFragment(after, old_paths + 1);
  }
}

//------------------------------------------------------------------------------

class FragmentMapTest : public ::testing::Test
{
public:
  Metadata metadata;
  std::vector<PathName> paths;
  std::vector<std::vector<FragmentMap::Fragment>> chains;

  void SetUp() override
  {
    this->addChain(0, 0, 0, 2);
    this->addChain(0, 0, 1, 1);
    this->addChain(0, 1, 0, 3);
    this->addChain(0, 1, 1, 2);
    this->addChain(1, 0, 0, 3);
    this->addChain(1, 0, 1, 3);
    this->addChain(1, 1, 0, 1);
    this->addChain(1, 1, 1, 4);
  }

  void addChain(size_type sample, size_type contig, size_type haplotype, size_type n)
  {
    std::vector<FragmentMap::Fragment> chain;
    for(size_type i = 0; i < n; i++)
    {
      this->paths.push_back(PathName(sample, contig, haplotype, i));
      size_type path = this->paths.size() - 1;
      size_type sequence = this->chains.size();
      size_type prev = (i == 0 ? invalid_sequence() : this->paths.size() - 2);
      size_type next = (i == n - 1 ? invalid_sequence() : this->paths.size());
      chain.push_back({ path, sequence, prev, next });
    }
    this->chains.push_back(chain);
  }

  void shufflePaths(std::uint64_t seed)
  {
    // Shuffle the paths.
    std::vector<std::pair<PathName, size_type>> pairs;
    for(size_type i = 0; i < this->paths.size(); i++) { pairs.push_back({ this->paths[i], i }); }
    std::shuffle(pairs.begin(), pairs.end(), std::mt19937_64(seed));

    // Determine path id mapping and update the paths.
    std::vector<size_type> old_to_new_path_id(this->paths.size());
    for(size_type i = 0; i < pairs.size(); i++)
    {
      this->paths[i] = pairs[i].first;
      old_to_new_path_id[pairs[i].second] = i;
    }

    // Update path ids in the chains.
    for(size_type i = 0; i < this->chains.size(); i++)
    {
      for(FragmentMap::Fragment& fragment : this->chains[i])
      {
        fragment.path = old_to_new_path_id[fragment.path];
        fragment.prev = (fragment.prev == invalid_sequence() ? invalid_sequence() : old_to_new_path_id[fragment.prev]);
        fragment.next = (fragment.next == invalid_sequence() ? invalid_sequence() : old_to_new_path_id[fragment.next]);
      }
    }
  }

  void initMetadata()
  {
    for(const PathName& path : this->paths) { this->metadata.addPath(path); }
  }

  void checkChains(const FragmentMap& fragments) const
  {
    ASSERT_EQ(fragments.size(), this->chains.size()) << "FragmentMap size is incorrect";
    for(size_type i = 0; i < this->chains.size(); i++)
    {
      const std::vector<FragmentMap::Fragment>& chain = this->chains[i];
      for(size_type j = 0; j < chain.size(); j++)
      {
        const FragmentMap::Fragment& fragment = chain[j];
        size_type forward = Path::encode(fragment.path, false);
        size_type next_fw = (fragment.next == invalid_sequence() ? invalid_sequence() : Path::encode(fragment.next, false));
        size_type prev_fw = (fragment.prev == invalid_sequence() ? invalid_sequence() : Path::encode(fragment.prev, false));
        size_type reverse = Path::encode(fragment.path, true);
        size_type next_rv = (fragment.prev == invalid_sequence() ? invalid_sequence() : Path::encode(fragment.prev, true));
        size_type prev_rv = (fragment.next == invalid_sequence() ? invalid_sequence() : Path::encode(fragment.next, true));
        ASSERT_EQ(fragments.next(fragment.path, false), fragment.next) << "Wrong next fragment id for path " << fragment.path << " (forward)";
        ASSERT_EQ(fragments.oriented_next(forward), next_fw) << "Wrong oriented next fragment id for path " << fragment.path << " (forward)";
        ASSERT_EQ(fragments.next(fragment.path, true), fragment.prev) << "Wrong next fragment id for path " << fragment.path << " (reverse)";
        ASSERT_EQ(fragments.oriented_next(reverse), next_rv) << "Wrong oriented next fragment id for path " << fragment.path << " (reverse)";
        ASSERT_EQ(fragments.prev(fragment.path, false), fragment.prev) << "Wrong previous fragment id for path " << fragment.path << " (forward)";
        ASSERT_EQ(fragments.oriented_prev(forward), prev_fw) << "Wrong oriented previous fragment id for path " << fragment.path << " (forward)";
        ASSERT_EQ(fragments.prev(fragment.path, true), fragment.next) << "Wrong previous fragment id for path " << fragment.path << " (reverse)";
        ASSERT_EQ(fragments.oriented_prev(reverse), prev_rv) << "Wrong oriented previous fragment id for path " << fragment.path << " (reverse)";
        ASSERT_EQ(fragments.chain(fragment.path), fragment.chain) << "Wrong chain id for path " << fragment.path;
      }
    }
  }

  void queryWithInvalidId(const FragmentMap& fragments) const
  {
    size_type invalid_path = this->metadata.paths();
    size_type invalid_gbwt_sequence = Path::encode(invalid_path, false);
    ASSERT_EQ(fragments.next(invalid_path), invalid_sequence()) << "Found a next fragment with an invalid path id";
    ASSERT_EQ(fragments.oriented_next(invalid_gbwt_sequence), invalid_sequence()) << "Found an oriented next fragment with an invalid path id";
    ASSERT_EQ(fragments.prev(invalid_path), invalid_sequence()) << "Found a previous fragment with an invalid path id";
    ASSERT_EQ(fragments.oriented_prev(invalid_gbwt_sequence), invalid_sequence()) << "Found an oriented previous fragment with an invalid path id";
    ASSERT_EQ(fragments.chain(invalid_path), invalid_sequence()) << "Found a chain id with an invalid path id";
  }
};

TEST_F(FragmentMapTest, Empty)
{
  FragmentMap fragments(this->metadata, false);
  ASSERT_TRUE(fragments.empty()) << "Empty FragmentMap is not empty";
  this->queryWithInvalidId(fragments);
}

TEST_F(FragmentMapTest, Ordered)
{
  this->initMetadata();
  FragmentMap fragments(this->metadata, false);
  this->checkChains(fragments);
  this->queryWithInvalidId(fragments);
}

TEST_F(FragmentMapTest, Shuffled)
{
  this->shufflePaths(0x0123456789ABCDEF);
  this->initMetadata();
  FragmentMap fragments(this->metadata, false);
  this->checkChains(fragments);
  this->queryWithInvalidId(fragments);
}

//------------------------------------------------------------------------------

} // namespace
