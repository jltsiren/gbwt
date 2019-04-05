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
      { static_cast<size_type>(0), static_cast<size_type>(0), static_cast<size_type>(0), static_cast<size_type>(0) },
      { static_cast<size_type>(0), static_cast<size_type>(0), static_cast<size_type>(1), static_cast<size_type>(0) },
      { static_cast<size_type>(1), static_cast<size_type>(0), static_cast<size_type>(0), static_cast<size_type>(0) },
      { static_cast<size_type>(1), static_cast<size_type>(0), static_cast<size_type>(1), static_cast<size_type>(0) },
      { static_cast<size_type>(1), static_cast<size_type>(1), static_cast<size_type>(0), static_cast<size_type>(0) },
      { static_cast<size_type>(1), static_cast<size_type>(1), static_cast<size_type>(1), static_cast<size_type>(0) },
      { static_cast<size_type>(2), static_cast<size_type>(0), static_cast<size_type>(0), static_cast<size_type>(0) },
      { static_cast<size_type>(2), static_cast<size_type>(0), static_cast<size_type>(0), static_cast<size_type>(1) },
      { static_cast<size_type>(2), static_cast<size_type>(0), static_cast<size_type>(1), static_cast<size_type>(0) },
      { static_cast<size_type>(2), static_cast<size_type>(0), static_cast<size_type>(1), static_cast<size_type>(1) }
    };
    this->path_samples = 3;
    this->path_haplotypes = 6;
    this->path_contigs = 2;
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
    ASSERT_TRUE(merged.check()) << "Merged metadata object is not in a valid state " << test_name;
    ASSERT_EQ(merged.samples(), correct_count) << "Expected " << correct_count << " samples, got " << merged.samples() << " " << test_name;

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
    ASSERT_TRUE(merged.check()) << "Merged metadata object is not in a valid state " << test_name;
    ASSERT_EQ(merged.contigs(), correct_count) << "Expected " << correct_count << " contigs, got " << merged.contigs() << " " << test_name;

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
    ASSERT_TRUE(merged.check()) << "Merged metadata object is not in a valid state " << test_name;
    ASSERT_TRUE(merged.hasPathNames()) << "Merged metadata object does not have path names" << test_name;
    ASSERT_EQ(merged.paths(), correct_count) << "Expected " << correct_count << " paths, got " << merged.paths() << " " << test_name;

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
  ASSERT_TRUE(empty.check()) << "Empty metadata object is not in a valid state";
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
  ASSERT_TRUE(nonempty.check()) << "Metadata object is not in a valid state";
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
  ASSERT_TRUE(metadata.check()) << "Metadata object with sample names is not in a valid state";
  EXPECT_TRUE(metadata.hasSampleNames()) << "Metadata object does not contain sample names";
  ASSERT_EQ(metadata.samples(), static_cast<size_type>(keys.size())) << "Sample count is incorrect";
  bool ok = true;
  for(size_type i = 0; i < keys.size(); i++)
  {
    ok &= (metadata.sample(i) == keys[i]);
    ok &= (metadata.sample(keys[i]) == i);
  }
  EXPECT_TRUE(ok) << "Sample names are incorrect";

  metadata.clearSampleNames();
  ASSERT_TRUE(metadata.check()) << "Metadata object is not in a valid state after clearing sample names";
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
  ASSERT_TRUE(metadata.check()) << "Metadata object with contig names is not in a valid state";
  EXPECT_TRUE(metadata.hasContigNames()) << "Metadata object does not contain contig names";
  ASSERT_EQ(metadata.contigs(), static_cast<size_type>(keys.size())) << "Contig count is incorrect";
  bool ok = true;
  for(size_type i = 0; i < keys.size(); i++)
  {
    ok &= (metadata.contig(i) == keys[i]);
    ok &= (metadata.contig(keys[i]) == i);
  }
  EXPECT_TRUE(ok) << "Contig names are incorrect";

  metadata.clearContigNames();
  ASSERT_TRUE(metadata.check()) << "Metadata object is not in a valid state after clearing contig names";
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
  ASSERT_TRUE(metadata.check()) << "Metadata object with path names is not in a valid state";
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
  ASSERT_TRUE(metadata.check()) << "Metadata object is not in a valid state after clearing path names";
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

  std::string filename = TempFile::getName("Metadata");
  sdsl::store_to_file(original, filename);
  Metadata copy;
  sdsl::load_from_file(copy, filename);
  TempFile::remove(filename);

  EXPECT_EQ(original, copy) << "Metadata serialization failed";
}

//------------------------------------------------------------------------------

} // namespace
