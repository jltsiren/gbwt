/*
  Copyright (c) 2018 Jouni Siren

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

#include <set>

#include <gtest/gtest.h>

#include <gbwt/variants.h>

using namespace gbwt;

namespace
{

//------------------------------------------------------------------------------

struct VariantTestData
{
  vector_type reference;
  std::vector<range_type> sites;
  std::vector<std::vector<vector_type>> alleles;
  size_type first_sample;
  std::vector<std::vector<Phasing>> phasing_information;
  bool skip_overlaps;
  std::vector<vector_type> true_haplotypes;
};

class VariantTest : public ::testing::TestWithParam<VariantTestData>
{
public:
  VariantPaths variants;
  PhasingInformation phasings;
  size_type allele_count, num_samples;

  VariantTest()
  {
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
  }

  void SetUp() override
  {
    const VariantTestData& data = GetParam();
    variants = VariantPaths(data.reference.size());
    for(auto node : data.reference) { variants.appendToReference(node); }

    allele_count = 0;
    for(size_type site = 0; site < data.sites.size(); site++)
    {
      variants.addSite(data.sites[site].first, data.sites[site].second);
      for(const vector_type& allele : data.alleles[site])
      {
        variants.addAllele(allele);
        allele_count++;
      }
    }

    num_samples = data.phasing_information.front().size();
    phasings = PhasingInformation(data.first_sample, num_samples);
    for(const std::vector<Phasing>& site : data.phasing_information)
    {
      phasings.append(site);
    }
  }
};

// Statistics of the VariantPaths object.
TEST_P(VariantTest, VariantPathsStatistics)
{
  const VariantTestData& data = GetParam();
  EXPECT_EQ(variants.size(), data.reference.size()) << "Wrong reference length";
  EXPECT_EQ(variants.paths(), allele_count) << "Wrong allele count";
  ASSERT_EQ(variants.sites(), data.sites.size()) << "Wrong site count";
  for(size_type site = 0; site < data.sites.size(); site++)
  {
    EXPECT_EQ(variants.alleles(site), data.alleles[site].size()) << "Wrong allele count at site " << site;
    EXPECT_EQ(variants.refStart(site), data.sites[site].first) << "Wrong refStart at site " << site;
    EXPECT_EQ(variants.refEnd(site), data.sites[site].second) << "Wrong refEnd at site " << site;
  }
}

// Reference index in the VariantPaths object.
TEST_P(VariantTest, VariantPathsReferenceIndex)
{
  const VariantTestData& data = GetParam();
  variants.indexReference();
  for(size_type i = 0; i < data.reference.size(); i++)
  {
    EXPECT_EQ(variants.firstOccurrence(data.reference[i]), i) << "First occurrence of " << data.reference[i] << " not at " << i;
  }
  EXPECT_EQ(variants.firstOccurrence(0), variants.invalid_position()) << "Found occurrence of an invalid node";
}

// Opening and closing the temporary file in PhasingInformation.
TEST_P(VariantTest, PhasingInformationFile)
{
  ASSERT_TRUE(phasings.isOpen()) << "The file was not open";
  phasings.close();
  ASSERT_FALSE(phasings.isOpen()) << "The file cannot be closed";
  phasings.open();
  ASSERT_TRUE(phasings.isOpen()) << "The file cannot be opened";
}

// Statistics of the PhasingInformation object.
TEST_P(VariantTest, PhasingInformationStatistics)
{
  const VariantTestData& data = GetParam();
  EXPECT_EQ(phasings.size(), num_samples) << "Wrong sample count";
  EXPECT_EQ(phasings.offset(), data.first_sample) << "Wrong sample range start";
  EXPECT_EQ(phasings.limit(), data.first_sample + num_samples) << "Wrong sample range limit";
  EXPECT_EQ(phasings.sites(), data.phasing_information.size()) << "Wrong site count";
}

// Haplotype generation.
TEST_P(VariantTest, HaplotypeGeneration)
{
  const VariantTestData& data = GetParam();
  std::vector<vector_type> haplotypes;
  generateHaplotypes(variants, phasings,
    [](size_type) -> bool { return true; },
    [&haplotypes](const Haplotype& haplotype) { haplotypes.push_back(haplotype.path); },
    [&data](size_type, size_type) -> bool { return data.skip_overlaps; });
  ASSERT_EQ(haplotypes.size(), data.true_haplotypes.size()) << "Wrong number of generated haplotypes";
  for(size_type haplotype = 0; haplotype < haplotypes.size(); haplotype++)
  {
    EXPECT_EQ(haplotypes[haplotype], data.true_haplotypes[haplotype]) << "Wrong haplotype " << haplotype;
  }
}

// Haplotype generation using an embedded file.
TEST_P(VariantTest, HaplotypeGenerationEmbedded)
{
  // Add the file.
  variants.addFile(phasings.name(), phasings.offset(), phasings.size());
  ASSERT_EQ(variants.files(), static_cast<size_type>(1)) << "Wrong number of associated files";
  EXPECT_EQ(variants.name(0), phasings.name()) << "Wrong file name";
  EXPECT_EQ(variants.offset(0), phasings.offset()) << "Wrong sample range start";
  EXPECT_EQ(variants.count(0), phasings.size()) << "Wrong sample count";

  // Serialize and load VariantPaths.
  std::string temp_file_name = TempFile::getName("variants");
  sdsl::store_to_file(variants, temp_file_name);
  sdsl::load_from_file(variants, temp_file_name);
  TempFile::remove(temp_file_name);

  // Close and reopen the phasings file to ensure that the buffer has been written to disk.
  phasings.close();
  phasings.open();

  // Generate haplotypes from the embedded file.
  const VariantTestData& data = GetParam();
  std::vector<vector_type> haplotypes;
  std::set<std::string> empty_set;
  generateHaplotypes(variants, empty_set,
    [](size_type) -> bool { return true; },
    [&haplotypes](const Haplotype& haplotype) { haplotypes.push_back(haplotype.path); },
    [&data](size_type, size_type) -> bool { return data.skip_overlaps; });
  ASSERT_EQ(haplotypes.size(), data.true_haplotypes.size()) << "Wrong number of generated haplotypes";
  for(size_type haplotype = 0; haplotype < haplotypes.size(); haplotype++)
  {
    EXPECT_EQ(haplotypes[haplotype], data.true_haplotypes[haplotype]) << "Wrong haplotype " << haplotype;
  }
}

//------------------------------------------------------------------------------

VariantTestData basic_data
{
  { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }, // reference
  { range_type(2, 4), range_type(6, 7), range_type(8, 9) }, // sites
  { // alleles
    { { 11, 12 } },
    { { 13 }, { 14, 15 } },
    { { 16 }, {} }
  },
  10, // first_sample
  { // phasing_information: diploid-haploid-diploid, phased-unphased-phased, haploid-phased-phased
    { Phasing(1, 0, true), Phasing(1, 1, true),  Phasing(0) }, // This site has a maximal value "1|1".
    { Phasing(0),          Phasing(1, 2, false), Phasing(2, 0, true) },
    { Phasing(2, 1, true), Phasing(1, 0, true),  Phasing(0, 1, true) }
  },
  false, // skip_overlaps,
  { // true_haplotypes: (sample, phase, sequence)
    { 1, 2, 11, 12, 5, 6 },     // (0, 0, 0)
    { 1, 2, 3, 4, 5, 6 },       // (0, 1, 0)
    { 1, 2, 11, 12, 5, 6 },     // (1, 0, 0)
    { 5, 6, 13, 8 },            // (1, 0, 1)
    { 1, 2, 11, 12, 5, 6 },     // (1, 1, 0)
    { 5, 6, 14, 15, 8 },        // (1, 1, 1)
    { 1, 2, 3, 4, 5, 6 },       // (2, 0, 0)
    { 5, 6, 7, 8 },             // (0, 0, 1)
    { 8, 10 },                  // (0, 0, 2)
    { 8, 16, 10 },              // (0, 1, 1)
    { 8, 16, 10 },              // (1, 0, 2)
    { 8, 9, 10 },               // (1, 1, 2)
    { 5, 6, 14, 15, 8, 9, 10 }, // (2, 0, 1)
    { 5, 6, 7, 8, 16, 10}       // (2, 1, 0)
  }
};

INSTANTIATE_TEST_CASE_P(BasicHaplotypes, VariantTest, ::testing::Values(basic_data));

//------------------------------------------------------------------------------

VariantTestData overlap_data
{
  { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }, // reference
  { range_type(2, 4), range_type(3, 6), range_type(4, 7) }, // sites
  { // alleles
    { { 3, 11 }, { 12, 4 } },
    { { 4, 13, 6 }, { 14, 5, 15 },{ 14, 5, 6 } },
    { { 5, 16, 7 }, { 17, 6, 7 }, { 5, 6, 18 } }
  },  
  0, // first_sample
  { // phasing_information: [remove ref nodes] 0-1: 1 left, 0-1: 1 right, 1-2: only left, 1-2: only right, 1-2: 1 both, 1-2: 2 left, 1-2: 2 right
    { Phasing(2), Phasing(1), Phasing(0), Phasing(0), Phasing(0), Phasing(0), Phasing(0) },
    { Phasing(2), Phasing(1), Phasing(1), Phasing(2), Phasing(1), Phasing(3), Phasing(2) },
    { Phasing(0), Phasing(0), Phasing(2), Phasing(1), Phasing(1), Phasing(2), Phasing(3) }
  },
  false, // skip_overlaps,
  { // true_haplotypes: (sample, phase, sequence)
    { 1, 2, 3, 4, 13, 6 },                // (2, 0, 0)
    { 1, 2, 3, 14, 5, 15 },               // (3, 0, 0)
    { 1, 2, 12, 14, 5, 15, 7, 8, 9, 10 }, // (0, 0, 0)
    { 1, 2, 3, 11, 13, 6, 7, 8, 9, 10 },  // (1, 0, 0)
    { 17, 6, 7, 8, 9, 10 },               // (2, 0, 1)
    { 5, 16, 7, 8, 9, 10 },               // (3, 0, 1)
    { 1, 2, 3, 4, 13, 16, 7, 8, 9, 10 },  // (4, 0, 0)
    { 1, 2, 3, 14, 17, 6, 7, 8, 9, 10 },  // (5, 0, 0)
    { 1, 2, 3, 14, 5, 15, 18, 8, 9, 10 }  // (6, 0, 0)
  }
};

INSTANTIATE_TEST_CASE_P(OverlappingVariants, VariantTest, ::testing::Values(overlap_data));

//------------------------------------------------------------------------------

VariantTestData overlap_skip_data
{
  { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }, // reference
  { range_type(2, 4), range_type(3, 6), range_type(4, 7) }, // sites
  { // alleles
    { { 3, 11 }, { 12, 4 } },
    { { 4, 13, 6 }, { 14, 5, 15 },{ 14, 5, 6 } },
    { { 5, 16, 7 }, { 17, 6, 7 }, { 5, 6, 18 } }
  },  
  0, // first_sample
  { // phasing_information: [remove ref nodes] 0-1: 1 left, 0-1: 1 right, 1-2: only left, 1-2: only right, 1-2: 1 both, 1-2: 2 left, 1-2: 2 right
    { Phasing(2), Phasing(1), Phasing(0), Phasing(0), Phasing(0), Phasing(0), Phasing(0) },
    { Phasing(2), Phasing(1), Phasing(1), Phasing(2), Phasing(1), Phasing(3), Phasing(2) },
    { Phasing(0), Phasing(0), Phasing(2), Phasing(1), Phasing(1), Phasing(2), Phasing(3) }
  },
  true, // skip_overlaps,
  { // true_haplotypes: (sample, phase, sequence)
    { 1, 2, 12, 14, 5, 15, 7, 8, 9, 10 }, // (0, 0, 0)
    { 1, 2, 3, 11, 13, 6, 7, 8, 9, 10 },  // (1, 0, 0)
    { 1, 2, 3, 4, 13, 6, 7, 8, 9, 10 },   // (2, 0, 0)
    { 1, 2, 3, 14, 5, 15, 7, 8, 9, 10},   // (3, 0, 0)
    { 1, 2, 3, 4, 13, 16, 7, 8, 9, 10 },  // (4, 0, 0)
    { 1, 2, 3, 14, 17, 6, 7, 8, 9, 10 },  // (5, 0, 0)
    { 1, 2, 3, 14, 5, 15, 18, 8, 9, 10 }  // (6, 0, 0)
  }
};

INSTANTIATE_TEST_CASE_P(SkipOverlaps, VariantTest, ::testing::Values(overlap_skip_data));

//------------------------------------------------------------------------------

TEST(PhasingTest, ParseGenotypes)
{
  std::vector<std::vector<std::string>> genotypes =
  {
    { "1|0", "1|1", "0" },
    { "0",   "1/2", "2|0" },
    { "2|1", "1|0", "0|1" },
    { "0/0", "1/1", "0/0" }
  };
  std::vector<std::vector<Phasing>> phasing_information =
  {
    { Phasing(1, 0, true), Phasing(1, 1, true),  Phasing(0) },
    { Phasing(0),          Phasing(1, 2, false), Phasing(2, 0, true) },
    { Phasing(2, 1, true), Phasing(1, 0, true),  Phasing(0, 1, true) },
    { Phasing(0, 0, true), Phasing(1, 1, true), Phasing(0, 0, false) }
  };
  std::vector<bool> phase_homozygous = { true, true, false };
  for(size_type site = 0; site < genotypes.size(); site++)
  {
    for(size_type sample = 0; sample < genotypes[site].size(); sample++)
    {
      EXPECT_EQ(Phasing(genotypes[site][sample], true, phase_homozygous[sample]), phasing_information[site][sample])
        << "Wrong genotype for site " << site << ", sample " << sample;
    }
  }
}

TEST(PhasingTest, ForcePhasing)
{
  Phasing unphased(0, 1, false);
  unphased.forcePhased([]() { return false; });
  EXPECT_TRUE(unphased.phased) << "Forced phasing failed";
  EXPECT_EQ(unphased.first, static_cast<size_type>(0)) << "Wrong first allele";
  EXPECT_EQ(unphased.second, static_cast<size_type>(1)) << "Wrong second allele";

  Phasing swapped(0, 1, false);
  swapped.forcePhased([]() { return true; });
  EXPECT_TRUE(swapped.phased) << "Forced phasing failed";
  EXPECT_EQ(swapped.first, static_cast<size_type>(1)) << "Wrong first allele";
  EXPECT_EQ(swapped.second, static_cast<size_type>(0)) << "Wrong second allele";
}

//------------------------------------------------------------------------------

} // namespace
