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


#ifndef GBWT_VARIANTS_H
#define GBWT_VARIANTS_H

#include <gbwt/utils.h>

namespace gbwt
{

/*
  variants.h: Support structures for transforming phased VCF files into paths.
*/

//------------------------------------------------------------------------------

struct Haplotype
{
  std::vector<node_type> path;
  size_type              offset; // In the reference.

  Haplotype() : offset(0) {}
  Haplotype(size_type reference_offset) : offset(reference_offset) {}

  size_type size() const { return this->path.size(); }
  bool empty() const { return this->path.empty(); }
};

//------------------------------------------------------------------------------

/*
  Paths corrensponding to the reference and the variants. Allele i corresponds to
  alternate allele i+1 in the VCF file.
*/

struct VariantPaths
{
  /*
    reference     reference path
    ref_starts    starting positions of each site in the reference
    ref_ends      the positions following each site in the reference

    alt_paths     concatenated allele paths
    path_starts   pointers to the start of each path in alt_paths with a sentinel at the end
    site_starts   pointers to the start of each site in path_starts with a sentinel at the end
  */
  std::vector<node_type> reference;
  std::vector<size_type> ref_starts, ref_ends;
  std::vector<node_type> alt_paths;
  std::vector<size_type> path_starts, site_starts;

  VariantPaths();

  size_type size() const { return this->reference.size(); }
  size_type paths() const { return this->path_starts.size() - 1; }
  size_type sites() const { return this->site_starts.size() - 1; }
  size_type alleles(size_type site) const { return this->site_starts[site + 1] - this->site_starts[site]; }

  void setReferenceSize(size_type size) { this->reference.reserve(size); }
  void appendReference(node_type node) { this->reference.push_back(node); }

  void addSite(size_type ref_start, size_type ref_end);
  void addAllele(const std::vector<node_type>& path);

  void appendReferenceUntil(Haplotype& haplotype, size_type site) const;
  void appendVariant(Haplotype& haplotype, size_type site, size_type allele) const;
};

//------------------------------------------------------------------------------

struct Phasing
{
  size_type first, second;
  bool      diploid, phased;

  const static size_type HAPLOID  = 0;
  const static size_type UNPHASED = 1;
  const static size_type PHASED   = 2;

  Phasing() {}
  explicit Phasing(size_type allele) : first(allele), second(0), diploid(false), phased(false) {}

  Phasing(size_type first_allele, size_type second_allele, bool is_phased = true) :
    first(first_allele), second(second_allele), diploid(true), phased(is_phased)
  {
  }

  size_type encode(size_type max_allele) const;
  void decode(size_type code, size_type max_allele);
  static size_type maxCode(size_type max_allele);
};

//------------------------------------------------------------------------------

/*
  Phasing information for a number of samples in a temporary file.
*/

struct PhasingInformation
{
  // Header
  size_type sample_count, sample_offset; 
  size_type sites;

  // File
  std::string                filename;
  sdsl::int_vector_buffer<8> data;
  const static std::string   TEMP_FILE_PREFIX;  // "phasing"

  // Iterator
  size_type            site, data_offset;
  std::vector<Phasing> phasings;

  explicit PhasingInformation(range_type sample_range);
  ~PhasingInformation();

  // Append the phasings for a new site.
  void append(const std::vector<Phasing>& new_site);

  // Iterate over the phasings.
  void begin() { this->site = 0; this->data_offset = 0; this->read(); }
  void operator++() { this->site++; this->read(); }
  const Phasing& operator[](size_type i) const { return this->phasings[i]; }
  size_type offset() const { return this->site; }

  // Statistics.
  size_type samples() const { return this->sample_count; }
  range_type range() const { return range_type(this->sample_offset, this->sample_offset + this->sample_count - 1); }
  size_type size() const { return this->sites; }

  // FIXME implement safe versions of these
  PhasingInformation(const PhasingInformation&) = delete;
  PhasingInformation& operator= (const PhasingInformation&) = delete;

private:
  void read();
};

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_VARIANTS_H
