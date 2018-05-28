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

  size_type paths() const { return this->path_starts.size() - 1; }
  size_type sites() const { return this->site_starts.size() - 1; }
  size_type alleles(size_type site) const { return this->site_starts[site + 1] - this->site_starts[site]; }

  void addSite(size_type ref_start, size_type ref_end);
  void addAllele(const std::vector<node_type>& path);

  void appendReferenceUntil(Haplotype& haplotype, size_type site) const;
  void appendVariant(Haplotype& haplotype, size_type site, size_type allele) const;
};

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_VARIANTS_H
