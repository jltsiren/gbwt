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

#include <gbwt/variants.h>

namespace gbwt
{

//------------------------------------------------------------------------------

VariantPaths::VariantPaths()
{
  this->path_starts.push_back(0);
  this->site_starts.push_back(0);
}

void
VariantPaths::addSite(size_type ref_start, size_type ref_end)
{
  this->ref_starts.push_back(ref_start);
  this->ref_ends.push_back(ref_end);
  this->site_starts.push_back(this->path_starts.size());
}

void
VariantPaths::addAllele(const std::vector<node_type>& path)
{
  this->alt_paths.insert(this->alt_paths.end(), path.begin(), path.end());
  this->path_starts.push_back(this->alt_paths.size());
  this->site_starts.back() = this->path_starts.size();
}

void
VariantPaths::appendReferenceUntil(Haplotype& haplotype, size_type site) const
{
  size_type until = this->ref_starts[site];
  if(haplotype.offset < until)
  {
    haplotype.path.insert(haplotype.path.end(), this->reference.begin() + haplotype.offset, this->reference.begin() + until);
    haplotype.offset = until;
  }
}

void
VariantPaths::appendVariant(Haplotype& haplotype, size_type site, size_type allele) const
{
  this->appendReferenceUntil(haplotype, site);

  size_type site_start = this->site_starts[site];
  size_type start = this->path_starts[site_start + allele], stop = this->path_starts[site_start + allele + 1];
  haplotype.path.insert(haplotype.path.end(), this->alt_paths.begin() + start, this->alt_paths.begin() + stop);
  haplotype.offset = this->ref_ends[site];
}

//------------------------------------------------------------------------------

} // namespace gbwt
