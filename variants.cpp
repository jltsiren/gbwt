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
#include <gbwt/internal.h>

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

size_type
Phasing::encode(size_type max_allele) const
{
  size_type code = HAPLOID;
  if(this->phased) { code = PHASED; }
  else if(this->diploid) { code = UNPHASED; }

  size_type radix = 3;
  code += radix * this->first;
  radix *= max_allele + 1;
  code += radix * this->second;

  return code;
}

void
Phasing::decode(size_type code, size_type max_allele)
{
  size_type radix = 3 * (max_allele + 1);
  this->second = code / radix;
  code -= radix * this->second;

  radix = 3;
  this->first = code / radix;
  code -= radix * this->first;

  this->diploid = (code != HAPLOID);
  this->phased = (code == PHASED);
}

size_type
Phasing::maxCode(size_type max_allele)
{
  Phasing phasing(max_allele, max_allele, true);
  return phasing.encode(max_allele);
}

//------------------------------------------------------------------------------

const std::string PhasingInformation::TEMP_FILE_PREFIX = "phasing";

PhasingInformation::PhasingInformation(range_type sample_range) :
  sample_count(Range::length(sample_range)), sample_offset(sample_range.first), sites(0),
  filename(TempFile::getName(TEMP_FILE_PREFIX)), data(filename, std::ios::out),
  site(0), data_offset(0), phasings(sample_count)
{
}

PhasingInformation::~PhasingInformation()
{
  this->data.close();
  TempFile::remove(this->filename);
}

void
PhasingInformation::append(const std::vector<Phasing>& new_site)
{
  if(new_site.empty()) { return; }
  if(new_site.size() != this->samples())
  {
    std::cerr << "PhasingInformation::append(): Expected " << this->samples() << " samples, got " << new_site.size() << std::endl;
    return;
  }

  // Determine and write the largest allele identifier.
  size_type max_allele = 0;
  for(const Phasing& phasing : new_site)
  {
    max_allele = std::max(max_allele, phasing.first);
    max_allele = std::max(max_allele, phasing.second);
  }
  ByteCode::write(this->data, max_allele);

  // Run-length encode the phasings.
  Run encoder(Phasing::maxCode(max_allele));
  size_type prev = new_site.front().encode(max_allele), run_length = 1;
  for(size_type i = 1; i < new_site.size(); i++)
  {
    size_type curr = new_site[i].encode(max_allele);
    if(curr == prev) { run_length++; }
    else
    {
      encoder.write(this->data, prev, run_length);
      prev = curr; run_length = 1;
    }
  }
  encoder.write(this->data, prev, run_length);  

  this->sites++;
}

void PhasingInformation::read()
{
  size_type max_allele = ByteCode::read(this->data, this->data_offset);

  Run decoder(Phasing::maxCode(max_allele));
  run_type run = decoder.read(this->data, this->data_offset);
  for(size_type i = 0; i < this->samples(); i++)
  {
    if(run.second == 0) { run = decoder.read(this->data, this->data_offset); }
    this->phasings[i].decode(run.first, max_allele); run.second--;
  }
}

//------------------------------------------------------------------------------

} // namespace gbwt
