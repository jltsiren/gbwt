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

#include <functional>

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
  this->site_starts.push_back(this->path_starts.size() - 1);
}

void
VariantPaths::addAllele(const std::vector<node_type>& path)
{
  this->alt_paths.insert(this->alt_paths.end(), path.begin(), path.end());
  this->path_starts.push_back(this->alt_paths.size());
  this->site_starts.back() = this->path_starts.size() - 1;
}

void
VariantPaths::appendReferenceUntil(Haplotype& haplotype, size_type site) const
{
  if(site >= this->sites())
  {
    std::cerr << "VariantPaths::appendReferenceUntil(): Invalid site: " << site << std::endl;
    return;
  }

  size_type until = this->ref_starts[site];
  if(haplotype.offset < until)
  {
    haplotype.path.insert(haplotype.path.end(), this->reference.begin() + haplotype.offset, this->reference.begin() + until);
    haplotype.offset = until;
  }
}

void
VariantPaths::appendReferenceUntilEnd(Haplotype& haplotype) const
{
  if(haplotype.offset < this->size())
  {
    haplotype.path.insert(haplotype.path.end(), this->reference.begin() + haplotype.offset, this->reference.end());
    haplotype.offset = this->size();
  }
}

void
VariantPaths::appendVariant(Haplotype& haplotype, size_type site, size_type allele) const
{
  if(site >= this->sites())
  {
    std::cerr << "VariantPaths::appendVariant(): Invalid site: " << site << std::endl;
    return;
  }
  if(allele == 0 || allele > this->alleles(site))
  {
    std::cerr << "VariantPaths::appendVariant(): Invalid allele at site " << site << ": " << allele << std::endl;
    return;
  }

  this->appendReferenceUntil(haplotype, site);
  size_type site_start = this->site_starts[site];
  size_type start = this->path_starts[site_start + allele - 1], stop = this->path_starts[site_start + allele];
  haplotype.path.insert(haplotype.path.end(), this->alt_paths.begin() + start, this->alt_paths.begin() + stop);
  haplotype.offset = this->ref_ends[site];
}

//------------------------------------------------------------------------------

size_type
Phasing::encode(size_type max_allele) const
{
  size_type code = HAPLOID;
  if(this->diploid)
  {
    if(this->phased) { code = PHASED; }
    else { code = UNPHASED; }
  }

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
  this->phased = (code != UNPHASED);
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
  sample_count(Range::length(sample_range)), sample_offset(sample_range.first), site_count(0),
  filename(TempFile::getName(TEMP_FILE_PREFIX)), data(filename, std::ios::out),
  site(0), data_offset(0), phasings(sample_count)
{
}

PhasingInformation::PhasingInformation(PhasingInformation&& source)
{
  *this = std::move(source);
}

PhasingInformation::~PhasingInformation()
{
  this->data.close();
  TempFile::remove(this->filename);
}

PhasingInformation&
PhasingInformation::operator=(PhasingInformation&& source)
{
  if(this != &source)
  {
    this->sample_count = source.sample_count;
    this->sample_offset = source.sample_offset;
    this->site_count = source.site_count;

    // Transfer the ownership of the file.
    this->filename = source.filename;
    source.filename.clear();
    source.data.close();
    this->data = sdsl::int_vector_buffer<8>(this->filename, std::ios::in);

    this->site = source.site;
    this->data_offset = source.data_offset;
    this->phasings = source.phasings;
  }
  return *this;
}

void
PhasingInformation::append(const std::vector<Phasing>& new_site)
{
  if(new_site.empty()) { return; }
  if(new_site.size() != this->size())
  {
    std::cerr << "PhasingInformation::append(): Expected " << this->size() << " samples, got " << new_site.size() << std::endl;
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

  this->site_count++;
}

void PhasingInformation::read()
{
  if(this->site >= this->sites()) { return; }

  size_type max_allele = ByteCode::read(this->data, this->data_offset);
  Run decoder(Phasing::maxCode(max_allele));
  run_type run = decoder.read(this->data, this->data_offset);
  for(size_type i = 0; i < this->size(); i++)
  {
    if(run.second == 0) { run = decoder.read(this->data, this->data_offset); }
    this->phasings[i].decode(run.first, max_allele); run.second--;
  }
}

//------------------------------------------------------------------------------

void
extractHaplotypes(const VariantPaths& variants, PhasingInformation& phasings, std::vector<Haplotype>& haplotypes, std::function<bool(const Phasing&)> condition)
{
  phasings.begin();
  while(phasings.offset() < phasings.sites())
  {
    for(size_type sample = 0; sample < phasings.size(); sample++)
    {
      if(!condition(phasings[sample])) { continue; }
      if(phasings[sample].first > 0 && variants.refStart(phasings.offset()) >= haplotypes[2 * sample].offset)
      {
        variants.appendVariant(haplotypes[2 * sample], phasings.offset(), phasings[sample].first);
      }
      if(phasings[sample].diploid && phasings[sample].second > 0 && variants.refStart(phasings.offset()) >= haplotypes[2 * sample + 1].offset)
      {
        variants.appendVariant(haplotypes[2 * sample + 1], phasings.offset(), phasings[sample].second);
      }
    }
    ++phasings;
  }
  for(size_type sample = 0; sample < phasings.size(); sample++)
  {
    if(!condition(phasings[sample])) { continue; }
    variants.appendReferenceUntilEnd(haplotypes[2 * sample]);
    if(phasings[sample].diploid)
    {
      variants.appendReferenceUntilEnd(haplotypes[2 * sample + 1]);
    }
  }
}

std::ostream&
operator<< (std::ostream& out, std::vector<node_type>& data)
{
  out << "{ ";
  for(node_type node : data) { out << node << " "; }
  out << "}";
  return out;
}

size_type
testVariants()
{
  size_type failures = 0;

  // Add a reference.
  std::vector<node_type> reference = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
  VariantPaths variants;
  variants.setReferenceSize(reference.size());
  for(node_type node : reference) { variants.appendToReference(node); }

  // Add sites and alternate alleles.
  std::vector<range_type> sites = { range_type(2, 4), range_type(6, 7) };
  std::vector<std::vector<std::vector<node_type>>> alleles = { { { 11, 12 } }, { { 13 }, { 14, 15 } } };
  size_type allele_count = 0;
  for(size_type site = 0; site < sites.size(); site++)
  {
    variants.addSite(sites[site].first, sites[site].second);
    for(const std::vector<node_type>& allele : alleles[site])
    {
      variants.addAllele(allele);
      allele_count++;
    }
  }

  // Test variant statistics.
  if(variants.size() != reference.size())
  {
    std::cerr << "testVariants(): VariantPaths: Reference size " << variants.size() << ", expected " << reference.size() << std::endl;
    failures++;
  }
  if(variants.paths() != allele_count)
  {
    std::cerr << "testVariants(): VariantPaths: Allele count " << variants.paths() << ", expected " << allele_count << std::endl;
    failures++;
  }
  if(variants.sites() != sites.size())
  {
    std::cerr << "testVariants(): VariantPaths: Site count " << variants.sites() << ", expected " << sites.size() << std::endl;
    failures++;
  }
  for(size_type site = 0; site < sites.size(); site++)
  {
    if(variants.alleles(site) != alleles[site].size())
    {
      std::cerr << "testVariants(): VariantPaths: Site " << site << ": Allele count " << variants.alleles(site) << ", expected " << alleles[site].size() << std::endl;
      failures++;
    }
    if(variants.refStart(site) != sites[site].first)
    {
      std::cerr << "testVariants(): VariantPaths: Site " << site << ": Start " << variants.refStart(site) << ", expected " << sites[site].first << std::endl;
      failures++;
    }
    if(variants.refEnd(site) != sites[site].second)
    {
      std::cerr << "testVariants(): VariantPaths: Site " << site << ": End " << variants.refEnd(site) << ", expected " << sites[site].second << std::endl;
      failures++;
    }
  }

  // Create 3 phased samples: haploid, unphased, phased,
  range_type sample_range(10, 12);
  std::vector<std::vector<Phasing>> phasing_information =
  {
    { Phasing(1), Phasing(0, 1, false), Phasing(1, 0, true) },
    { Phasing(0), Phasing(1, 2, false), Phasing(2, 0, true) }
  };
  PhasingInformation phasings(sample_range);
  for(const std::vector<Phasing>& site : phasing_information)
  {
    phasings.append(site);
  }

  // Test phasing statistics.
  if(phasings.size() != Range::length(sample_range))
  {
    std::cerr << "testVariants(): PhasingInformation: Sample count " << phasings.size() << ", expected " << Range::length(sample_range) << std::endl;
    failures++;
  }
  if(phasings.range() != sample_range)
  {
    std::cerr << "testVariants(): PhasingInformation: Sample range " << phasings.range() << ", expected " << sample_range << std::endl;
    failures++;
  }
  if(phasings.sites() != phasing_information.size())
  {
    std::cerr << "testVariants(): PhasingInformation: Site count " << phasings.sites() << ", expected " << phasing_information.size() << std::endl;
    failures++;
  }

  // Extract full haplotypes.
  std::vector<std::vector<node_type>> full_haplotypes =
  {
    { 1, 2, 11, 12, 5, 6, 7, 8, 9, 10 },
    {},
    { 1, 2, 3, 4, 5, 6, 13, 8, 9, 10 },
    { 1, 2, 11, 12, 5, 6, 14, 15, 8, 9, 10 },
    { 1, 2, 11, 12, 5, 6, 14, 15, 8, 9, 10 },
    { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 }
  };
  std::vector<Haplotype> extracted_haplotypes(2 * phasings.size());
  extractHaplotypes(variants, phasings, extracted_haplotypes, [](const Phasing&) -> bool { return true; } );
  for(size_type haplotype = 0; haplotype < extracted_haplotypes.size(); haplotype++)
  {
    if(extracted_haplotypes[haplotype].path != full_haplotypes[haplotype])
    {
      std::cerr << "testVariants(): Full haplotype " << haplotype << ": " << extracted_haplotypes[haplotype].path << ", expected " << full_haplotypes[haplotype] << std::endl;
      failures++;
    }
  }

  // Extract haplotypes starting after the first site for phased samples.
  std::vector<std::vector<node_type>> partial_haplotypes =
  {
    { 5, 6, 7, 8, 9, 10 },
    {},
    {},
    {},
    { 5, 6, 14, 15, 8, 9, 10 },
    { 5, 6, 7, 8, 9, 10 }
  };
  extracted_haplotypes = std::vector<Haplotype>(2 * phasings.size(), Haplotype(variants.refEnd(0)));
  extractHaplotypes(variants, phasings, extracted_haplotypes, [](const Phasing& phasing) -> bool { return phasing.phased; } );
  for(size_type haplotype = 0; haplotype < extracted_haplotypes.size(); haplotype++)
  {
    if(extracted_haplotypes[haplotype].path != partial_haplotypes[haplotype])
    {
      std::cerr << "testVariants(): Partial haplotype " << haplotype << ": " << extracted_haplotypes[haplotype].path << ", expected " << partial_haplotypes[haplotype] << std::endl;
      failures++;
    }
  }

  return failures;
}

//------------------------------------------------------------------------------

} // namespace gbwt
