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

#include <string>

namespace gbwt
{

//------------------------------------------------------------------------------

VariantPaths::VariantPaths(size_type reference_size) :
  ref_index(16, wang_hash_64)
{
  this->reference.reserve(reference_size);
  this->path_starts.push_back(0);
  this->site_starts.push_back(0);
}

void
VariantPaths::indexReference()
{
  if(!(this->ref_index.empty())) { this->ref_index.clear(); }
  for(size_type i = 0; i < this->size(); i++)
  {
    if(this->ref_index.find(this->reference[i]) == this->ref_index.end()) { this->ref_index[this->reference[i]] = i; }
  }
}

size_type
VariantPaths::firstOccurrence(node_type node)
{
  auto iter = this->ref_index.find(node);
  if(iter == this->ref_index.end()) { return this->invalid_position(); }
  return iter->second;
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

Phasing::Phasing(const std::string& genotype, bool was_diploid) :
  first(0), second(0), diploid(was_diploid), phased(true)
{
  if(genotype.empty()) { return; }

  size_type separator_pos = genotype.find('|');
  this->phased = (separator_pos != std::string::npos);
  if(this->phased) { this->diploid = true; }
  else
  {
    separator_pos = genotype.find('/');
    this->diploid = (separator_pos != std::string::npos);
  }
  if(!(this->diploid)) { this->phased = true; }

  this->first = std::stoul(genotype.substr(0, separator_pos));
  if(this->diploid) { this->second = std::stoul(genotype.substr(separator_pos + 1)); }
}

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

PhasingInformation::PhasingInformation(size_type first_sample, size_type num_samples) :
  sample_count(num_samples), sample_offset(first_sample), site_count(0),
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
  this->close();
  TempFile::remove(this->filename);
}

PhasingInformation&
PhasingInformation::operator=(PhasingInformation&& source)
{
  if(this != &source)
  {
    this->close();
    bool open_file = source.isOpen();
    source.close();

    this->sample_count = source.sample_count;
    this->sample_offset = source.sample_offset;
    this->site_count = source.site_count;

    // Transfer the ownership of the file.
    this->filename = source.filename;
    source.filename.clear();
    if(open_file) { this->open(); }

    this->site = source.site;
    this->data_offset = source.data_offset;
    this->phasings = source.phasings;
  }
  return *this;
}

void
PhasingInformation::open()
{
  if(!(this->isOpen())) { this->data = sdsl::int_vector_buffer<8>(this->filename, std::ios::in); }
}

void PhasingInformation::close()
{
  if(this->isOpen()) { this->data = sdsl::int_vector_buffer<8>(); }
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
  if(!(this->isOpen()))
  {
    std::cerr << "PhasingInformation::append(): Cannot append to a closed file" << std::endl;
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
  Run encoder(Phasing::maxCode(max_allele) + 1);
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
  if(!(this->isOpen()))
  {
    std::cerr << "PhasingInformation::read(): Cannot read from a closed file" << std::endl;
    return;
  }

  size_type max_allele = ByteCode::read(this->data, this->data_offset);
  Run decoder(Phasing::maxCode(max_allele) + 1);
  run_type run = decoder.read(this->data, this->data_offset);
  for(size_type i = 0; i < this->size(); i++)
  {
    if(run.second == 0) { run = decoder.read(this->data, this->data_offset); }
    this->phasings[i].decode(run.first, max_allele); run.second--;
  }
}

//------------------------------------------------------------------------------

void
finishHaplotype(Haplotype& haplotype, const VariantPaths& variants, size_type site,
                std::function<void(const Haplotype&)> output)
{
  if(!(haplotype.active)) { return; }
  if(site < variants.sites()) { variants.appendReferenceUntil(haplotype, site); }
  else { variants.appendReferenceUntilEnd(haplotype); }
  output(haplotype);
  haplotype.deactivate();
}

void
generateHaplotypes(const VariantPaths& variants, PhasingInformation& phasings,
                   std::function<bool(size_type)> process_sample, std::function<void(const Haplotype&)> output)
{
  phasings.open();

  // Initialize diploid haplotypes.
  std::vector<Haplotype> haplotypes;
  haplotypes.reserve(2 * phasings.size());
  for(size_type sample = 0; sample < phasings.size(); sample++)
  {
    haplotypes.emplace_back(phasings.offset() + sample, 0);
    haplotypes.emplace_back(phasings.offset() + sample, 1);
  }

  // Determine which samples to process.
  sdsl::bit_vector active_samples(phasings.size(), 0);
  for(size_type sample = 0; sample < phasings.size(); sample++) { active_samples[sample] = process_sample(sample); }

  phasings.begin();
  while(phasings.current() < phasings.sites())
  {
    for(size_type sample = 0; sample < phasings.size(); sample++)
    {
      if(!active_samples[sample]) { continue; }
      const Phasing& phasing = phasings[sample];
      Haplotype& first = haplotypes[2 * sample];
      Haplotype& second = haplotypes[2 * sample + 1];

      /*
        First phase:
        - If the current site is unphased or the ploidy changes, finish the existing haplotype.
        - If the current haplotype is inactive, activate it.
        - If the current haplotype has an alternate allele, append it.
        - If the current site is unphased, finish the current haplotype.
      */
      if(!(phasing.phased) || first.diploid != phasing.diploid)
      {
        finishHaplotype(first, variants, phasings.current(), output);
      }
      if(!(first.active)) { first.activate(variants.refPrev(phasings.current()), phasing.diploid); }
      if(phasing.first > 0)
      {
        variants.appendVariant(first, phasings.current(), phasing.first);
      }
      if(!(phasing.phased))
      {
        finishHaplotype(first, variants, phasings.current() + 1, output);
      }

      /*
        Second phase:
        - If the current site is unphased or haploid, finish the existing haplotype.
        - If the current site is haploid, skip the rest.
        - If the current haplotype is inactive, activate it.
        - If the current haplotype has an alternate allele, append it.
        - If the current site is unphased, finish the current haplotype.
      */
      if(!(phasing.phased && phasing.diploid))
      {
        finishHaplotype(second, variants, phasings.current(), output);
      }
      if(phasing.diploid)
      {
        if(!(second.active)) { second.activate(variants.refPrev(phasings.current()), phasing.diploid); }
        if(phasing.second > 0)
        {
          variants.appendVariant(second, phasings.current(), phasing.second);
        }
        if(!(phasing.phased))
        {
          finishHaplotype(second, variants, phasings.current() + 1, output);
        }
      }
    }
    ++phasings;
  }

  // Finish the active haplotypes.
  for(size_type sample = 0; sample < phasings.size(); sample++)
  {
    if(!process_sample(sample)) { continue; }
    finishHaplotype(haplotypes[2 * sample], variants, phasings.sites(), output);
    finishHaplotype(haplotypes[2 * sample + 1], variants, phasings.sites(), output);
  }

  phasings.close();
}

//------------------------------------------------------------------------------

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
  VariantPaths variants(reference.size());
  for(node_type node : reference) { variants.appendToReference(node); }

  // Add sites and alternate alleles.
  std::vector<range_type> sites = { range_type(2, 4), range_type(6, 7), range_type(8, 9) };
  std::vector<std::vector<std::vector<node_type>>> alleles = { { { 11, 12 } }, { { 13 }, { 14, 15 } }, { { 16 }, {} } };
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

  // Test the reference index.
  variants.indexReference();
  for(size_type i = 0; i < reference.size(); i++)
  {
    if(variants.firstOccurrence(reference[i]) != i)
    {
      std::cerr << "testVariants(): VariantPaths: First occurrence of " << reference[i] << " at " << variants.firstOccurrence(reference[i]) << ", expected " << i << std::endl;
      failures++;
    }
  }
  if(variants.firstOccurrence(0) != variants.invalid_position())
  {
    std::cerr << "testVariants(): VariantPaths: Reference position found for an invalid node" << std::endl;
    failures++;
  }

  // Create 3 samples: diploid-haploid-diploid, phased-unphased-phased, haploid-phased-phased.
  size_type first_sample = 10, num_samples = 3;
  std::vector<std::vector<Phasing>> phasing_information =
  {
    { Phasing(1, 0, true), Phasing(1, 1, true),  Phasing(0) },
    { Phasing(0),          Phasing(1, 2, false), Phasing(2, 0, true) },
    { Phasing(2, 1, true), Phasing(1, 0, true),  Phasing(0, 1, true) }
  };
  PhasingInformation phasings(first_sample, num_samples);
  for(const std::vector<Phasing>& site : phasing_information)
  {
    phasings.append(site);
  }

  // Close and reopen the phasings.
  if(!(phasings.isOpen()))
  {
    std::cerr << "testVariants(): PhasingInformation: The file is not open after construction" << std::endl;
    failures++;
  }
  phasings.close();
  if(phasings.isOpen())
  {
    std::cerr << "testVariants(): PhasingInformation: The file was not properly closed" << std::endl;
    failures++;
  }
  phasings.open();
  if(!(phasings.isOpen()))
  {
    std::cerr << "testVariants(): PhasingInformation: The file was not properly reopened" << std::endl;
    failures++;
  }

  // Test phasing statistics.
  if(phasings.size() != num_samples)
  {
    std::cerr << "testVariants(): PhasingInformation: Sample count " << phasings.size() << ", expected " << num_samples << std::endl;
    failures++;
  }
  if(phasings.offset() != first_sample)
  {
    std::cerr << "testVariants(): PhasingInformation: Sample offset " << phasings.offset() << ", expected " << first_sample << std::endl;
    failures++;
  }
  if(phasings.limit() != first_sample + num_samples)
  {
    std::cerr << "testVariants(): PhasingInformation: Sample limit " << phasings.limit() << ", expected " << (first_sample + num_samples) << std::endl;
    failures++;
  }
  if(phasings.sites() != phasing_information.size())
  {
    std::cerr << "testVariants(): PhasingInformation: Site count " << phasings.sites() << ", expected " << phasing_information.size() << std::endl;
    failures++;
  }

  // Test genotype string parsing.
  std::vector<std::vector<std::string>> genotypes =
  {
    { "1|0", "1|1", "0" },
    { "0",   "1/2", "2|0" },
    { "2|1", "1|0", "0|1" }
  };
  for(size_type site = 0; site < genotypes.size(); site++)
  {
    for(size_type sample = 0; sample < genotypes[site].size(); sample++)
    {
      if(Phasing(genotypes[site][sample]) != phasing_information[site][sample])
      {
        std::cerr << "testVariants(): Phasing: Parsing failure: " << genotypes[site][sample] << std::endl;
        failures++;
      }
    }
  }

  /*
  1 2  3  4 5 6  7    8  9 10
      11 12     13      16
                14 15    -
  */

  // Generate haplotypes.
  std::vector<std::vector<node_type>> true_haplotypes =
  {
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
  };
  std::vector<std::vector<node_type>> haplotypes;
  generateHaplotypes(variants, phasings,
    [](size_type) -> bool { return true; },
    [&haplotypes](const Haplotype& haplotype) { haplotypes.push_back(haplotype.path); });
  if(haplotypes.size() != true_haplotypes.size())
  {
    std::cerr << "testVariants(): generateHaplotypes(): Expected " << true_haplotypes.size() << " haplotypes, got " << haplotypes.size() << std::endl;
    failures++;
  }
  else
  {
    for(size_type haplotype = 0; haplotype < haplotypes.size(); haplotype++)
    {
      if(haplotypes[haplotype] != true_haplotypes[haplotype])
      {
        std::cerr << "testVariants(): generateHaplotypes(): Haplotype " << haplotype << ": expected " << true_haplotypes[haplotype] << ", got " << haplotypes[haplotype] << std::endl;
        failures++;
      }
    }
  }

  return failures;
}

//------------------------------------------------------------------------------

} // namespace gbwt
