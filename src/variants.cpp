/*
  Copyright (c) 2018, 2019, 2021 Jouni Siren

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

// Numerical class constants.

constexpr std::uint32_t VariantPaths::TAG;
constexpr std::uint32_t VariantPaths::VERSION;
constexpr std::uint64_t VariantPaths::FLAG_MASK;
constexpr std::uint64_t VariantPaths::FLAG_SAMPLE_NAMES;
constexpr std::uint64_t VariantPaths::FLAG_CONTIG_NAME;

constexpr size_type Phasing::HAPLOID;
constexpr size_type Phasing::UNPHASED;
constexpr size_type Phasing::PHASED;

//------------------------------------------------------------------------------

// Other class variables.

const std::string PhasingInformation::TEMP_FILE_PREFIX = "phasing";

//------------------------------------------------------------------------------

VariantPaths::VariantPaths() :
  tag(TAG), version(VERSION), flags(0),
  ref_index(16, wang_hash_64)
{
  this->path_starts.push_back(0);
  this->site_starts.push_back(0);
}

VariantPaths::VariantPaths(size_type reference_size) :
  tag(TAG), version(VERSION), flags(0),
  ref_index(16, wang_hash_64)
{
  this->reference.reserve(reference_size);
  this->path_starts.push_back(0);
  this->site_starts.push_back(0);
}

VariantPaths::VariantPaths(const VariantPaths& source)
{
  this->copy(source);
}

VariantPaths::VariantPaths(VariantPaths&& source)
{
  *this = std::move(source);
}

VariantPaths::~VariantPaths()
{
}

void
VariantPaths::swap(VariantPaths& another)
{
  if(this != &another)
  {
    std::swap(this->tag, another.tag);
    std::swap(this->version, another.version);
    std::swap(this->flags, another.flags);

    this->reference.swap(another.reference);
    this->ref_starts.swap(another.ref_starts);
    this->ref_ends.swap(another.ref_ends);
    this->alt_paths.swap(another.alt_paths);
    this->path_starts.swap(another.path_starts);
    this->site_starts.swap(another.site_starts);

    this->phasing_files.swap(another.phasing_files);
    this->file_offsets.swap(another.file_offsets);
    this->file_counts.swap(another.file_counts);

    this->ref_index.swap(another.ref_index);
  }
}

VariantPaths&
VariantPaths::operator=(const VariantPaths& source)
{
  if(this != &source) { this->copy(source); }
  return *this;
}

VariantPaths&
VariantPaths::operator=(VariantPaths&& source)
{
  if(this != &source)
  {
    this->tag = std::move(source.tag);
    this->version = std::move(source.version);
    this->flags = std::move(source.flags);

    this->reference = std::move(source.reference);
    this->ref_starts = std::move(source.ref_starts);
    this->ref_ends = std::move(source.ref_ends);
    this->alt_paths = std::move(source.alt_paths);
    this->path_starts = std::move(source.path_starts);
    this->site_starts = std::move(source.site_starts);

    this->phasing_files = std::move(source.phasing_files);
    this->file_offsets = std::move(source.file_offsets);
    this->file_counts = std::move(source.file_counts);

    this->ref_index = std::move(source.ref_index);
  }
  return *this;
}

size_type
VariantPaths::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += sdsl::write_member(this->tag, out, child, "tag");
  written_bytes += sdsl::write_member(this->version, out, child, "version");
  written_bytes += sdsl::write_member(this->flags, out, child, "flags");

  written_bytes += serializeVector(this->reference, out, child, "reference");
  written_bytes += serializeVector(this->ref_starts, out, child, "ref_starts");
  written_bytes += serializeVector(this->ref_ends, out, child, "ref_ends");
  written_bytes += serializeVector(this->alt_paths, out, child, "alt_paths");
  written_bytes += serializeVector(this->path_starts, out, child, "path_starts");
  written_bytes += serializeVector(this->site_starts, out, child, "site_starts");

  written_bytes += serializeVector(this->phasing_files, out, child, "phasing_files");
  written_bytes += serializeVector(this->file_offsets, out, child, "file_offsets");
  written_bytes += serializeVector(this->file_counts, out, child, "file_counts");

  if(this->hasSampleNames())
  {
    written_bytes += serializeVector(this->sample_names, out, child, "sample_names");
  }
  if(this->hasContigName())
  {
    written_bytes += sdsl::write_member(this->contig_name, out, child, "contig_name");
  }

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
VariantPaths::load(std::istream& in)
{
  sdsl::read_member(this->tag, in);
  sdsl::read_member(this->version, in);
  sdsl::read_member(this->flags, in);

  // Check the header and update to the current version.
  this->check();
  this->setVersion();

  loadVector(this->reference, in);
  loadVector(this->ref_starts, in);
  loadVector(this->ref_ends, in);
  loadVector(this->alt_paths, in);
  loadVector(this->path_starts, in);
  loadVector(this->site_starts, in);

  loadVector(this->phasing_files, in);
  loadVector(this->file_offsets, in);
  loadVector(this->file_counts, in);

  if(this->hasSampleNames()) { loadVector(this->sample_names, in); }
  if(this->hasContigName()) { sdsl::read_member(this->contig_name, in); }
}

void
VariantPaths::check() const
{
  if(this->tag != TAG)
  {
    throw sdsl::simple_sds::InvalidData("VariantPaths: Invalid tag");
  }

  if(this->version != VERSION)
  {
    std::string msg = "VariantPaths: Expected v" + std::to_string(VERSION) + ", got v" + std::to_string(this->version);
    throw sdsl::simple_sds::InvalidData(msg);
  }

  std::uint64_t mask = 0;
  switch(this->version)
  {
  case VERSION:
    mask = FLAG_MASK; break;
  }
  if((this->flags & mask) != this->flags)
  {
    throw sdsl::simple_sds::InvalidData("VariantPaths: Invalid flags");
  }
}

void
VariantPaths::copy(const VariantPaths& source)
{
  this->tag = source.tag;
  this->version = source.version;
  this->flags = source.flags;

  this->reference = source.reference;
  this->ref_starts = source.ref_starts;
  this->ref_ends = source.ref_ends;
  this->alt_paths = source.alt_paths;
  this->path_starts = source.path_starts;
  this->site_starts = source.site_starts;

  this->phasing_files = source.phasing_files;
  this->file_offsets = source.file_offsets;
  this->file_counts = source.file_counts;

  this->ref_index = source.ref_index;
}

//------------------------------------------------------------------------------

size_type
VariantPaths::nodeWidth(bool both_orientations) const
{
  node_type max_node = 0;
  for(auto node : this->reference)
  {
    node_type temp = node;
    max_node = std::max(temp, max_node);
    if(both_orientations) { max_node = std::max(Node::reverse(temp), max_node); }
  }
  for(auto node : this->alt_paths)
  {
    node_type temp = node;
    max_node = std::max(temp, max_node);
    if(both_orientations) { max_node = std::max(Node::reverse(temp), max_node); }
  }
  return sdsl::bits::length(max_node);
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
VariantPaths::addAllele(const vector_type& path)
{
  this->alt_paths.insert(this->alt_paths.end(), path.begin(), path.end());
  this->path_starts.push_back(this->alt_paths.size());
  this->site_starts.back() = this->path_starts.size() - 1;
}

void
VariantPaths::addFile(const std::string& filename, size_type sample_offset, size_type sample_count)
{
  this->phasing_files.push_back(filename);
  this->file_offsets.push_back(sample_offset);
  this->file_counts.push_back(sample_count);
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
VariantPaths::appendVariant(Haplotype& haplotype, size_type site, size_type allele,
  std::function<void(const Haplotype&)> output, std::function<bool(size_type, size_type)> report_overlap) const
{
  if(allele == 0) { return; }
  if(site >= this->sites())
  {
    std::cerr << "VariantPaths::appendVariant(): Invalid site: " << site << std::endl;
    return;
  }
  if(allele > this->alleles(site))
  {
    std::cerr << "VariantPaths::appendVariant(): Invalid allele at site " << site << ": " << allele << std::endl;
    return;
  }

  /*
    The alternate allele may overlap with the previous alternate allele. Try to resolve
    the overlap by removing reference nodes from the end of the current path and by
    ignoring reference nodes at the start of the alt path. If this fails, create a phase
    break or ignore the alternate allele.
  */
  size_type start = this->path_starts[this->site_starts[site] + allele - 1];
  size_type stop = this->path_starts[this->site_starts[site] + allele];
  if(haplotype.offset > this->ref_starts[site])
  {
    size_type overlap = haplotype.offset - this->ref_starts[site];
    size_type ends_with_reference = 0, starts_with_reference = 0;
    size_type potential_ref_nodes = std::min(haplotype.offset, static_cast<size_type>(haplotype.path.size()));
    while(overlap > 0 && ends_with_reference < potential_ref_nodes)
    {
      node_type haplotype_node = haplotype.path[haplotype.path.size() - ends_with_reference - 1];
      node_type reference_node = this->reference[haplotype.offset - ends_with_reference - 1];
      if(haplotype_node == reference_node) { ends_with_reference++; overlap--; }
      else { break; }
    }
    while(overlap > 0 && start + starts_with_reference < stop)
    {
      node_type haplotype_node = this->alt_paths[start + starts_with_reference];
      node_type reference_node = this->reference[this->ref_starts[site] + starts_with_reference];
      if(haplotype_node == reference_node) { starts_with_reference++; overlap--; }
      else { break; }
    }
    if(overlap == 0) // We can resolve the overlap.
    {
      while(ends_with_reference > 0) { haplotype.path.pop_back(); haplotype.offset--; ends_with_reference--; }
      start += starts_with_reference;
    }
    else
    {
      if(report_overlap(site, allele)) { return; }
      else // Phase break.
      {
        output(haplotype);
        haplotype.deactivate();
        haplotype.activate(this->ref_starts[site], haplotype.diploid);
      }
    }
  }
  else // No overlap, just append reference nodes.
  {
    this->appendReferenceUntil(haplotype, site);
  }

  // Append the alternate allele.
  haplotype.path.insert(haplotype.path.end(), this->alt_paths.begin() + start, this->alt_paths.begin() + stop);
  haplotype.offset = this->ref_ends[site];
}

vector_type
VariantPaths::getAllele(size_type site, size_type allele) const
{
  if(site >= this->sites()) { return vector_type(); }

  if(allele == 0)
  {
    size_type start = this->refStart(site);
    size_type stop = this->refEnd(site);
    return vector_type(this->reference.begin() + start, this->reference.begin() + stop);
  }
  else
  {
    size_type start = this->path_starts[this->site_starts[site] + allele - 1];
    size_type stop = this->path_starts[this->site_starts[site] + allele];
    return vector_type(this->alt_paths.begin() + start, this->alt_paths.begin() + stop);
  }
}

//------------------------------------------------------------------------------

void
VariantPaths::setSampleNames(const std::vector<std::string>& new_names)
{
  this->set(FLAG_SAMPLE_NAMES);
  this->sample_names = new_names;
}

void
VariantPaths::setContigName(const std::string& new_name)
{
  this->set(FLAG_CONTIG_NAME);
  this->contig_name = new_name;
}

//------------------------------------------------------------------------------

vector_type
getVector(const vector_type& source, bool get_ids)
{
  if(!get_ids) { return source; }
  vector_type result = source;
  for(auto& node : result) { node = Node::id(node); }
  return result;
}

void
checkOverlaps(const VariantPaths& variants, std::ostream& out, bool print_ids)
{
  size_type prev_start = variants.refStart(0), prev_end = variants.refEnd(0);
  size_type overlapping_sites = 0, unresolved_sites = 0;
  std::map<size_type, size_type> overlaps_by_length;
  for(size_type site = 1; site < variants.sites(); site++)
  {
    size_type start = variants.refStart(site), end = variants.refEnd(site);
    bool site_resolved = true;
    if(start < prev_end) // We have an overlap.
    {
      overlapping_sites++;
      size_type overlap = prev_end - start;
      overlaps_by_length[overlap]++;
      for(size_type prev_allele = 1; prev_allele <= variants.alleles(site - 1); prev_allele++)
      {
        vector_type prev_path = variants.getAllele(site - 1, prev_allele);
        size_type ends_with_reference = 0;
        while(ends_with_reference < prev_path.size() && ends_with_reference < prev_end)
        {
          node_type alt_node = prev_path[prev_path.size() - ends_with_reference - 1];
          node_type ref_node = variants.refAt(prev_end - ends_with_reference - 1);
          if(alt_node != ref_node) { break; }
          ends_with_reference++;
        }
        for(size_type allele = 1; allele <= variants.alleles(site); allele++)
        {
          vector_type path = variants.getAllele(site, allele);
          size_type starts_with_reference = 0;
          while(starts_with_reference < path.size() && start + starts_with_reference < end)
          {
            node_type alt_node = path[starts_with_reference];
            node_type ref_node = variants.refAt(start + starts_with_reference);
            if(alt_node != ref_node) { break; }
            starts_with_reference++;
          }
          if(ends_with_reference + starts_with_reference < overlap)
          {
            vector_type prev_ref = variants.getAllele(site - 1, 0);
            vector_type curr_ref = variants.getAllele(site, 0);
            out << "Site " << site << ", alleles " << range_type(prev_allele, allele) << ":" << std::endl;
            out << "  overlap " << overlap << " nodes: "
                << range_type(prev_start, prev_end) << " followed by " << range_type(start, end) << std::endl;
            out << "  prev ref: " << getVector(prev_ref, print_ids) << std::endl;
            out << "  prev alt: " << getVector(prev_path, print_ids) << std::endl;
            out << "  curr ref: " << getVector(curr_ref, print_ids) << std::endl;
            out << "  curr alt: " << getVector(path, print_ids) << std::endl;
            site_resolved = false;
          }
        }
      }
    }
    prev_start = start; prev_end = end;
    if(!site_resolved) { unresolved_sites++; }
  }

  out << "Sites: " << variants.sites() << " total, "
      << overlapping_sites << " overlapping, "
      << unresolved_sites << " unresolved" << std::endl;
  out << "Overlap lengths:" << std::endl;
  for(auto iter = overlaps_by_length.begin(); iter != overlaps_by_length.end(); ++iter)
  {
    out << "  " << iter->first << " nodes: " << iter->second << " overlaps" << std::endl;
  }
}

//------------------------------------------------------------------------------

Phasing::Phasing(const std::string& genotype, bool was_diploid, bool phase_homozygous) :
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

  try
  {
    this->first = std::stoul(genotype.substr(0, separator_pos));
  }
  catch(...)
  {
    this->first = 0; // Fall back to reference.
  }
  if(this->diploid)
  {
    try
    {
      this->second = std::stoul(genotype.substr(separator_pos + 1));
    }
    catch(...)
    {
      this->second = 0; // Fall back to reference.
    }
    if(phase_homozygous && this->first == this->second)
    {
      this->phased = true;
    }
  }
}

void
Phasing::forcePhased(std::function<bool()> rng)
{
  if(this->phased) { return; }
  this->phased = true;
  if(rng()) { std::swap(this->first, this->second); }
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

std::ostream&
operator<<(std::ostream& out, Phasing phasing)
{
  out << phasing.first;
  if(phasing.diploid)
  {
    out << (phasing.phased ? "|" : "/") << phasing.second;
  }
  return out;
}

//------------------------------------------------------------------------------

PhasingInformation::PhasingInformation() :
  sample_count(0), sample_offset(0), site_count(0),
  filename(), data(), temp_file(false),
  site(0), data_offset(0), phasings(0)
{
}

PhasingInformation::PhasingInformation(size_type first_sample, size_type num_samples) :
  sample_count(num_samples), sample_offset(first_sample), site_count(0),
  filename(TempFile::getName(TEMP_FILE_PREFIX)), data(filename, std::ios::out), temp_file(true),
  site(0), data_offset(0), phasings(sample_count)
{
}

PhasingInformation::PhasingInformation(const std::string& base_name, size_type first_sample, size_type num_samples) :
  sample_count(num_samples), sample_offset(first_sample), site_count(0),
  filename(base_name + '_' + std::to_string(first_sample) + '_' + std::to_string(num_samples)),
  data(filename, std::ios::out), temp_file(false),
  site(0), data_offset(0), phasings(sample_count)
{
}

PhasingInformation::PhasingInformation(const VariantPaths& variants, size_type file)
{
  if(file >= variants.files())
  {
    std::cerr << "PhasingInformation::PhasingInformation(): Invalid file number: " << file << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if(Verbosity::level >= Verbosity::BASIC)
  {
    std::cerr << "PhasingInformation::PhasingInformation(): Opening file " << variants.name(file) << std::endl;
  }

  this->sample_count = variants.count(file);
  this->sample_offset = variants.offset(file);
  this->site_count = variants.sites();

  this->filename = variants.name(file);
  this->data = sdsl::int_vector_buffer<8>(this->filename, std::ios::in);
  this->temp_file = false;

  this->site = 0;
  this->data_offset = 0;
  this->phasings = std::vector<Phasing>(this->sample_count);
}

PhasingInformation::PhasingInformation(PhasingInformation&& source)
{
  *this = std::move(source);
}

PhasingInformation::~PhasingInformation()
{
  this->close();
  if(this->temp_file) { TempFile::remove(this->filename); }
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
    this->temp_file = source.temp_file;
    source.filename.clear();
    source.temp_file = false;
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
  if(this->isOpen()) { this->data.close(false); this->data = sdsl::int_vector_buffer<8>(); }
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
                   std::function<bool(size_type)> process_sample, std::function<void(const Haplotype&)> output,
                   std::function<bool(size_type, size_type)> report_overlap)
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
        variants.appendVariant(first, phasings.current(), phasing.first, output, report_overlap);
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
          variants.appendVariant(second, phasings.current(), phasing.second, output, report_overlap);
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

void
generateHaplotypes(const VariantPaths& variants, const std::set<std::string>& phasing_files,
                   std::function<bool(size_type)> process_sample, std::function<void(const Haplotype&)> output,
                   std::function<bool(size_type, size_type)> report_overlap)
{
  for(size_type file = 0; file < variants.files(); file++)
  {
    if(phasing_files.empty() || phasing_files.find(variants.name(file)) != phasing_files.end())
    {
      PhasingInformation phasings(variants, file);
      generateHaplotypes(variants, phasings, process_sample, output, report_overlap);
    }
  }
}

//------------------------------------------------------------------------------

} // namespace gbwt
