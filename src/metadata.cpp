/*
  Copyright (c) 2019, 2020, 2021, 2024, 2025 Jouni Siren

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

#include <gbwt/internal.h>
#include <gbwt/metadata.h>

#include <set>
#include <unordered_set>

namespace gbwt
{

//------------------------------------------------------------------------------

Metadata::Metadata()
{
}

Metadata::Metadata(std::vector<const Metadata*> sources, bool same_samples, bool same_contigs)
{
  if(sources.empty()) { return; }

  *this = *(sources.front());
  for(size_type i = 1; i < sources.size(); i++)
  {
    this->merge(*(sources[i]), same_samples, same_contigs);
  }
}

size_type
Metadata::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += this->header.serialize(out, child, "header");

  if(this->hasPathNames())
  {
    written_bytes += serializeVector(this->path_names, out, child, "path_names");
  }
  if(this->hasSampleNames())
  {
    written_bytes += this->sample_names.serialize(out, child, "sample_names");
  }
  if(this->hasContigNames())
  {
    written_bytes += this->contig_names.serialize(out, child, "contig_names");
  }

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
Metadata::load(std::istream& in)
{
  // Load and check the header.
  this->header.load(in);
  this->header.check();
  bool old_version = (this->header.version < MetadataHeader::VERSION);
  this->header.setVersion(); // Update to the current version.

  if(this->hasPathNames())
  {
    loadVector(this->path_names, in);
  }
  if(this->hasSampleNames())
  {
    if(old_version) { this->sample_names.load_v1(in); }
    else { this->sample_names.load(in); }
  }
  if(this->hasContigNames())
  {
    if(old_version) { this->contig_names.load_v1(in); }
    else { this->contig_names.load(in); }
  }

  this->sanityChecks();
}

void
Metadata::simple_sds_serialize(std::ostream& out) const
{
  sdsl::simple_sds::serialize_value(this->header, out);
  sdsl::simple_sds::serialize_vector(this->path_names, out);
  this->sample_names.simple_sds_serialize(out);
  this->contig_names.simple_sds_serialize(out);
}

void
Metadata::simple_sds_load(std::istream& in)
{
  // Header.
  this->header = sdsl::simple_sds::load_value<MetadataHeader>(in);
  this->header.check_simple_sds();
  this->header.setVersion(); // Update to the current version.

  // Path / sample / contig names.
  this->path_names = sdsl::simple_sds::load_vector<PathName>(in);
  this->sample_names.simple_sds_load(in);
  this->contig_names.simple_sds_load(in);

  this->sanityChecks();
}

size_t
Metadata::simple_sds_size() const
{
  return sdsl::simple_sds::value_size(this->header) +
    sdsl::simple_sds::vector_size(this->path_names) +
    this->sample_names.simple_sds_size() +
    this->contig_names.simple_sds_size();
}

void
Metadata::swap(Metadata& another)
{
  if(this != &another)
  {
    this->header.swap(another.header);
    this->path_names.swap(another.path_names);
    this->sample_names.swap(another.sample_names);
    this->contig_names.swap(another.contig_names);
  }
}

bool
Metadata::operator==(const Metadata& another) const
{
  return (this->header == another.header &&
          this->path_names == another.path_names &&
          this->sample_names == another.sample_names &&
          this->contig_names == another.contig_names);
}

void
Metadata::sanityChecks() const
{
  if(!(this->hasPathNames()) && this->path_names.size() > 0)
  {
    throw sdsl::simple_sds::InvalidData("Metadata: Invalid path name flag in the header");
  }

  if(this->hasSampleNames())
  {
    if(this->header.sample_count != this->sample_names.size())
    {
      throw sdsl::simple_sds::InvalidData("Metadata: Sample / sample name count mismatch");
    }
  }
  else if(this->sample_names.size() > 0)
  {
    throw sdsl::simple_sds::InvalidData("Metadata: Invalid sample name flag in the header");
  }

  if(this->hasContigNames())
  {
    if(this->header.contig_count != this->contig_names.size())
    {
      throw sdsl::simple_sds::InvalidData("Metadata: Contig / contig name count mismatch");
    }
  }
  else if(this->contig_names.size() > 0)
  {
    throw sdsl::simple_sds::InvalidData("Metadata: Invalid contig name flag in the header");
  }
}

//------------------------------------------------------------------------------

void
Metadata::setSamples(size_type n)
{
  if(this->hasSampleNames() && n != this->sample_names.size() && Verbosity::level >= Verbosity::FULL)
  {
    std::cerr << "Metadata::setSamples(): Warning: Changing sample count without changing sample names" << std::endl;
  }
  this->header.sample_count = n;
}

void
Metadata::setHaplotypes(size_type n)
{
  this->header.haplotype_count = n;
}

void
Metadata::setContigs(size_type n)
{
  if(this->hasContigNames() && n != this->contig_names.size() && Verbosity::level >= Verbosity::FULL)
  {
    std::cerr << "Metadata::setContigs(): Warning: Changing contig count without changing contig names" << std::endl;
  }
  this->header.contig_count = n;
}

//------------------------------------------------------------------------------

PathName
Metadata::path(const FullPathName& name) const
{
  size_type sample = (this->hasSampleNames() ? this->sample(name.sample_name) : this->samples());
  size_type contig = (this->hasContigNames() ? this->contig(name.contig_name) : this->contigs());
  return PathName(sample, contig, name.haplotype, name.offset);
}

FullPathName
Metadata::fullPath(size_type i) const
{
  const PathName& path = this->path(i);
  std::string sample_name = (this->hasSampleNames() ? this->sample(path.sample) : std::to_string(path.sample));
  std::string contig_name = (this->hasContigNames() ? this->contig(path.contig) : std::to_string(path.contig));
  size_t haplotype = path.phase;
  size_t offset = path.count;
  return FullPathName { sample_name, contig_name, haplotype, offset };
}

size_type
Metadata::findFragment(const PathName& name) const
{
  size_type result = this->paths();
  if(!(this->hasPathNames())) { return result; }

  for(size_type i = 0; i < this->paths(); i++)
  {
    const PathName& path = this->path(i);
    if(path.sample == name.sample && path.contig == name.contig && path.phase == name.phase && path.count <= name.count)
    {
      if(result >= this->paths() || path.count > this->path(result).count) { result = i; }
    }
  }

  return result;
}

size_type
Metadata::findFragment(const FullPathName& name) const
{
  return this->findFragment(this->path(name));
}

std::vector<size_type>
Metadata::findPaths(size_type sample_id, size_type contig_id) const
{
  std::vector<size_type> result;
  for(size_type i = 0; i < this->paths(); i++)
  {
    if(this->path(i).sample == sample_id && this->path(i).contig == contig_id)
    {
      result.push_back(i);
    }
  }
  return result;
}

std::vector<size_type>
Metadata::pathsForSample(size_type sample_id) const
{
  std::vector<size_type> result;
  for(size_type i = 0; i < this->paths(); i++)
  {
    if(this->path(i).sample == sample_id)
    {
      result.push_back(i);
    }
  }

  return result;
}

std::vector<size_type>
Metadata::pathsForContig(size_type contig_id) const
{
  std::vector<size_type> result;
  for(size_type i = 0; i < this->paths(); i++)
  {
    if(this->path(i).contig == contig_id)
    {
      result.push_back(i);
    }
  }
  return result;
}

void
Metadata::addPath(const PathName& path)
{
  this->header.set(MetadataHeader::FLAG_PATH_NAMES);
  this->path_names.push_back(path);
}

void
Metadata::addPath(size_type sample, size_type contig, size_type phase, size_type count)
{
  this->addPath(PathName(sample, contig, phase, count));
}

void
Metadata::clearPathNames()
{
  this->header.unset(MetadataHeader::FLAG_PATH_NAMES);
  this->path_names = std::vector<PathName>();
}

//------------------------------------------------------------------------------

void
Metadata::setSamples(const std::vector<std::string>& names)
{
  if(names.empty())
  {
    this->clearSampleNames();
    return;
  }

  this->setSamples(names.size());
  this->header.set(MetadataHeader::FLAG_SAMPLE_NAMES);
  this->sample_names = Dictionary(names);
}

void
Metadata::addSamples(const std::vector<std::string>& names)
{
  if(names.empty()) { return; }

  Dictionary additional_names(names);
  this->sample_names.append(additional_names);
  this->setSamples(this->sample_names.size());
  this->header.set(MetadataHeader::FLAG_SAMPLE_NAMES);
}

void
Metadata::clearSampleNames()
{
  this->header.unset(MetadataHeader::FLAG_SAMPLE_NAMES);
  this->sample_names.clear();
}

//------------------------------------------------------------------------------

void
Metadata::setContigs(const std::vector<std::string>& names)
{
  if(names.empty())
  {
    this->clearContigNames();
    return;
  }

  this->setContigs(names.size());
  this->header.set(MetadataHeader::FLAG_CONTIG_NAMES);
  this->contig_names = Dictionary(names);
}

void
Metadata::addContigs(const std::vector<std::string>& names)
{
  if(names.empty()) { return; }

  Dictionary additional_names(names);
  this->contig_names.append(additional_names);
  this->setContigs(this->contig_names.size());
  this->header.set(MetadataHeader::FLAG_CONTIG_NAMES);
}

void
Metadata::clearContigNames()
{
  this->header.unset(MetadataHeader::FLAG_CONTIG_NAMES);
  this->contig_names.clear();
}

//------------------------------------------------------------------------------

// Remove the path when callee() returns true.
// Returns the set of removed path identifiers.
// Also gives an opportunity to update the remaining paths.
std::vector<size_type>
removePaths(Metadata& metadata, std::function<bool(PathName&)> callee)
{
  std::vector<size_type> result;
  size_type tail = 0;
  for(size_type i = 0; i < metadata.path_names.size(); i++)
  {
    if(!callee(metadata.path_names[i]))
    {
      metadata.path_names[tail] = metadata.path_names[i];
      tail++;
    }
    else
    {
      result.push_back(i);
    }
  }

  if(tail > 0) { metadata.path_names.resize(tail); }
  else { metadata.clearPathNames(); }
  return result;
}

std::vector<size_type>
Metadata::removeSample(size_type sample_id)
{
  std::vector<size_type> result;
  if(sample_id >= this->samples()) { return result; }

  // Update paths, determine the number of removed haplotypes.
  size_type haplotypes_to_remove = 0;
  if(this->hasPathNames())
  {
    std::unordered_set<size_type> phases;
    result = gbwt::removePaths(*this, [sample_id, &phases, &result](PathName& path) -> bool {
      if(path.sample == sample_id)
      {
        phases.insert(path.phase);
        return true;
      }
      else
      {
        if(path.sample > sample_id) { path.sample--; }
        return false;
      }
    });
    haplotypes_to_remove = phases.size();
  }
  else
  {
    if(Verbosity::level >= Verbosity::FULL)
    {
      std::cerr << "Metadata::removeSample(): Warning: Estimating new haplotype count" << std::endl;
    }
    haplotypes_to_remove = static_cast<double>(this->haplotypes()) / this->samples();
  }

  // Update samples and haplotypes.
  if(this->hasSampleNames()) { this->sample_names.remove(sample_id); }
  this->header.sample_count--;
  this->header.haplotype_count -= haplotypes_to_remove;

  return result;
}

std::vector<size_type>
Metadata::removeContig(size_type contig_id)
{
  std::vector<size_type> result;
  if(contig_id >= this->contigs()) { return result; }

  // Update paths.
  if(this->hasPathNames())
  {
    result = gbwt::removePaths(*this, [contig_id, &result](PathName& path) -> bool {
      if(path.contig == contig_id) { return true; }
      else
      {
        if(path.contig > contig_id) { path.contig--; }
        return false;
      }
    });
  }

  // Update contigs.
  if(this->hasContigNames()) { this->contig_names.remove(contig_id); }
  this->header.contig_count--;

  return result;
}

//------------------------------------------------------------------------------

void
Metadata::merge(const Metadata& source, bool same_samples, bool same_contigs)
{
  size_type source_sample_offset = 0, source_contig_offset = 0;
  bool merge_sample_names = (this->hasSampleNames() & source.hasSampleNames());
  bool merge_contig_names = (this->hasContigNames() & source.hasContigNames());
  bool merge_path_names = (this->hasPathNames() & source.hasPathNames());

  // Merge samples and haplotypes.
  if(merge_sample_names)
  {
    this->sample_names = Dictionary(this->sample_names, source.sample_names);
    if(!merge_path_names)
    {
      if(Verbosity::level >= Verbosity::FULL)
      {
        std::cerr << "Metadata::merge(): Warning: Estimating new haplotype count" << std::endl;
      }
      double added_samples = this->sample_names.size() - this->header.sample_count;
      this->header.haplotype_count += (added_samples * source.haplotypes()) / source.samples();
    }
    this->header.sample_count = this->sample_names.size();
  }
  else if(same_samples)
  {
    if(this->samples() != source.samples() || this->haplotypes() != source.haplotypes())
    {
      if(Verbosity::level >= Verbosity::FULL)
      {
        std::cerr << "Metadata::merge(): Warning: Sample/haplotype counts do not match" << std::endl;
      }
    }
    if(!(this->hasSampleNames()) && source.hasSampleNames())
    {
      if(Verbosity::level >= Verbosity::FULL)
      {
        std::cerr << "Metadata::merge(): Warning: Taking sample names from the source" << std::endl;
      }
      this->sample_names = source.sample_names;
      this->header.set(MetadataHeader::FLAG_SAMPLE_NAMES);
    }
  }
  else
  {
    source_sample_offset = this->samples();
    this->header.sample_count += source.samples();
    this->header.haplotype_count += source.haplotypes();
    if(this->hasSampleNames())
    {
      if(Verbosity::level >= Verbosity::FULL)
      {
        std::cerr << "Metadata::merge(): Warning: Clearing sample names; the source has no sample names" << std::endl;
      }
      this->clearSampleNames();
    }
  }

  // Merge contigs.
  if(merge_contig_names)
  {
    this->contig_names = Dictionary(this->contig_names, source.contig_names);
    this->header.contig_count = this->contig_names.size();
  }
  else if(same_contigs)
  {
    if(this->contigs() != source.contigs() && Verbosity::level >= Verbosity::FULL)
    {
      std::cerr << "Metadata::merge(): Warning: Contig counts do not match" << std::endl;
    }
    if(!(this->hasContigNames()) && source.hasContigNames())
    {
      if(Verbosity::level >= Verbosity::FULL)
      {
        std::cerr << "Metadata::merge(): Warning: Taking contig names from the source" << std::endl;
      }
      this->contig_names = source.contig_names;
      this->header.set(MetadataHeader::FLAG_CONTIG_NAMES);
    }
  }
  else
  {
    source_contig_offset = this->contigs();
    this->header.contig_count += source.contigs();
    if(this->hasContigNames())
    {
      if(Verbosity::level >= Verbosity::FULL)
      {
        std::cerr << "Metadata::merge(): Warning: Clearing contig names; the source has no contig names" << std::endl;
      }
      this->clearContigNames();
    }
  }

  // Merge paths.
  if(merge_path_names)
  {
    size_type source_path_offset = this->paths();
    this->path_names.insert(this->path_names.end(), source.path_names.begin(), source.path_names.end());
    for(size_type i = source_path_offset; i < this->path_names.size(); i++)
    {
      if(merge_sample_names)
      {
        this->path_names[i].sample = this->sample(source.sample(this->path_names[i].sample));
      }
      else { this->path_names[i].sample += source_sample_offset; }
      if(merge_contig_names)
      {
        this->path_names[i].contig = this->contig(source.contig(this->path_names[i].contig));
      }
      else { this->path_names[i].contig += source_contig_offset; }
    }
    // Determine the new haplotype count.
    if(merge_sample_names)
    {
      std::set<std::pair<PathName::path_name_type, PathName::path_name_type>> found_haplotypes;
      for(const PathName& path : this->path_names)
      {
        found_haplotypes.emplace(path.sample, path.phase);
      }
      this->header.haplotype_count = found_haplotypes.size();
    }
  }
  else if(this->hasPathNames())
  {
    if(Verbosity::level >= Verbosity::FULL)
    {
      std::cerr << "Metadata::merge(): Warning: Clearing path names; the source has no path names" << std::endl;
    }
    this->clearPathNames();
  }
}

void
Metadata::clear()
{
  *this = Metadata();
}

//------------------------------------------------------------------------------

std::ostream& operator<<(std::ostream& stream, const Metadata& metadata)
{
  if(metadata.hasPathNames())
  {
    stream << metadata.paths() << " paths with names, ";
  }

  stream << metadata.samples() << " samples";
  if(metadata.hasSampleNames()) { stream << " with names"; }
  stream << ", ";

  stream << metadata.haplotypes() << " haplotypes, ";

  stream << metadata.contigs() << " contigs";
  if(metadata.hasContigNames()) { stream << " with names"; }

  return stream;
}

//------------------------------------------------------------------------------

FragmentMap::FragmentMap(const Metadata& metadata, bool verbose) :
  chains(0)
{
  std::vector<std::pair<PathName, size_type>> sorted_paths;
  for(size_type i = 0; i < metadata.paths(); i++)
  {
    sorted_paths.emplace_back(metadata.path(i), i);
  }
  sequentialSort(sorted_paths.begin(), sorted_paths.end());

  auto same_sequence = [&sorted_paths](size_type i, size_type j) -> bool
  {
    const PathName& a = sorted_paths[i].first;
    const PathName& b = sorted_paths[j].first;
    return (a.sample == b.sample && a.contig == b.contig && a.phase == b.phase);
  };

  for(size_type i = 0; i < sorted_paths.size(); i++)
  {
    if(i == 0 || !same_sequence(i, i - 1))
    {
      this->chains++;
    }
    Fragment curr { sorted_paths[i].second, this->chains - 1, invalid_sequence(), invalid_sequence() };
    if(i > 0 && same_sequence(i, i - 1))
    {
      curr.prev = sorted_paths[i - 1].second;
    }
    if(i + 1 < sorted_paths.size() && same_sequence(i, i + 1))
    {
      curr.next = sorted_paths[i + 1].second;
    }
    this->fragments[sorted_paths[i].second] = curr;
  }

  if(verbose)
  {
    size_type head = 0, tail = 1;
    while(head < sorted_paths.size())
    {
      while(tail < sorted_paths.size() && same_sequence(head, tail)) { tail++; }
      if(tail - head > 1)
      {
        const FullPathName name = metadata.fullPath(sorted_paths[head].second);
        std::cerr << name.sample_name << "#" << name.haplotype << "#" << name.contig_name << ": " << (tail - head) << " fragments" << std::endl;
      }
      head = tail; tail = head + 1;
    }
  }
}

size_type
FragmentMap::next(size_type path_id, bool is_reverse) const
{
  auto iter = this->fragments.find(path_id);
  if(iter == this->fragments.end()) { return invalid_sequence(); }
  return (is_reverse ? iter->second.prev : iter->second.next);
}

size_type
FragmentMap::oriented_next(size_type sequence_id) const
{
  size_type path_id = Path::id(sequence_id);
  bool is_reverse = Path::is_reverse(sequence_id);
  size_type result = this->next(path_id, is_reverse);
  return (result == invalid_sequence() ? invalid_sequence() : Path::encode(result, is_reverse));
}

size_type
FragmentMap::prev(size_type path_id, bool is_reverse) const
{
  auto iter = this->fragments.find(path_id);
  if(iter == this->fragments.end()) { return invalid_sequence(); }
  return (is_reverse ? iter->second.next : iter->second.prev);
}

size_type
FragmentMap::oriented_prev(size_type sequence_id) const
{
  size_type path_id = Path::id(sequence_id);
  bool is_reverse = Path::is_reverse(sequence_id);
  size_type result = this->prev(path_id, is_reverse);
  return (result == invalid_sequence() ? invalid_sequence() : Path::encode(result, is_reverse));
}

size_type
FragmentMap::chain(size_type path_id) const
{
  auto iter = this->fragments.find(path_id);
  if(iter == this->fragments.end()) { return invalid_sequence(); }
  return iter->second.chain;
}

//------------------------------------------------------------------------------

} // namespace gbwt
