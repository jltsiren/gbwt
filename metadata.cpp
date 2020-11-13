/*
  Copyright (c) 2019, 2020 Jouni Siren

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

namespace gbwt
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr std::uint32_t Metadata::TAG;
constexpr std::uint32_t Metadata::VERSION;
constexpr std::uint64_t Metadata::FLAG_MASK;
constexpr std::uint64_t Metadata::FLAG_PATH_NAMES;
constexpr std::uint64_t Metadata::FLAG_SAMPLE_NAMES;
constexpr std::uint64_t Metadata::FLAG_CONTIG_NAMES;
constexpr std::uint32_t Metadata::INITIAL_VERSION;
constexpr std::uint64_t Metadata::INITIAL_FLAG_MASK;

//------------------------------------------------------------------------------

Metadata::Metadata() :
  tag(TAG), version(VERSION),
  sample_count(0), haplotype_count(0), contig_count(0),
  flags(0)
{
}

Metadata::Metadata(std::vector<const Metadata*> sources, bool same_samples, bool same_contigs)  :
  tag(TAG), version(VERSION),
  sample_count(0), haplotype_count(0), contig_count(0),
  flags(0)
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

  written_bytes += sdsl::write_member(this->tag, out, child, "tag");
  written_bytes += sdsl::write_member(this->version, out, child, "version");
  written_bytes += sdsl::write_member(this->sample_count, out, child, "sample_count");
  written_bytes += sdsl::write_member(this->haplotype_count, out, child, "haplotype_count");
  written_bytes += sdsl::write_member(this->contig_count, out, child, "contig_count");
  written_bytes += sdsl::write_member(this->flags, out, child, "flags");

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
  sdsl::read_member(this->tag, in);
  sdsl::read_member(this->version, in);
  sdsl::read_member(this->sample_count, in);
  sdsl::read_member(this->haplotype_count, in);
  sdsl::read_member(this->contig_count, in);
  sdsl::read_member(this->flags, in);

  // Check the header.
  if(!(this->check()))
  {
    std::cerr << "Metadata::load(): Invalid metadata: " << *this << std::endl;
  }
  this->setVersion(); // Update to the current version.

  if(this->hasPathNames())
  {
    loadVector(this->path_names, in);
  }
  if(this->hasSampleNames())
  {
    this->sample_names.load(in);
  }
  if(this->hasContigNames())
  {
    this->contig_names.load(in);
  }
}

bool
Metadata::check() const
{
  if(this->tag != TAG) { return false; }
  switch(this->version)
  {
  case VERSION:
    return ((this->flags & FLAG_MASK) == this->flags);
  case INITIAL_VERSION:
    return ((this->flags & INITIAL_FLAG_MASK) == this->flags);
  default:
    return false;
  }
}

void
Metadata::swap(Metadata& another)
{
  if(this != &another)
  {
    std::swap(this->tag, another.tag);
    std::swap(this->version, another.version);
    std::swap(this->sample_count, another.sample_count);
    std::swap(this->haplotype_count, another.haplotype_count);
    std::swap(this->contig_count, another.contig_count);
    std::swap(this->flags, another.flags);

    this->path_names.swap(another.path_names);
    this->sample_names.swap(another.sample_names);
    this->contig_names.swap(another.contig_names);
  }
}

bool
Metadata::operator==(const Metadata& another) const
{
  return (this->tag == another.tag &&
          this->version == another.version &&
          this->sample_count == another.sample_count &&
          this->haplotype_count == another.haplotype_count &&
          this->contig_count == another.contig_count &&
          this->flags == another.flags &&
          this->path_names == another.path_names &&
          this->sample_names == another.sample_names &&
          this->contig_names == another.contig_names);
}

//------------------------------------------------------------------------------

void
Metadata::setSamples(size_type n)
{
  if(this->hasSampleNames() && Verbosity::level >= Verbosity::FULL)
  {
    std::cerr << "Metadata::setSamples(): Warning: Changing sample count without changing sample names" << std::endl;
  }
  this->sample_count = n;
}

void
Metadata::setHaplotypes(size_type n)
{
  this->haplotype_count = n;
}

void
Metadata::setContigs(size_type n)
{
  if(this->hasContigNames() && Verbosity::level >= Verbosity::FULL)
  {
    std::cerr << "Metadata::setContigs(): Warning: Changing contig count without changing contig names" << std::endl;
  }
  this->contig_count = n;
}

//------------------------------------------------------------------------------

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
  this->set(FLAG_PATH_NAMES);
  this->path_names.emplace_back(path);
}

void
Metadata::addPath(size_type sample, size_type contig, size_type phase, size_type count)
{
  this->set(FLAG_PATH_NAMES);
  PathName path =
  {
    static_cast<PathName::path_name_type>(sample),
    static_cast<PathName::path_name_type>(contig),
    static_cast<PathName::path_name_type>(phase),
    static_cast<PathName::path_name_type>(count)
  };
  this->path_names.emplace_back(path);
}

void
Metadata::clearPathNames()
{
  this->unset(FLAG_PATH_NAMES);
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
  this->set(FLAG_SAMPLE_NAMES);
  this->sample_names = Dictionary(names);
}

void
Metadata::addSamples(const std::vector<std::string>& names)
{
  if(names.empty()) { return; }

  Dictionary additional_names(names);
  this->sample_names.append(additional_names);
  this->setSamples(this->sample_names.size());
  this->set(FLAG_SAMPLE_NAMES);
}

void
Metadata::clearSampleNames()
{
  this->unset(FLAG_SAMPLE_NAMES);
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
  this->set(FLAG_CONTIG_NAMES);
  this->contig_names = Dictionary(names);
}

void
Metadata::addContigs(const std::vector<std::string>& names)
{
  if(names.empty()) { return; }

  Dictionary additional_names(names);
  this->contig_names.append(additional_names);
  this->setContigs(this->contig_names.size());
  this->set(FLAG_CONTIG_NAMES);
}

void
Metadata::clearContigNames()
{
  this->unset(FLAG_CONTIG_NAMES);
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
    std::set<size_type> phases;
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
  this->sample_count--;
  this->haplotype_count -= haplotypes_to_remove;

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
  this->contig_count--;

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
      double added_samples = this->sample_names.size() - this->sample_count;
      this->haplotype_count += (added_samples * source.haplotypes()) / source.samples();
    }
    this->sample_count = this->sample_names.size();
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
      this->set(FLAG_SAMPLE_NAMES);
    }
  }
  else
  {
    source_sample_offset = this->samples();
    this->sample_count += source.samples();
    this->haplotype_count += source.haplotypes();
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
    this->contig_count = this->contig_names.size();
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
      this->set(FLAG_CONTIG_NAMES);
    }
  }
  else
  {
    source_contig_offset = this->contigs();
    this->contig_count += source.contigs();
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
      std::set<std::pair<size_type, size_type>> found_haplotypes;
      for(const PathName& path : this->path_names)
      {
        found_haplotypes.emplace(path.sample, path.phase);
      }
      this->haplotype_count = found_haplotypes.size();
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
  if(metadata.get(Metadata::FLAG_PATH_NAMES))
  {
    stream << metadata.paths() << " paths with names, ";
  }

  stream << metadata.samples() << " samples";
  if(metadata.get(Metadata::FLAG_SAMPLE_NAMES)) { stream << " with names"; }
  stream << ", ";

  stream << metadata.haplotypes() << " haplotypes, ";

  stream << metadata.contigs() << " contigs";
  if(metadata.get(Metadata::FLAG_CONTIG_NAMES)) { stream << " with names"; }

  return stream;
}

//------------------------------------------------------------------------------

} // namespace gbwt
