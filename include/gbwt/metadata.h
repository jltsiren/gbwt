/*
  Copyright (c) 2019 Jouni Siren

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

#ifndef GBWT_METADATA_H
#define GBWT_METADATA_H

#include <gbwt/support.h>

namespace gbwt
{

/*
  metadata.h: Metadata structure.
*/

//------------------------------------------------------------------------------

struct PathName
{
#ifdef GBWT_SAVE_MEMORY
  typedef short_type path_name_type;
#else
  typedef size_type path_name_type;
#endif

  path_name_type sample;
  path_name_type contig;
  path_name_type phase;
  path_name_type count;

  bool operator==(const PathName& another) const
  {
    return (this->sample == another.sample && this->contig == another.contig && this->phase == another.phase && this->count == another.count);
  }

  bool operator!=(const PathName& another) const { return !(this->operator==(another)); }

  bool operator<(const PathName& another) const
  {
    if(this->sample != another.sample) { return (this->sample < another.sample); }
    if(this->contig != another.contig) { return (this->contig < another.contig); }
    if(this->phase != another.phase) { return (this->phase < another.phase); }
    return (this->count < another.count);
  }
};

//------------------------------------------------------------------------------

/*
  Version 1:
  - Sample names, contig names, path names.
  - Compatible with version 0.

  Version 0:
  - Preliminary version with sample/haplotype/contig counts.

  Future versions:
  - Haplotype coverage over ranges of node ids (run-length encode with sd_vector).
  - Assign contig ids to ranges of node ids (as above).
*/

class Metadata
{
public:
  typedef gbwt::size_type size_type; // Needed for SDSL serialization.

  // Header.
  std::uint32_t tag;
  std::uint32_t version;
  std::uint64_t sample_count;
  std::uint64_t haplotype_count;
  std::uint64_t contig_count;
  std::uint64_t flags;

  // Path / sample / contig names.
  std::vector<PathName> path_names;
  Dictionary            sample_names;
  Dictionary            contig_names;

  constexpr static std::uint32_t TAG = 0x6B375E7A;
  constexpr static std::uint32_t VERSION = Version::METADATA_VERSION;

  constexpr static std::uint64_t FLAG_MASK         = 0x0007;
  constexpr static std::uint64_t FLAG_PATH_NAMES   = 0x0001;
  constexpr static std::uint64_t FLAG_SAMPLE_NAMES = 0x0002;
  constexpr static std::uint64_t FLAG_CONTIG_NAMES = 0x0004;

  // Flag masks for old compatible versions.
  constexpr static std::uint32_t INITIAL_VERSION   = 0;
  constexpr static std::uint64_t INITIAL_FLAG_MASK = 0x0000;

  Metadata();
  Metadata(std::vector<const Metadata*> sources, bool same_samples, bool same_contigs);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);
  bool check() const;

  void setVersion() { this->version = VERSION; }

  void set(std::uint64_t flag) { this->flags |= flag; }
  void unset(std::uint64_t flag) { this->flags &= ~flag; }
  bool get(std::uint64_t flag) const { return (this->flags & flag); }

  void swap(Metadata& another);

  bool operator==(const Metadata& another) const;
  bool operator!=(const Metadata& another) const { return !(this->operator==(another)); }

  // Header operations.
  size_type samples() const { return this->sample_count; }
  size_type haplotypes() const { return this->haplotype_count; }
  size_type contigs() const { return this->contig_count; }
  void setSamples(size_type n);
  void setHaplotypes(size_type n);
  void setContigs(size_type n);

  // Path operations.
  bool hasPathNames() const { return this->get(FLAG_PATH_NAMES); }
  size_type paths() const { return this->path_names.size(); }
  const PathName& path(size_type i) const { return this->path_names[i]; }
  std::vector<size_type> findPaths(size_type sample_id, size_type contig_id) const;
  std::vector<size_type> pathsForSample(size_type sample_id) const;
  std::vector<size_type> pathsForContig(size_type contig_id) const;
  void addPath(const PathName& path);
  void clearPathNames();

  // Sample operations.
  bool hasSampleNames() const { return this->get(FLAG_SAMPLE_NAMES); }
  std::string sample(size_type i) const { return this->sample_names[i]; }
  size_type sample(const std::string& name) const { return this->sample_names.find(name); }
  void setSamples(const std::vector<std::string>& names);
  void addSamples(const std::vector<std::string>& names);
  void clearSampleNames();

  // Contig operations.
  bool hasContigNames() const { return this->get(FLAG_CONTIG_NAMES); }
  std::string contig(size_type i) const { return this->contig_names[i]; }
  size_type contig(const std::string& name) const { return this->contig_names.find(name); }
  void setContigs(const std::vector<std::string>& names);
  void addContigs(const std::vector<std::string>& names);
  void clearContigNames();

  // Remove metadata corresponding to the sample/contig and return the set of removed
  // path identifiers.
  std::vector<size_type> removeSample(size_type sample_id);
  std::vector<size_type> removeContig(size_type contig_id);

  // Merge the metadata from the source into this object.
  // If the objects to be merged both contain sample / contig names, this overides the
  // same_samples / same_contigs flags.
  void merge(const Metadata& source, bool same_samples, bool same_contigs);

  void clear();
};

std::ostream& operator<<(std::ostream& stream, const Metadata& metadata);

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_METADATA_H
