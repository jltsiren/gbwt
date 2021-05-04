/*
  Copyright (c) 2019, 2020, 2021 Jouni Siren

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

#include <gbwt/files.h>
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

class Metadata
{
public:
  typedef gbwt::size_type size_type; // Needed for SDSL serialization.

  // Header.
  MetadataHeader        header;

  // Path / sample / contig names.
  std::vector<PathName> path_names;
  Dictionary            sample_names;
  Dictionary            contig_names;

  Metadata();
  Metadata(std::vector<const Metadata*> sources, bool same_samples, bool same_contigs);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  void simple_sds_serialize(std::ostream& out) const;
  void simple_sds_load(std::istream& in);
  size_t simple_sds_size() const;

  void swap(Metadata& another);

  bool operator==(const Metadata& another) const;
  bool operator!=(const Metadata& another) const { return !(this->operator==(another)); }

  // Header operations.
  size_type samples() const { return this->header.sample_count; }
  size_type haplotypes() const { return this->header.haplotype_count; }
  size_type contigs() const { return this->header.contig_count; }
  void setSamples(size_type n);
  void setHaplotypes(size_type n);
  void setContigs(size_type n);

  // Path operations.
  bool hasPathNames() const { return this->header.get(MetadataHeader::FLAG_PATH_NAMES); }
  size_type paths() const { return this->path_names.size(); }
  const PathName& path(size_type i) const { return this->path_names[i]; }
  std::vector<size_type> findPaths(size_type sample_id, size_type contig_id) const;
  std::vector<size_type> pathsForSample(size_type sample_id) const;
  std::vector<size_type> pathsForContig(size_type contig_id) const;
  void addPath(const PathName& path);
  void addPath(size_type sample, size_type contig, size_type phase, size_type count);
  void clearPathNames();

  // Sample operations.
  bool hasSampleNames() const { return this->header.get(MetadataHeader::FLAG_SAMPLE_NAMES); }
  std::string sample(size_type i) const { return this->sample_names[i]; }
  size_type sample(const std::string& name) const { return this->sample_names.find(name); }
  void setSamples(const std::vector<std::string>& names);
  void addSamples(const std::vector<std::string>& names);
  void clearSampleNames();

  // Contig operations.
  bool hasContigNames() const { return this->header.get(MetadataHeader::FLAG_CONTIG_NAMES); }
  std::string contig(size_type i) const { return this->contig_names[i]; }
  size_type contig(const std::string& name) const { return this->contig_names.find(name); }
  void setContigs(const std::vector<std::string>& names);
  void addContigs(const std::vector<std::string>& names);
  void clearContigNames();

  // Remove metadata corresponding to the sample/contig and return the set of removed
  // path identifiers.
  // NOTE: This may invalidate existing sample/contig ids and path ids. When removing
  // multiple samples/contigs, remove the corresponding sequences from the GBWT
  // before proceeding to the next sample/contig.
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
