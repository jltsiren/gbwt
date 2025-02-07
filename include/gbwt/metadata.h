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

#ifndef GBWT_METADATA_H
#define GBWT_METADATA_H

#include <gbwt/files.h>
#include <gbwt/support.h>

#include <unordered_map>

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

  PathName() : sample(0), contig(0), phase(0), count(0) {}
  PathName(path_name_type sample, path_name_type contig, path_name_type phase, path_name_type count) :
    sample(sample), contig(contig), phase(phase), count(count)
  {
  }

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

// This is a standalone version of PathName. The fields are named in a way that is
// compatible with GFA walk lines.
struct FullPathName
{
  std::string sample_name;
  std::string contig_name;
  size_t      haplotype;
  size_t      offset;
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
  PathName path(const FullPathName& name) const;
  FullPathName fullPath(size_type i) const;
  size_type findFragment(const PathName& name) const;
  size_type findFragment(const FullPathName& name) const;
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

private:
  // Throws `sdsl::simple_sds::InvalidData` if the checks fail.
  void sanityChecks() const;
};

std::ostream& operator<<(std::ostream& stream, const Metadata& metadata);

//------------------------------------------------------------------------------

// Helper data structure for fragmented haplotype sequences.
// The basic unit is a chain of paths. Each chain corresponds to a distinct
// (sample, contig, phase) combination, with the paths ordered by the count field.
struct FragmentMap
{
  struct Fragment
  {
    size_type path; // Path identifier for the fragment.
    size_type chain; // Running identifier for the chain.
    size_type prev, next; // Path identifiers for the previous and next fragment in the chain.
  };

  // Maps a path identifier to fragment information.
  std::unordered_map<size_type, Fragment> fragments;

  // Number chains in the metadata.
  // Chain identifiers are in the range [0, chains).
  size_type chains;

  // Builds a FragmentMap from the path names in the metadata.
  explicit FragmentMap(const Metadata& metadata, bool verbose);

  FragmentMap() = default;
  FragmentMap(const FragmentMap&) = default;
  FragmentMap(FragmentMap&&) = default;
  FragmentMap& operator=(const FragmentMap&) = default;
  FragmentMap& operator=(FragmentMap&&) = default;

  // Returns the number of chains.
  size_type size() const { return this->chains; }

  // Returns true if the map is empty.
  bool empty() const { return (this->chains == 0); }

  // Returns the path identifier for the next fragment in the given orientation.
  // Returns invalid_sequence() if there is no next fragment.
  size_type next(size_type path_id, bool is_reverse = false) const;

  // Returns the sequence identifier for the next fragment.
  // Returns invalid_sequence() if there is no next fragment.
  size_type oriented_next(size_type sequence_id) const;

  // Returns the path identifier for the previous fragment in the given orientation.
  // Returns invalid_sequence() if there is no previous fragment.
  size_type prev(size_type path_id, bool is_reverse = false) const;

  // Returns the sequence identifier for the previous fragment.
  // Returns invalid_sequence() if there is no previous fragment.
  size_type oriented_prev(size_type sequence_id) const;

  // Returns the chain identifier for the fragment.
  // Returns invalid_sequence() if the fragment is not in the map.
  size_type chain(size_type path_id) const;
};

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_METADATA_H
