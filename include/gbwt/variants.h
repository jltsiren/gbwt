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


#ifndef GBWT_VARIANTS_H
#define GBWT_VARIANTS_H

#include <gbwt/utils.h>

#include <functional>
#include <map>
#include <set>

namespace gbwt
{

/*
  variants.h: Support structures for transforming phased VCF files into paths.
*/

//------------------------------------------------------------------------------

struct Haplotype
{
  vector_type path;
  size_type   offset; // In the reference.
  bool        active, diploid;

  size_type   sample, phase, count;

  Haplotype() : offset(0), active(false), diploid(false), sample(0), phase(0), count(0) {}
  Haplotype(size_type sample_id, size_type phase_id) : offset(0), active(false), diploid(false), sample(sample_id), phase(phase_id), count(0) {}

  size_type size() const { return this->path.size(); }
  bool empty() const { return this->path.empty(); }

  void activate(size_type start_offset, bool is_diploid) { this->offset = start_offset; this->active = true; this->diploid = is_diploid; }
  void deactivate() { this->path = vector_type(); this->active = false; this->count++; }
};

//------------------------------------------------------------------------------

/*
  Paths corrensponding to the reference and the variants. Allele identifiers are 1-based;
  0 corresponds to the reference, which is not stored as an allele.

  Each site corresponds to a semiopen interval [start, end) on the reference path. When
  we append a variant to a haplotype, we first append reference until start, then append
  the variant path, and finally update the reference position to end.

  The intervals may overlap for nearby/adjacent sites. The variant paths may have
  flanking reference nodes that cause this.

  - If there are enough reference nodes at the end of the haplotype path and/or at the
    start of the variant path, we resolve the overlap by skipping these nodes.

  - Otherwise if report_overlap() returns true, we revert to the reference allele for
    the overlapping variant.

  - Otherwise we have a phase break.

  FILE FORMAT

  Version 1:
  - Proper header.
  - Sample and contig names.

  Version 0:
  - Original version without a header.
*/

class VariantPaths
{
public:
  typedef gbwt::size_type size_type;

//------------------------------------------------------------------------------

  VariantPaths();
  explicit VariantPaths(size_type reference_size);
  VariantPaths(const VariantPaths& source);
  VariantPaths(VariantPaths&& source);
  ~VariantPaths();

  void swap(VariantPaths& another);
  VariantPaths& operator=(const VariantPaths& source);
  VariantPaths& operator=(VariantPaths&& source);

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

//------------------------------------------------------------------------------

  // Header.
  std::uint32_t tag;
  std::uint32_t version;
  std::uint64_t flags;

  // Throws `sdsl::simple_sds::InvalidData` if the header is invalid.
  void check() const;

  void setVersion() { this->version = VERSION; }
  void set(std::uint64_t flag) { this->flags |= flag; }
  void unset(std::uint64_t flag) { this->flags &= ~flag; }
  bool get(std::uint64_t flag) const { return (this->flags & flag); }

  constexpr static std::uint32_t TAG = 0x6B37AA81;
  constexpr static std::uint32_t VERSION = Version::VARIANT_VERSION;

  constexpr static std::uint64_t FLAG_MASK         = 0x0003;
  constexpr static std::uint64_t FLAG_SAMPLE_NAMES = 0x0001;
  constexpr static std::uint64_t FLAG_CONTIG_NAME  = 0x0002;

//------------------------------------------------------------------------------

  /*
    reference      reference path
    ref_starts     starting positions of each site in the reference
    ref_ends       the positions following each site in the reference

    alt_paths      concatenated allele paths
    path_starts    pointers to the start of each path in alt_paths with a sentinel at the end
    site_starts    pointers to the start of each site in path_starts with a sentinel at the end

    phasing_files  file names for permanent PhasingInformation files
    file_offsets   the identifier of the first sample in each file
    file_counts    the number of samples in each file

    sample_names   names of the samples
    contig_name    name of the contig

    ref_index      the first occurrence of each node in the reference (not serialized)
  */
  vector_type              reference;
  std::vector<size_type>   ref_starts, ref_ends;
  vector_type              alt_paths;
  std::vector<size_type>   path_starts, site_starts;

  std::vector<std::string> phasing_files;
  std::vector<size_type>   file_offsets, file_counts;

  std::vector<std::string> sample_names;
  std::string              contig_name;

  std::unordered_map<node_type, node_type, size_type(*)(size_type)> ref_index;

//------------------------------------------------------------------------------

  // Queries.
  size_type size() const { return this->reference.size(); }
  size_type invalid_position() const { return this->size() + 1; }
  size_type paths() const { return this->path_starts.size() - 1; }
  size_type sites() const { return this->site_starts.size() - 1; }

  size_type alleles(size_type site) const { return this->site_starts[site + 1] - this->site_starts[site]; }

  size_type refPrev(size_type site) const { return (site > 0 ? this->refEnd(site - 1) : 0); }
  size_type refStart(size_type site) const { return this->ref_starts[site]; }
  size_type refEnd(size_type site) const { return this->ref_ends[site]; }

  node_type refAt(size_type offset) const { return this->reference[offset]; }
  vector_type getAllele(size_type site, size_type allele) const;

  size_type files() const { return this->phasing_files.size(); }
  const std::string& name(size_type file) const { return this->phasing_files[file]; }
  size_type offset(size_type file) const { return this->file_offsets[file]; }
  size_type count(size_type file) const { return this->file_counts[file]; }

  // Node width in bits. If both_orientations is set, check also the reverse of each node.
  size_type nodeWidth(bool both_orientations) const;

//------------------------------------------------------------------------------

  // Construction.
  void appendToReference(node_type node) { this->reference.push_back(node); }

  void indexReference();
  size_type firstOccurrence(node_type node);

  void addSite(size_type ref_start, size_type ref_end);
  void addAllele(const vector_type& path);

  void addFile(const std::string& filename, size_type sample_offset, size_type sample_count);

//------------------------------------------------------------------------------

  // Sample/contig names.
  bool hasSampleNames() const { return this->get(FLAG_SAMPLE_NAMES); }
  const std::vector<std::string>& getSampleNames() const { return this->sample_names; }
  void setSampleNames(const std::vector<std::string>& new_names);

  bool hasContigName() const { return this->get(FLAG_CONTIG_NAME); }
  const std::string& getContigName() const { return this->contig_name; }
  void setContigName(const std::string& new_name);

//------------------------------------------------------------------------------

  // Path generation.
  void appendReferenceUntil(Haplotype& haplotype, size_type site) const;
  void appendReferenceUntilEnd(Haplotype& haplotype) const;
  void appendVariant(Haplotype& haplotype, size_type site, size_type allele,
    std::function<void(const Haplotype&)> output, std::function<bool(size_type, size_type)> report_overlap) const;

private:
  void copy(const VariantPaths& source);
};

void checkOverlaps(const VariantPaths& variants, std::ostream& out, bool print_ids = false);

//------------------------------------------------------------------------------

struct Phasing
{
  size_type first, second;
  bool      diploid, phased;

  constexpr static size_type HAPLOID  = 0;
  constexpr static size_type UNPHASED = 1;
  constexpr static size_type PHASED   = 2;

  Phasing() {}

  explicit Phasing(size_type allele) : first(allele), second(0), diploid(false), phased(true) {}

  Phasing(size_type first_allele, size_type second_allele, bool is_phased = true) :
    first(first_allele), second(second_allele), diploid(true), phased(is_phased)
  {
  }

  Phasing(const std::string& genotype, bool was_diploid = true, bool phase_homozygous = true);

  void forcePhased(std::function<bool()> rng);

  size_type encode(size_type max_allele) const;
  void decode(size_type code, size_type max_allele);
  static size_type maxCode(size_type max_allele);

  bool operator== (const Phasing& another) const
  {
    return (this->first == another.first) && (this->second == another.second) && (this->diploid == another.diploid) && (this->phased == another.phased);
  }

  bool operator!= (const Phasing& another) const { return !(this->operator==(another)); }
};

std::ostream& operator<<(std::ostream& out, Phasing phasing);

//------------------------------------------------------------------------------

/*
  Phasing information for a number of samples in a temporary file.
*/

struct PhasingInformation
{
  // Header
  size_type sample_count, sample_offset; 
  size_type site_count;

  // File
  std::string                filename;
  sdsl::int_vector_buffer<8> data;
  bool                       temp_file;
  const static std::string   TEMP_FILE_PREFIX;  // "phasing"

  // Iterator
  size_type            site, data_offset;
  std::vector<Phasing> phasings;

  PhasingInformation();

  // Supported by a temporary file.
  PhasingInformation(size_type first_sample, size_type num_samples);

  // Supported by a permanent file [base_name]_[first_sample]_[num_samples].
  PhasingInformation(const std::string& base_name, size_type first_sample, size_type num_samples);

  // Use an existing permanent file.
  PhasingInformation(const VariantPaths& variants, size_type file);

  PhasingInformation(PhasingInformation&& source);
  ~PhasingInformation();

  PhasingInformation& operator= (PhasingInformation&& source);

  // Closing inactive files can save memory when there are many batches.
  void open();
  void close();
  bool isOpen() { return this->data.is_open(); }

  // The user may want to get the filename.
  const std::string& name() const { return this->filename; }
  bool isTemporary() const { return this->temp_file; }

  // If the structure is backed by a temporary file, make it persistent so that it
  // will not be deleted with the object. The file will still be deleted when the
  // program exits. No effect if the structure was backed by a permenent file.
  void makePersistent() { this->temp_file = false; }

  // Append the phasings for a new site.
  void append(const std::vector<Phasing>& new_site);

  // Iterate over the phasings.
  void begin() { this->site = 0; this->data_offset = 0; this->read(); }
  void operator++() { this->site++; this->read(); }
  const Phasing& operator[](size_type i) const { return this->phasings[i]; }
  size_type current() const { return this->site; }

  // Statistics.
  size_type size() const { return this->sample_count; }
  size_type offset() const { return this->sample_offset; }
  size_type limit() const { return this->offset() + this->size(); }
  size_type sites() const { return this->site_count; }
  size_type bytes() const { return this->data.size(); }

private:
  void read();

  PhasingInformation(const PhasingInformation&) = delete;
  PhasingInformation& operator= (const PhasingInformation&) = delete;
};

//------------------------------------------------------------------------------

// Opens and closes the phasings.
void generateHaplotypes(const VariantPaths& variants, PhasingInformation& phasings,
                        std::function<bool(size_type)> process_sample, std::function<void(const Haplotype&)> output,
                        std::function<bool(size_type, size_type)> report_overlap);

// As above, but use the PhasingInformation files stored in VariantPaths.
// Select a subset of files or use an empty set to indicate all files.
void generateHaplotypes(const VariantPaths& variants, const std::set<std::string>& phasing_files,
                        std::function<bool(size_type)> process_sample, std::function<void(const Haplotype&)> output,
                        std::function<bool(size_type, size_type)> report_overlap);

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_VARIANTS_H
