#ifndef GBWT_FILES_H
#define GBWT_FILES_H

#include "utils.h"

namespace gbwt
{

/*
  files.h: File format headers.
*/

//------------------------------------------------------------------------------

/*
  GBWT file header.

  Version 6:
  - BWT is compressed using Zstandard.
  - Versions 1 to 5 can still be read.

  Version 5:
  - Uses metadata version 2.
  - Includes tags.
  - SDSL and simple-sds formats.
  - Versions 1 to 4 can still be read.

  Version 4:
  - Uses metadata version 1.
  - Versions 1 to 3 can still be read.

  Version 3:
  - Includes a flag for metadata.
  - Versions 1 and 2 can still be read.

  Version 2:
  - Includes a flag for a bidirectional index.
  - Version 1 can still be read.

  Version 1:
  - The first proper version.
  - Identical to version 0.

  Version 0:
  - Preliminary version.
*/

struct GBWTHeader
{
  typedef gbwt::size_type size_type;  // Needed for SDSL serialization.

  std::uint32_t tag;
  std::uint32_t version;
  std::uint64_t sequences;
  std::uint64_t size;           // Including the endmarkers.
  std::uint64_t offset;         // Range [1..offset] of the alphabet is empty.
  std::uint64_t alphabet_size;  // Largest node id + 1.
  std::uint64_t flags;

  constexpr static std::uint32_t TAG = 0x6B376B37;
  constexpr static std::uint32_t VERSION = Version::GBWT_VERSION;

  constexpr static std::uint64_t FLAG_MASK          = 0x0007;
  constexpr static std::uint64_t FLAG_BIDIRECTIONAL = 0x0001; // The index is bidirectional.
  constexpr static std::uint64_t FLAG_METADATA      = 0x0002; // The index contains metadata.
  constexpr static std::uint64_t FLAG_SIMPLE_SDS    = 0x0004; // simple-sds file format.

  // Symbolic names for versions that may be relevant when examining files,
  // even when the version is current or obsolete.
  constexpr static std::uint32_t ZSTD_VERSION = 6; // The current version with Zstandard compression.

  // Flag masks for old compatible versions.
  constexpr static std::uint32_t TAGS_VERSION = 5;
  constexpr static std::uint64_t TAGS_FLAG_MASK = 0x0007;

  constexpr static std::uint32_t META2_VERSION = 4;
  constexpr static std::uint64_t META2_FLAG_MASK = 0x0003;

  constexpr static std::uint32_t META_VERSION = 3;
  constexpr static std::uint64_t META_FLAG_MASK = 0x0003;

  constexpr static std::uint32_t BD_VERSION = 2;
  constexpr static std::uint64_t BD_FLAG_MASK = 0x0001;

  constexpr static std::uint32_t OLD_VERSION = 1;
  constexpr static std::uint64_t OLD_FLAG_MASK = 0x0000;

  constexpr static std::uint32_t MIN_VERSION = OLD_VERSION; // The oldest version we can read.
  constexpr static std::uint32_t MIN_SERIALIZE_VERSION = TAGS_VERSION; // The oldest version we can write.
  constexpr static std::uint32_t DEFAULT_VERSION = TAGS_VERSION; // The version we write by default.

  GBWTHeader();

  // simple-sds serialization sees the header as a serializable value.
  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  // Throws `sdsl::simple_sds::InvalidData` if the header is invalid.
  void check() const;

  void setVersion() { this->version = VERSION; }

  void set(std::uint64_t flag) { this->flags |= flag; }
  void unset(std::uint64_t flag) { this->flags &= ~flag; }
  bool get(std::uint64_t flag) const { return (this->flags & flag); }

  void swap(GBWTHeader& another) noexcept;

  bool operator==(const GBWTHeader& another) const;
  bool operator!=(const GBWTHeader& another) const { return !(this->operator==(another)); }
};

std::ostream& operator<<(std::ostream& stream, const GBWTHeader& header);

//------------------------------------------------------------------------------

/*
  Metadata structure header.

  Version 2:
  - Uses Dictionary based on StringArray.
  - Compatible with versions 0 to 1.

  Version 1:
  - Sample names, contig names, path names.
  - Compatible with version 0.

  Version 0:
  - Preliminary version with sample/haplotype/contig counts.
*/

struct MetadataHeader
{
  typedef gbwt::size_type size_type; // Needed for SDSL serialization.

  std::uint32_t tag;
  std::uint32_t version;
  std::uint64_t sample_count;
  std::uint64_t haplotype_count;
  std::uint64_t contig_count;
  std::uint64_t flags;

  constexpr static std::uint32_t TAG = 0x6B375E7A;
  constexpr static std::uint32_t VERSION = Version::METADATA_VERSION;

  constexpr static std::uint64_t FLAG_MASK         = 0x0007;
  constexpr static std::uint64_t FLAG_PATH_NAMES   = 0x0001;
  constexpr static std::uint64_t FLAG_SAMPLE_NAMES = 0x0002;
  constexpr static std::uint64_t FLAG_CONTIG_NAMES = 0x0004;

  // Flag masks for old compatible versions.
  constexpr static std::uint32_t NAMES_VERSION     = 1;
  constexpr static std::uint64_t NAMES_FLAG_MASK   = 0x0007;
  constexpr static std::uint32_t INITIAL_VERSION   = 0;
  constexpr static std::uint64_t INITIAL_FLAG_MASK = 0x0000;

  MetadataHeader();

  // simple-sds serialization sees the header as a serializable value.
  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);

  // Throws `sdsl::simple_sds::InvalidData` if the header is invalid.
  void check() const;
  void check_simple_sds() const;

  void setVersion() { this->version = VERSION; }

  void set(std::uint64_t flag) { this->flags |= flag; }
  void unset(std::uint64_t flag) { this->flags &= ~flag; }
  bool get(std::uint64_t flag) const { return (this->flags & flag); }

  void swap(MetadataHeader& another) noexcept;

  bool operator==(const MetadataHeader& another) const;
  bool operator!=(const MetadataHeader& another) const { return !(this->operator==(another)); }
};

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_FILES_H
