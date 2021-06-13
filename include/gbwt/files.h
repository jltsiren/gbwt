/*
  Copyright (c) 2017, 2018, 2019, 2021 Jouni Siren
  Copyright (c) 2015, 2016, 2017 Genome Research Ltd.

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

#ifndef GBWT_FILES_H
#define GBWT_FILES_H

#include <gbwt/utils.h>

namespace gbwt
{

/*
  files.h: File format headers.
*/

//------------------------------------------------------------------------------

/*
  GBWT file header.

  Version 5:
  - Uses metadata version 2.
  - Includes tags.
  - SDSL and simple-sds formats.
  - Compatible with versions 1 to 4.

  Version 4:
  - Uses metadata version 1.
  - Compatible with versions 1 to 3.

  Version 3:
  - Includes a flag for metadata.
  - Compatible with versions 1 and 2.

  Version 2:
  - Includes a flag for a bidirectional index.
  - Compatible with version 1.

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

  // Flag masks for old compatible versions.
  constexpr static std::uint32_t META2_VERSION      = 4;
  constexpr static std::uint64_t META2_FLAG_MASK    = 0x0003;

  constexpr static std::uint32_t META_VERSION       = 3;
  constexpr static std::uint64_t META_FLAG_MASK     = 0x0003;

  constexpr static std::uint32_t BD_VERSION         = 2;
  constexpr static std::uint64_t BD_FLAG_MASK       = 0x0001;

  constexpr static std::uint32_t OLD_VERSION        = 1;
  constexpr static std::uint64_t OLD_FLAG_MASK      = 0x0000;

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

  void swap(GBWTHeader& another);

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

  void swap(MetadataHeader& another);

  bool operator==(const MetadataHeader& another) const;
  bool operator!=(const MetadataHeader& another) const { return !(this->operator==(another)); }
};

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_FILES_H
