/*
  Copyright (c) 2017, 2018, 2019, 2021 Jouni Sirén
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

#include <gbwt/files.h>

namespace gbwt
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr std::uint32_t GBWTHeader::TAG;
constexpr std::uint32_t GBWTHeader::VERSION;
constexpr std::uint64_t GBWTHeader::FLAG_MASK;
constexpr std::uint64_t GBWTHeader::FLAG_BIDIRECTIONAL;
constexpr std::uint64_t GBWTHeader::FLAG_METADATA;
constexpr std::uint64_t GBWTHeader::FLAG_SIMPLE_SDS;
constexpr std::uint32_t GBWTHeader::META2_VERSION;
constexpr std::uint64_t GBWTHeader::META2_FLAG_MASK;
constexpr std::uint32_t GBWTHeader::META_VERSION;
constexpr std::uint64_t GBWTHeader::META_FLAG_MASK;
constexpr std::uint32_t GBWTHeader::BD_VERSION;
constexpr std::uint64_t GBWTHeader::BD_FLAG_MASK;
constexpr std::uint32_t GBWTHeader::OLD_VERSION;
constexpr std::uint64_t GBWTHeader::OLD_FLAG_MASK;

constexpr std::uint32_t MetadataHeader::TAG;
constexpr std::uint32_t MetadataHeader::VERSION;
constexpr std::uint64_t MetadataHeader::FLAG_MASK;
constexpr std::uint64_t MetadataHeader::FLAG_PATH_NAMES;
constexpr std::uint64_t MetadataHeader::FLAG_SAMPLE_NAMES;
constexpr std::uint64_t MetadataHeader::FLAG_CONTIG_NAMES;
constexpr std::uint32_t MetadataHeader::NAMES_VERSION;
constexpr std::uint64_t MetadataHeader::NAMES_FLAG_MASK;
constexpr std::uint32_t MetadataHeader::INITIAL_VERSION;
constexpr std::uint64_t MetadataHeader::INITIAL_FLAG_MASK;

//------------------------------------------------------------------------------

/*
  An empty index is bidirectional.
*/
GBWTHeader::GBWTHeader() :
  tag(TAG), version(VERSION),
  sequences(0), size(0),
  offset(0), alphabet_size(0),
  flags(FLAG_BIDIRECTIONAL)
{
}

size_type
GBWTHeader::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;
  written_bytes += sdsl::write_member(this->tag, out, child, "tag");
  written_bytes += sdsl::write_member(this->version, out, child, "version");
  written_bytes += sdsl::write_member(this->sequences, out, child, "sequences");
  written_bytes += sdsl::write_member(this->size, out, child, "size");
  written_bytes += sdsl::write_member(this->offset, out, child, "offset");
  written_bytes += sdsl::write_member(this->alphabet_size, out, child, "alphabet_size");
  written_bytes += sdsl::write_member(this->flags, out, child, "flags");
  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
GBWTHeader::load(std::istream& in)
{
  sdsl::read_member(this->tag, in);
  sdsl::read_member(this->version, in);
  sdsl::read_member(this->sequences, in);
  sdsl::read_member(this->size, in);
  sdsl::read_member(this->offset, in);
  sdsl::read_member(this->alphabet_size, in);
  sdsl::read_member(this->flags, in);
}

void
GBWTHeader::check() const
{
  if(this->tag != TAG)
  {
    throw sdsl::simple_sds::InvalidData("GBWTHeader: Invalid tag");
  }

  if(this->version > VERSION || this->version < OLD_VERSION)
  {
    std::string msg = "GBWTHeader: Expected v" + std::to_string(OLD_VERSION) + " to v" + std::to_string(VERSION) + ", got v" + std::to_string(this->version);
    throw sdsl::simple_sds::InvalidData(msg);
  }

  std::uint64_t mask = 0;
  switch(this->version)
  {
  case VERSION:
    mask = FLAG_MASK; break;
  case META2_VERSION:
    mask = META2_FLAG_MASK; break;
  case META_VERSION:
    mask = META_FLAG_MASK; break;
  case BD_VERSION:
    mask = BD_FLAG_MASK; break;
  case OLD_VERSION:
    mask = OLD_FLAG_MASK; break;
  }
  if((this->flags & mask) != this->flags)
  {
    throw sdsl::simple_sds::InvalidData("GBWTHeader: Invalid flags");
  }
}

void
GBWTHeader::swap(GBWTHeader& another)
{
  if(this != &another)
  {
    std::swap(this->tag, another.tag);
    std::swap(this->version, another.version);
    std::swap(this->sequences, another.sequences);
    std::swap(this->size, another.size);
    std::swap(this->offset, another.offset);
    std::swap(this->alphabet_size, another.alphabet_size);
    std::swap(this->flags, another.flags);
  }
}

bool
GBWTHeader::operator==(const GBWTHeader& another) const
{
  return (this->tag == another.tag &&
          this->version == another.version &&
          this->sequences == another.sequences &&
          this->size == another.size &&
          this->offset == another.offset &&
          this->alphabet_size == another.alphabet_size &&
          this->flags == another.flags);
}

std::ostream& operator<<(std::ostream& stream, const GBWTHeader& header)
{
  return stream << "GBWT v" << header.version << ": "
                << header.sequences << " sequences of total length " << header.size
                << ", alphabet size " << header.alphabet_size << " with offset " << header.offset;
}

//------------------------------------------------------------------------------

MetadataHeader::MetadataHeader() :
  tag(TAG), version(VERSION),
  sample_count(0), haplotype_count(0), contig_count(0),
  flags(0)
{
}

size_type
MetadataHeader::serialize(std::ostream& out, sdsl::structure_tree_node* v, std::string name) const
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(*this));
  size_type written_bytes = 0;

  written_bytes += sdsl::write_member(this->tag, out, child, "tag");
  written_bytes += sdsl::write_member(this->version, out, child, "version");
  written_bytes += sdsl::write_member(this->sample_count, out, child, "sample_count");
  written_bytes += sdsl::write_member(this->haplotype_count, out, child, "haplotype_count");
  written_bytes += sdsl::write_member(this->contig_count, out, child, "contig_count");
  written_bytes += sdsl::write_member(this->flags, out, child, "flags");

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

void
MetadataHeader::load(std::istream& in)
{
  sdsl::read_member(this->tag, in);
  sdsl::read_member(this->version, in);
  sdsl::read_member(this->sample_count, in);
  sdsl::read_member(this->haplotype_count, in);
  sdsl::read_member(this->contig_count, in);
  sdsl::read_member(this->flags, in);
}

void
MetadataHeader::check() const
{
  if(this->tag != TAG)
  {
    throw sdsl::simple_sds::InvalidData("MetadataHeader: Invalid tag");
  }

  if(this->version > VERSION || this->version < INITIAL_VERSION)
  {
    std::string msg = "MetadataHeader: Expected v" + std::to_string(INITIAL_VERSION) + " to v" + std::to_string(VERSION) + ", got v" + std::to_string(this->version);
    throw sdsl::simple_sds::InvalidData(msg);
  }

  std::uint64_t mask = 0;
  switch(this->version)
  {
  case VERSION:
    mask = FLAG_MASK; break;
  case NAMES_VERSION:
    mask = NAMES_FLAG_MASK; break;
  case INITIAL_VERSION:
    mask = INITIAL_FLAG_MASK; break;
  }
  if((this->flags & mask) != this->flags)
  {
    throw sdsl::simple_sds::InvalidData("MetadataHeader: Invalid flags");
  }
}

void
MetadataHeader::check_simple_sds() const
{
  if(this->tag != TAG)
  {
    throw sdsl::simple_sds::InvalidData("MetadataHeader: Invalid tag");
  }

  if(this->version != VERSION)
  {
    std::string msg = "MetadataHeader: Expected v" + std::to_string(VERSION) + ", got v" + std::to_string(this->version);
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
    throw sdsl::simple_sds::InvalidData("MetadataHeader: Invalid flags");
  }
}

void
MetadataHeader::swap(MetadataHeader& another)
{
  if(this != &another)
  {
    std::swap(this->tag, another.tag);
    std::swap(this->version, another.version);
    std::swap(this->sample_count, another.sample_count);
    std::swap(this->haplotype_count, another.haplotype_count);
    std::swap(this->contig_count, another.contig_count);
    std::swap(this->flags, another.flags);
  }
}

bool
MetadataHeader::operator==(const MetadataHeader& another) const
{
  return (this->tag == another.tag &&
          this->version == another.version &&
          this->sample_count == another.sample_count &&
          this->haplotype_count == another.haplotype_count &&
          this->contig_count == another.contig_count &&
          this->flags == another.flags);
}

//------------------------------------------------------------------------------

} // namespace gbwt
