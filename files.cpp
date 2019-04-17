/*
  Copyright (c) 2017, 2018, 2019 Jouni Sirén
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
constexpr std::uint32_t GBWTHeader::META_VERSION;
constexpr std::uint64_t GBWTHeader::META_FLAG_MASK;
constexpr std::uint32_t GBWTHeader::BD_VERSION;
constexpr std::uint64_t GBWTHeader::BD_FLAG_MASK;
constexpr std::uint32_t GBWTHeader::OLD_VERSION;
constexpr std::uint64_t GBWTHeader::OLD_FLAG_MASK;

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

bool
GBWTHeader::check() const
{
  if(this->tag != TAG) { return false; }
  switch(this->version)
  {
  case VERSION:
    return ((this->flags & FLAG_MASK) == this->flags);
  case META_VERSION:
    return ((this->flags & META_FLAG_MASK) == this->flags);
  case BD_VERSION:
    return ((this->flags & BD_FLAG_MASK) == this->flags);
  case OLD_VERSION:
    return ((this->flags & OLD_FLAG_MASK) == this->flags);
  default:
    return false;
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

} // namespace gbwt
