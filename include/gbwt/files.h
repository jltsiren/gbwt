/*
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
  files.h: Public interface for file formats.
*/

//------------------------------------------------------------------------------

/*
  GBWT file header.

  Version 0:
  - Current version.
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

  const static std::uint32_t TAG = 0x6B376B37;
  const static std::uint32_t VERSION = 0;
  const static std::uint32_t MIN_VERSION = 0;

  GBWTHeader();

  size_type serialize(std::ostream& out, sdsl::structure_tree_node* v = nullptr, std::string name = "") const;
  void load(std::istream& in);
  bool check(std::uint32_t expected_version = VERSION) const;
  bool checkNew() const;

  void swap(GBWTHeader& another);

  bool operator==(const GBWTHeader& another) const;
  bool operator!=(const GBWTHeader& another) const { return !(this->operator==(another)); }
};

std::ostream& operator<<(std::ostream& stream, const GBWTHeader& header);

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_FILES_H
