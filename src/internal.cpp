/*
  Copyright (c) 2017, 2018, 2019 Jouni Siren
  Copyright (c) 2017 Genome Research Ltd.

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

namespace gbwt
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_type DiskIO::block_size;

constexpr size_type ByteCode::DATA_BITS;
constexpr ByteCode::code_type ByteCode::DATA_MASK;
constexpr ByteCode::code_type ByteCode::NEXT_BYTE;

//------------------------------------------------------------------------------

template<>
size_type
serializeVector<std::string>(const std::vector<std::string>& data, std::ostream& out, sdsl::structure_tree_node* v, std::string name)
{
  sdsl::structure_tree_node* child = sdsl::structure_tree::add_child(v, name, sdsl::util::class_name(data));
  size_type written_bytes = 0;

  size_type data_size = data.size();
  written_bytes += sdsl::write_member(data_size, out, child, "size");

  for(const std::string& value : data)
  {
    written_bytes += sdsl::write_member(value, out, child, "data");
  }

  sdsl::structure_tree::add_size(child, written_bytes);
  return written_bytes;
}

template<>
void
loadVector<std::string>(std::vector<std::string>& data, std::istream& in)
{
  size_type data_size = 0;
  sdsl::read_member(data_size, in);

  data.resize(data_size);
  for(std::string& value : data)
  {
    sdsl::read_member(value, in);
  }
}

//------------------------------------------------------------------------------

Run::Run(value_type alphabet_size) :
  sigma(alphabet_size),
  run_continues(0)
{
  value_type max_code = std::numeric_limits<code_type>::max();
  if(this->sigma < max_code)
  {
    this->run_continues = (max_code + 1) / this->sigma;
  }
}

//------------------------------------------------------------------------------

Sequence::Sequence() :
  id(0), curr(ENDMARKER), next(ENDMARKER), offset(0), pos(0)
{
}

Sequence::Sequence(const text_type& text, size_type i, size_type seq_id) :
  id(seq_id), curr(ENDMARKER), next(text[i]), offset(seq_id), pos(i)
{
}

Sequence::Sequence(const vector_type& text, size_type i, size_type seq_id) :
  id(seq_id), curr(ENDMARKER), next(text[i]), offset(seq_id), pos(i)
{
}

Sequence::Sequence(node_type node, size_type seq_id, size_type source_pos) :
  id(seq_id), curr(ENDMARKER), next(node), offset(seq_id), pos(source_pos)
{
}

//------------------------------------------------------------------------------

} // namespace gbwt
