/*
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

#include "internal.h"

namespace gbwt
{

//------------------------------------------------------------------------------

Run::Run(size_type alphabet_size) :
  sigma(alphabet_size),
  run_continues(0)
{
  size_type max_code = std::numeric_limits<code_type>::max();
  if(this->sigma < max_code)
  {
    this->run_continues = (max_code + 1) / this->sigma;
  }
}

//------------------------------------------------------------------------------

void
RunMerger::flush()
{
  if(this->accumulator.second > 0)
  {
    this->runs.push_back(this->accumulator);
    this->accumulator.second = 0;
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

//------------------------------------------------------------------------------

} // namespace gbwt
