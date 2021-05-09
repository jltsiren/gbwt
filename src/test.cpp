/*
  Copyright (c) 2018, 2020 Jouni Siren
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

#include <algorithm>
#include <random>

#include <gbwt/test.h>

namespace gbwt
{

//------------------------------------------------------------------------------

std::vector<vector_type>
generateQueries(const std::string& base_name, bool print, size_type queries, size_type query_length, size_type random_seed)
{
  std::mt19937_64 rng(random_seed);
  std::vector<vector_type> result;
  text_buffer_type text(base_name);
  if(text.size() <= query_length)
  {
    std::cerr << "generateQueries(): Text length " << text.size() << " is too small for query length " << query_length << std::endl;
    return result;
  }

  size_type attempts = 0;
  while(result.size() < queries && attempts < 2 * queries)
  {
    size_type start_offset = rng() % (text.size() - query_length);  // We do not want queries containing the final endmarker.
    vector_type candidate(text.begin() + start_offset, text.begin() + start_offset + query_length);
    if(std::find(candidate.begin(), candidate.end(), ENDMARKER) == candidate.end())
    {
      result.push_back(candidate);
    }
    attempts++;
  }

  if(print)
  {
    std::cout << "Generated " << result.size() << " queries of total length " << (result.size() * query_length) << std::endl;
  }
  return result;
}

//------------------------------------------------------------------------------

std::vector<size_type>
startOffsets(const std::string& base_name, bool print)
{
  std::vector<size_type> offsets;
  text_buffer_type text(base_name);
  bool seq_start = true;
  for(size_type i = 0; i < text.size(); i++)
  {
    if(seq_start) { offsets.push_back(i); seq_start = false; }
    if(text[i] == ENDMARKER) { seq_start = true; }
  }

  if(print)
  {
    std::cout << "Found the starting offsets for " << offsets.size() << " sequences" << std::endl;
  }
  return offsets;
}

//------------------------------------------------------------------------------

} // namespace gbwt
