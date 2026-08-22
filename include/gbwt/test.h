/*
  Copyright (c) 2020 Jouni Siren

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

#ifndef GBWT_TEST_H
#define GBWT_TEST_H

#include "utils.h"

namespace gbwt
{

/*
  test.h: Utility functions used in various tests.
*/

//------------------------------------------------------------------------------

constexpr size_type DEFAULT_QUERIES      = 20000;
constexpr size_type DEFAULT_QUERY_LENGTH = 60;
constexpr size_type DEFAULT_RANDOM_SEED  = 0xDEADBEEF;

// Generate queries for find() from the input file.
std::vector<vector_type> generateQueries(const std::string& base_name, bool print, size_type queries = DEFAULT_QUERIES, size_type query_length = DEFAULT_QUERY_LENGTH, size_type random_seed = DEFAULT_RANDOM_SEED);

//------------------------------------------------------------------------------

// Find the starting offsets of the sequences in the input file.
std::vector<size_type> startOffsets(const std::string& base_name, bool print);

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_TEST_H
