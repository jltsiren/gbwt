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

#include <sdsl/construct.hpp>

#include "files.h"

using namespace gbwt;

//------------------------------------------------------------------------------

size_type buildBWT(sdsl::cache_config& config, GBWTHeader& header);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: index_text base_name" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Builds the suffix array and the BWT." << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  std::string base_name = argv[1];
  std::string header_name = base_name + HEADER_EXTENSION;
  std::string text_name = base_name + TEXT_EXTENSION;
  std::string alphabet_name = base_name + ALPHABET_EXTENSION;
  std::string document_name = base_name + DOCUMENT_EXTENSION;
  std::string sa_name = base_name + SA_EXTENSION;
  std::string bwt_name = base_name + BWT_EXTENSION;

  std::cout << "Indexing the text" << std::endl;
  std::cout << std::endl;

  printHeader("Base name"); std::cout << base_name << std::endl;
  std::cout << std::endl;

  double start = readTimer();

  GBWTHeader header;
  sdsl::load_from_file(header, header_name);
  std::cout << header << std::endl;
  std::cout << std::endl;

  sdsl::cache_config config(false);
  config.file_map[sdsl::conf::KEY_TEXT_INT] = text_name;

  config.file_map[sdsl::conf::KEY_SA] = sa_name;
  sdsl::construct_sa<0>(config);
  std::cout << "Suffix array built" << std::endl;

  config.file_map[sdsl::conf::KEY_BWT_INT] = bwt_name;
  size_type runs = buildBWT(config, header);
  std::cout << "BWT built, " << runs << " runs" << std::endl;

  std::cout << std::endl;
  double seconds = readTimer() - start;
  size_type data_size = header.sequences + header.total_length;

  std::cout << "Indexed " << data_size << " nodes in " << seconds << " seconds (" << (data_size / seconds) << " nodes/second)" << std::endl;
  std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

size_type
buildBWT(sdsl::cache_config& config, GBWTHeader& header)
{
  text_type text;
  sdsl::load_from_cache(text, sdsl::conf::KEY_TEXT_INT, config);

  sdsl::int_vector_buffer<0> sa(cache_file_name(sdsl::conf::KEY_SA, config));
  text_buffer_type bwt(cache_file_name(sdsl::conf::KEY_BWT_INT, config), std::ios::out, MEGABYTE, bit_length(header.alphabet_size - 1));

  // Return to the original alphabet, where each terminator is 0.
  // We skip the global terminator at SA[0].
  size_type to_add[2] = { ~(size_type)0, text.size() - 1 };
  size_type runs = 0; value_type prev = ~(value_type)0;
  for(size_type i = 1; i < text.size(); i++)
  {
    value_type value = text[sa[i] + to_add[sa[i] == 0]];
    value = (value > header.sequences ? value - header.sequences : 0);
    if(value != prev) { runs++; prev = value; }
    bwt.push_back(value);
  }

  sa.close(); bwt.close();
  return runs;
}

//------------------------------------------------------------------------------
