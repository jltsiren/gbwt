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

int
main(int argc, char** argv)
{
  if(argc < 3)
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
  sdsl::construct_bwt<0>(config); // FIXME should use 0 for all terminators and ignore the global terminator
  std::cout << "BWT built" << std::endl;

  std::cout << std::endl;
  double seconds = readTimer() - start;
  size_type data_size = header.sequences + header.total_length;

  std::cout << "Indexed " << data_size << " nodes in " << seconds << " (" << (data_size / seconds) << " nodes/second)" << std::endl;
  std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------
