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

#include "gbwt.h"

using namespace gbwt;

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: build_gbwt base_name" << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  std::string base_name = argv[1];
  std::string gbwt_name = base_name + DynamicGBWT::EXTENSION;

  std::cout << "GBWT construction" << std::endl;
  std::cout << std::endl;

  printHeader("Base name"); std::cout << base_name << std::endl;
  std::cout << std::endl;

  double start = readTimer();

  text_type text;
  sdsl::load_from_file(text, base_name);

  DynamicGBWT gbwt;
  gbwt.insert(text);
  sdsl::store_to_file(gbwt, gbwt_name);

  double seconds = readTimer() - start;

  std::cout << gbwt.header << std::endl;
  std::cout << std::endl;

  std::cout << "Indexed " << gbwt.size() << " nodes in " << seconds << " seconds (" << (gbwt.size() / seconds) << " nodes/second)" << std::endl;
  std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------
