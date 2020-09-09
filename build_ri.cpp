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

#include <unistd.h>

#include <gbwt/fast_locate.h>

using namespace gbwt;

//------------------------------------------------------------------------------

const std::string tool_name = "R-index construction";

void printUsage(int exit_code = EXIT_SUCCESS);

// FIXME index verification

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2) { printUsage(); }

  std::string base_name;
  bool verify_index = false;
  int c = 0;
  while((c = getopt(argc, argv, "v")) != -1)
  {
    switch(c)
    {
    case 'v':
      verify_index = true; break;
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }
  if(optind + 1 != argc)
  {
    printUsage(EXIT_FAILURE);
  }
  base_name = argv[optind];

  Version::print(std::cout, tool_name);
  printHeader("Base name"); std::cout << base_name << std::endl;
  std::cout << std::endl;

  double start = readTimer();

  std::string gbwt_name = base_name + GBWT::EXTENSION;
  GBWT gbwt_index;
  if(!sdsl::load_from_file(gbwt_index, gbwt_name))
  {
    std::cerr << "build_ri: Cannot load the GBWT index from " << gbwt_name << std::endl;
    std::exit(EXIT_FAILURE);
  }
  printStatistics(gbwt_index, base_name);

  FastLocate r_index(gbwt_index);
  std::string r_index_name = base_name + FastLocate::EXTENSION;
  if(!sdsl::store_to_file(r_index, r_index_name))
  {
    std::cerr << "build_ri: Cannot write the r-index to " << r_index_name << std::endl;
    std::exit(EXIT_FAILURE);
  }
  printStatistics(r_index, base_name);

  double seconds = readTimer() - start;

  std::cout << "Indexed " << gbwt_index.size() << " nodes in " << seconds << " seconds (" << (gbwt_index.size() / seconds) << " nodes/second)" << std::endl;
  std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  if(verify_index)
  {
    std::cerr << "build_ri: Verification has not been implemented yet" << std::endl;
  }

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: build_ri [options] base_name" << std::endl;
  std::cerr << "  -v    Verify the index after construction" << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

//------------------------------------------------------------------------------
