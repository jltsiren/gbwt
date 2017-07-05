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

#include <unistd.h>

#include "gbwt.h"

using namespace gbwt;

//------------------------------------------------------------------------------

void printUsage(int exit_code = EXIT_SUCCESS);

double build(DynamicGBWT& gbwt, std::string& base_name, size_type batch_size);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2) { printUsage(); }

  size_type batch_size = 0;
  int c = 0;
  while((c = getopt(argc, argv, "b:")) != -1)
  {
    switch(c)
    {
    case 'b':
      batch_size = std::stoul(optarg); break;
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }
  if(optind >= argc) { printUsage(EXIT_FAILURE); }
  std::string base_name = argv[optind];

  std::cout << "GBWT construction" << std::endl;
  std::cout << std::endl;

  printHeader("Base name"); std::cout << base_name << std::endl;
  if(batch_size != 0) { printHeader("Batch size"); std::cout << batch_size << " million" << std::endl; }
  std::cout << std::endl;

  DynamicGBWT gbwt;
  double seconds = build(gbwt, base_name, batch_size);

  printHeader("Total length"); std::cout << gbwt.size() << std::endl;
  printHeader("Sequences"); std::cout << gbwt.sequences() << std::endl;
  printHeader("Alphabet size"); std::cout << gbwt.sigma() << std::endl;
  printHeader("Effective"); std::cout << gbwt.effective() << std::endl;
  printHeader("Runs"); std::cout << gbwt.runs() << std::endl;
  std::cout << std::endl;

  std::cout << "Indexed " << gbwt.size() << " nodes in " << seconds << " seconds (" << (gbwt.size() / seconds) << " nodes/second)" << std::endl;
  std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  std::cerr << "Usage: build_gbwt [options] base_name" << std::endl;
  std::cerr << "  -b N  Insert in batches of N million nodes" << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

//------------------------------------------------------------------------------

double
build(DynamicGBWT& gbwt, std::string& base_name, size_type batch_size)
{
  double start = readTimer();

  text_buffer_type input(base_name);
  for(size_type i = 0; i < input.size(); )
  {
    size_type limit = (batch_size == 0 ? input.size() : std::min(input.size(), i + batch_size * MILLION));
    while(limit > i)
    {
      if(input[limit - 1] == ENDMARKER) { break; }
      limit--;
    }
    if(limit <= i)
    {
      std::cerr << "build_gbwt: Cannot find an endmarker in the batch starting from offset " << i << std::endl;
      std::exit(EXIT_FAILURE);
    }
    text_type batch(limit - i, 0, input.width());
    for(size_type j = i; j < limit; j++) { batch[j - i] = input[j]; }
    gbwt.insert(batch);
    i = limit;
  }

  std::string gbwt_name = base_name + DynamicGBWT::EXTENSION;
  sdsl::store_to_file(gbwt, gbwt_name);

  return readTimer() - start;
}

//------------------------------------------------------------------------------
