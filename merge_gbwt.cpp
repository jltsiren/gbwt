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

#include "dynamic_gbwt.h"

using namespace gbwt;

//------------------------------------------------------------------------------

void printUsage(int exit_code = EXIT_SUCCESS);

void printGBWT(const DynamicGBWT& gbwt, const std::string& name);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 4) { printUsage(); }

  size_type batch_size = DynamicGBWT::BATCH_SIZE;
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
  if(optind + 2 >= argc) { printUsage(EXIT_FAILURE); }
  std::string input1 = argv[optind], input2 = argv[optind + 1], output = argv[optind + 2];

  std::cout << "GBWT merging" << std::endl;
  std::cout << std::endl;

  printHeader("Input 1"); std::cout << input1 << std::endl;
  printHeader("Input 2"); std::cout << input2 << std::endl;
  printHeader("Output"); std::cout << output << std::endl;
  printHeader("Batch size"); std::cout << batch_size << std::endl;
  std::cout << std::endl;

  double start = readTimer();

  DynamicGBWT left;
  sdsl::load_from_file(left, input1 + DynamicGBWT::EXTENSION);
  printGBWT(left, input1);

  DynamicGBWT right;
  sdsl::load_from_file(right, input2 + DynamicGBWT::EXTENSION);
  printGBWT(right, input2);

  left.merge(right, batch_size);
  sdsl::store_to_file(left, output + DynamicGBWT::EXTENSION);
  printGBWT(left, output);

  double seconds = readTimer() - start;

  std::cout << "Inserted " << right.size() << " nodes in " << seconds << " seconds (" << (right.size() / seconds) << " nodes/second)" << std::endl;
  std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  std::cerr << "Usage: merge_gbwt input1 input2 output" << std::endl;
  std::cerr << "  -b N  Use batches of N sequences for merging (default " << DynamicGBWT::BATCH_SIZE << ")" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Use base names for the inputs and the output." << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

void
printGBWT(const DynamicGBWT& gbwt, const std::string& name)
{
  printHeader("GBWT"); std::cout << name << std::endl;
  printHeader("Total length"); std::cout << gbwt.size() << std::endl;
  printHeader("Sequences"); std::cout << gbwt.sequences() << std::endl;
  printHeader("Alphabet size"); std::cout << gbwt.sigma() << std::endl;
  printHeader("Effective"); std::cout << gbwt.effective() << std::endl;
  printHeader("Runs"); std::cout << gbwt.runs() << std::endl;
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
