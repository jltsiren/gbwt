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

void printGBWT(const GBWT& gbwt, const std::string& name);
void printGBWT(const DynamicGBWT& gbwt, const std::string& name);

size_type insert(DynamicGBWT& index, const std::string input_name, size_type batch_size);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 3) { printUsage(); }

  size_type batch_size = DynamicGBWT::MERGE_BATCH_SIZE;
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
  if(optind + 1 >= argc) { printUsage(EXIT_FAILURE); }
  std::string first_input = argv[optind]; optind++;
  std::string output = argv[argc - 1];

  std::cout << "GBWT merging" << std::endl;
  std::cout << std::endl;

  printHeader("Output"); std::cout << output << std::endl;
  printHeader("Batch size"); std::cout << batch_size << std::endl;
  std::cout << std::endl;

  double start = readTimer();

  DynamicGBWT index;
  sdsl::load_from_file(index, first_input + DynamicGBWT::EXTENSION);
  printGBWT(index, first_input);

  size_type total_inserted = 0;
  while(optind + 1 < argc)
  {
    std::string input_name = argv[optind]; optind++;
    total_inserted += insert(index, input_name, batch_size);
  }

  sdsl::store_to_file(index, output + DynamicGBWT::EXTENSION);
  printGBWT(index, output);

  double seconds = readTimer() - start;

  std::cout << "Inserted " << total_inserted << " nodes in " << seconds << " seconds ("
            << (total_inserted / seconds) << " nodes/second)" << std::endl;
  std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  std::cerr << "Usage: merge_gbwt [options] input1 [input2 ...] output" << std::endl;
  std::cerr << "  -b N  Use batches of N sequences for merging (default "
            << DynamicGBWT::MERGE_BATCH_SIZE << ")" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Use base names for the inputs and the output. Using compressed GBWTs from input2" << std::endl;
  std::cerr << "onwards saves memory but is slower." << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

void
printGBWT(const GBWT& gbwt, const std::string& name)
{
  printHeader("Compressed GBWT"); std::cout << name << std::endl;
  printHeader("Total length"); std::cout << gbwt.size() << std::endl;
  printHeader("Sequences"); std::cout << gbwt.sequences() << std::endl;
  printHeader("Alphabet size"); std::cout << gbwt.sigma() << std::endl;
  printHeader("Effective"); std::cout << gbwt.effective() << std::endl;
  printHeader("Runs"); std::cout << gbwt.runs() << std::endl;
  std::cout << std::endl;
}

void
printGBWT(const DynamicGBWT& gbwt, const std::string& name)
{
  printHeader("Dynamic GBWT"); std::cout << name << std::endl;
  printHeader("Total length"); std::cout << gbwt.size() << std::endl;
  printHeader("Sequences"); std::cout << gbwt.sequences() << std::endl;
  printHeader("Alphabet size"); std::cout << gbwt.sigma() << std::endl;
  printHeader("Effective"); std::cout << gbwt.effective() << std::endl;
  printHeader("Runs"); std::cout << gbwt.runs() << std::endl;
  std::cout << std::endl;
}

size_type
insert(DynamicGBWT& index, const std::string input_name, size_type batch_size)
{
  GBWT next;
  sdsl::load_from_file(next, input_name + GBWT::EXTENSION);
  printGBWT(next, input_name);
  index.merge(next, batch_size);
  return next.size();
}

//------------------------------------------------------------------------------
