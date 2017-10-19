/*
  Copyright (c) 2017 Jouni Siren
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

#include <gbwt/dynamic_gbwt.h>

using namespace gbwt;

//------------------------------------------------------------------------------

const std::string tool_name = "GBWT merging";

void printUsage(int exit_code = EXIT_SUCCESS);

size_type insert(DynamicGBWT& index, const std::string input_name, size_type batch_size);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 5) { printUsage(); }

  size_type batch_size = DynamicGBWT::MERGE_BATCH_SIZE;
  std::string output;
  int c = 0;
  while((c = getopt(argc, argv, "b:o:")) != -1)
  {
    switch(c)
    {
    case 'b':
      batch_size = std::stoul(optarg); break;
    case 'o':
      output = optarg; break;
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }

  size_type input_files = argc - optind;
  size_type total_inserted = 0;
  std::string first_input = argv[optind]; optind++;
  if(input_files <= 1 || output.empty()) { printUsage(EXIT_FAILURE); }

  Version::print(std::cout, tool_name);

  printHeader("Input files"); std::cout << input_files << std::endl;
  printHeader("Output name"); std::cout << output << std::endl;
  printHeader("Batch size"); std::cout << batch_size << std::endl;
  std::cout << std::endl;

  double start = readTimer();

  DynamicGBWT index;
  sdsl::load_from_file(index, first_input + DynamicGBWT::EXTENSION);
  printStatistics(index, first_input);

  while(optind < argc)
  {
    std::string input_name = argv[optind];
    total_inserted += insert(index, input_name, batch_size);
    optind++;
  }

  sdsl::store_to_file(index, output + DynamicGBWT::EXTENSION);
  printStatistics(index, output);

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
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: merge_gbwt [options] -o output input1 input2 [input3 ...]" << std::endl;
  std::cerr << "  -b N  Use batches of N sequences for merging (default: " << DynamicGBWT::MERGE_BATCH_SIZE << ")" << std::endl;
  std::cerr << "  -o X  Use X as the base name for output (required)" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Use base names for the inputs and the output." << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

size_type
insert(DynamicGBWT& index, const std::string input_name, size_type batch_size)
{
  GBWT next;
  sdsl::load_from_file(next, input_name + GBWT::EXTENSION);
  printStatistics(next, input_name);
  index.merge(next, batch_size);
  return next.size();
}

//------------------------------------------------------------------------------
