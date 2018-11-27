/*
  Copyright (c) 2018 Jouni Siren

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

#include <string>
#include <unistd.h>

#include <gbwt/dynamic_gbwt.h>

using namespace gbwt;

//------------------------------------------------------------------------------

const std::string tool_name = "GBWT sequence remove";

void printUsage(int exit_code = EXIT_SUCCESS);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 3) { printUsage(); }

  int c = 0;
  std::string output;
  size_type chunk_size = DynamicGBWT::REMOVE_CHUNK_SIZE;
  bool range = false;
  while((c = getopt(argc, argv, "c:o:r")) != -1)
  {
    switch(c)
    {
    case 'c':
      chunk_size = std::max(1ul, std::stoul(optarg));
      break;
    case 'o':
      output = optarg; break;
    case 'r':
      range = true; break;
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }

  if(optind + 1 >= argc) { printUsage(EXIT_FAILURE); }
  std::string base_name = argv[optind]; optind++;
  if(output.empty()) { output = base_name; }
  std::vector<size_type> seq_ids;
  if(range)
  {
    if(argc != optind + 2) { printUsage(EXIT_FAILURE); }
    size_type start = std::stoul(argv[optind]); optind++;
    size_type stop = std::stoul(argv[optind]); optind++;
    if(stop < start) { printUsage(EXIT_FAILURE); }
    for(size_type seq_id = start; seq_id <= stop; seq_id++)
    {
      seq_ids.push_back(seq_id);
    }
  }
  else
  {
    while(optind < argc)
    {
      seq_ids.push_back(std::stoul(argv[optind]));
      optind++;
    }
  }

  Version::print(std::cout, tool_name);

  printHeader("Input"); std::cout << base_name << std::endl;
  printHeader("Output"); std::cout << output << std::endl;
  printHeader("Sequences"); std::cout << seq_ids.size() << std::endl;
  printHeader("Chunk size"); std::cout << chunk_size << std::endl;
  std::cout << std::endl;

  double start = readTimer();

  DynamicGBWT index;
  sdsl::load_from_file(index, base_name + DynamicGBWT::EXTENSION);
  printStatistics(index, base_name);

  size_type total_length = index.remove(seq_ids, chunk_size);
  sdsl::store_to_file(index, output + DynamicGBWT::EXTENSION);
  printStatistics(index, output);

  double seconds = readTimer() - start;

  std::cout << "Removed " << total_length << " nodes in " << seconds << " seconds ("
            << (total_length / seconds) << " nodes/second)" << std::endl;
  std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: remove_seq [options] base_name seq1 [seq2 ...]" << std::endl;
  std::cerr << std::endl;
  std::cerr << "  -c N  Build the R in chunks of N sequences per thread (default: " << DynamicGBWT::REMOVE_CHUNK_SIZE << ")" << std::endl;
  std::cerr << "  -o X  Use X as the base name for output" << std::endl;
  std::cerr << "  -r    Remove a range of sequences (inclusive; requires 2 sequence ids)" << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

//------------------------------------------------------------------------------
