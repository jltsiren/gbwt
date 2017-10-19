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

#include <limits>

#include <unistd.h>

#include <gbwt/files.h>

using namespace gbwt;

//------------------------------------------------------------------------------

const std::string tool_name = "Input transformation";

void printUsage(int exit_code = EXIT_SUCCESS);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 3) { printUsage(); }

  size_type max_sequences = std::numeric_limits<size_type>::max();
  int c = 0;
  while((c = getopt(argc, argv, "m:")) != -1)
  {
    switch(c)
    {
    case 'm':
      max_sequences = std::stoul(optarg); break;
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }
  if(optind + 1 >= argc) { printUsage(EXIT_FAILURE); }
  std::string input_name = argv[optind], output_name = argv[optind + 1];

  Version::print(std::cout, tool_name);

  printHeader("Input"); std::cout << input_name << std::endl;
  printHeader("Output"); std::cout << output_name << std::endl;
  printHeader("Max sequences"); std::cout << max_sequences << std::endl;
  std::cout << std::endl;

  double start = readTimer();

  // First pass: determine the largest node identifier.
  node_type max_node = 0;
  size_type total_length = 0;
  sdsl::int_vector_buffer<64> infile(input_name, std::ios::in, MEGABYTE, 64, true);
  for(node_type node : infile) { max_node = std::max(node, max_node); }
  if(infile.size() > 0 && infile[infile.size() - 1] != ENDMARKER)
  {
    std::cerr << "prepare_text: The text does not end with an endmarker" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // Second pass: transform the text.
  size_type sequences = 0;
  text_buffer_type outfile(output_name, std::ios::out, MEGABYTE, bit_length(max_node));
  for(node_type node : infile)
  {
    outfile.push_back(node); total_length++;
    if(node == ENDMARKER)
    {
      sequences++;
      if(sequences >= max_sequences) { break; }
    }
  }

  infile.close(); outfile.close();
  double seconds = readTimer() - start;

  printHeader("Text length"); std::cout << total_length << std::endl;
  printHeader("Alphabet size"); std::cout << (max_node + 1) << std::endl;
  printHeader("Sequences"); std::cout << sequences << std::endl;
  std::cout << std::endl;

  std::cout << "Processed " << total_length << " nodes in " << seconds << " seconds (" << (total_length / seconds) << " nodes/second)" << std::endl;
  std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: prepare_text input output" << std::endl;
  std::cerr << "  -m N  Read up to N sequences" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Transforms a sequence of 64-bit integers into the GBWT input format." << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

//------------------------------------------------------------------------------
