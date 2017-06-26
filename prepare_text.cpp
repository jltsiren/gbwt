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

#include <limits>

#include <unistd.h>

#include "files.h"

using namespace gbwt;

//------------------------------------------------------------------------------

void printUsage(int exit_code = EXIT_SUCCESS);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2) { printUsage(); }

  size_type max_sequences = std::numeric_limits<size_type>::max();
  int c = 0;
  std::string input_name;
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
  if(optind >= argc) { printUsage(EXIT_FAILURE); }

  std::string base_name = argv[optind];
  std::string text_name = base_name + TEXT_EXTENSION;
  std::string header_name = base_name + HEADER_EXTENSION;
  std::string alphabet_name = base_name + ALPHABET_EXTENSION;
  std::string document_name = base_name + DOCUMENT_EXTENSION;

  std::cout << "Preparing the text for indexing" << std::endl;
  std::cout << std::endl;

  if(!(input_name.empty())) { printHeader("Input"); std::cout << input_name << std::endl; }
  printHeader("Base name"); std::cout << base_name << std::endl;
  printHeader("Max sequences"); std::cout << max_sequences << std::endl;
  std::cout << std::endl;

  double start = readTimer();

  // Pass 1: Determine data size, alphabet size, and the number of sequences.
  GBWTHeader header;
  text_buffer_type text(text_name);
  std::vector<size_type> terminators;
  std::cout << "First pass: data size, alphabet size, number of sequences" << std::endl;
  for(auto iter = text.begin(); iter != text.end(); ++iter)
  {
    node_type value = *iter;
    header.alphabet_size = std::max(value + 1, header.alphabet_size);
    if(value == 0)
    {
      header.sequences++; terminators.push_back(iter - text.begin());
      if(header.sequences >= max_sequences) { break; }
    }
    else { header.total_length++; }
  }
  std::cout << header << std::endl;
  std::cout << std::endl;

  // Sanity checks.
  if(header.sequences == 0 || header.total_length == 0)
  {
    std::cerr << "prepare_text: The input is empty" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  if(header.alphabet_size > MAX_NODES)
  {
    std::cerr << "prepare_text: The alphabet is too large" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  for(size_type i = 1; i < terminators.size(); i++)
  {
    if(terminators[i] == terminators[i - 1] + 1)
    {
      std::cerr << "prepare_text: One of the sequences is empty" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  if(terminators.back() + 1 != header.sequences + header.total_length)
  {
    std::cerr << "prepare_text: The input does not end with a terminator" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  sdsl::sd_vector<> document_borders(terminators.begin(), terminators.end());
  sdsl::store_to_file(document_borders, document_name);

  // Pass 2: Determine the node identifiers that are present.
  std::vector<size_type> counts(header.alphabet_size, 0);
  size_type current = 0;
  std::cout << "Second pass: active alphabet" << std::endl;
  for(node_type value : text)
  {
    counts[value]++;
    if(value == 0)
    {
      current++;
      if(current >= max_sequences) { break; }
    }
  }
  text.close();

  for(size_type i = 1; i < counts.size(); i++) { if(counts[i] > 0) { header.nodes++; } }
  sdsl::store_to_file(header, header_name);

  std::cout << header << std::endl;
  std::cout << std::endl;

  // Transform the cumulative counts into the C array and mark bits i + C[i].
  // Add a guard to the end to support determining C[sigma].
  size_type cumulative = 0;
  counts.push_back(0);
  for(size_type i = 0; i < counts.size(); i++)
  {
    size_type temp = counts[i]; counts[i] = cumulative + i; cumulative += temp;
  }
  sdsl::sd_vector<> C(counts.begin(), counts.end());
  sdsl::store_to_file(C, alphabet_name);

  double seconds = readTimer() - start;
  size_type data_size = header.sequences + header.total_length;

  std::cout << "Processed " << data_size << " nodes in " << seconds << " seconds (" << (data_size / seconds) << " nodes/second)" << std::endl;
  std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  std::cerr << "Usage: prepare_text [options] base_name" << std::endl;
  std::cerr << "  -m N  Read up to N sequences" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Build header, sequence borders, and alphabet from the sequences." << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

//------------------------------------------------------------------------------
