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

double build(DynamicGBWT& gbwt, const std::string& base_name, size_type batch_size);

void verify(DynamicGBWT& gbwt, const std::string& base_name);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2) { printUsage(); }

  size_type batch_size = 0;
  bool verify_index = false;
  int c = 0;
  while((c = getopt(argc, argv, "b:v")) != -1)
  {
    switch(c)
    {
    case 'b':
      batch_size = std::stoul(optarg); break;
    case 'v':
      verify_index = true; break;
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

  if(verify_index) { verify(gbwt, base_name); }

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  std::cerr << "Usage: build_gbwt [options] base_name" << std::endl;
  std::cerr << "  -b N  Insert in batches of N million nodes" << std::endl;
  std::cerr << "  -v    Verify the index after construction" << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

//------------------------------------------------------------------------------

double
build(DynamicGBWT& gbwt, const std::string& base_name, size_type batch_size)
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

void
verify(DynamicGBWT& gbwt, const std::string& base_name)
{
  double start = readTimer();

  // Read the input and find the starting offsets.
  std::vector<size_type> offsets;
  {
    text_buffer_type text(base_name);
    bool seq_start = true;
    for(size_type i = 0; i < text.size(); i++)
    {
      if(seq_start) { offsets.push_back(i); seq_start = false; }
      if(text[i] == ENDMARKER) { seq_start = true; }
    }
  }
  if(offsets.empty()) { return; }
  std::vector<range_type> blocks = Range::partition(range_type(0, offsets.size() - 1), 4 * omp_get_max_threads());

  bool failed = false;
  #pragma omp parallel for schedule(dynamic, 1)
  for(size_type block = 0; block < blocks.size(); block++)
  {
    text_buffer_type text(base_name);
    for(size_type sequence = blocks[block].first; sequence <= blocks[block].second; sequence++)
    {
      edge_type current(ENDMARKER, sequence);
      for(size_type offset = offsets[sequence]; text[offset] != ENDMARKER; offset++)
      {
        edge_type next = gbwt.LF(current.first, current.second);
        if(next.first != text[offset])
        {
          #pragma omp critical
          {
            std::cerr << "build_gbwt: Index verification failed with sequence " << sequence << ", offset "
                      << (offset - offsets[sequence]) << std::endl;
            std::cerr << "build_gbwt: Expected an edge from " << current.first << " to " << text[offset]
                      << ", ended up in " << next.first << std::endl;
            failed = true;
          }
          break;
        }
        current = next;
      }
    }
  }

  double seconds = readTimer() - start;

  if(failed) { std::cout << "Index verification failed" << std::endl; }
  else { std::cout << "Index verified in " << seconds << " seconds" << std::endl; }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
