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

#include <random>
#include <unistd.h>

#include <gbwt/dynamic_gbwt.h>

using namespace gbwt;

//------------------------------------------------------------------------------

void printUsage(int exit_code = EXIT_SUCCESS);

template<class GBWTType>
void verify(const std::string& base_name);

void verifyLocate(const std::string& base_name);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2) { printUsage(); }

  size_type batch_size = DynamicGBWT::INSERT_BATCH_SIZE / MILLION;
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

  double start = readTimer();

  DynamicGBWT gbwt;
  text_buffer_type input(base_name);
  gbwt.insert(input, batch_size * MILLION);

  std::string gbwt_name = base_name + DynamicGBWT::EXTENSION;
  sdsl::store_to_file(gbwt, gbwt_name);

  double seconds = readTimer() - start;

  printStatistics(gbwt, base_name);

  std::cout << "Indexed " << gbwt.size() << " nodes in " << seconds << " seconds (" << (gbwt.size() / seconds) << " nodes/second)" << std::endl;
  std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  sdsl::util::clear(gbwt);
  if(verify_index)
  {
    std::cout << "Verifying compressed GBWT..." << std::endl;
    verify<GBWT>(base_name);

    std::cout << "Verifying dynamic GBWT..." << std::endl;
    verify<DynamicGBWT>(base_name);

    std::cout << "Verifying locate()..." << std::endl;
    verifyLocate(base_name);
  }

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  std::cerr << "Usage: build_gbwt [options] base_name" << std::endl;
  std::cerr << "  -b N  Insert in batches of N million nodes (default "
            << (DynamicGBWT::INSERT_BATCH_SIZE / MILLION) << ")" << std::endl;
  std::cerr << "  -v    Verify the index after construction" << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

//------------------------------------------------------------------------------

template<class GBWTType>
void
verify(const std::string& base_name)
{
  double start = readTimer();

  GBWTType gbwt;
  sdsl::load_from_file(gbwt, base_name + GBWTType::EXTENSION);

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
  std::atomic<size_type> samples_found(0);
  #pragma omp parallel for schedule(dynamic, 1)
  for(size_type block = 0; block < blocks.size(); block++)
  {
    text_buffer_type text(base_name);
    for(size_type sequence = blocks[block].first; sequence <= blocks[block].second; sequence++)
    {
      edge_type current(ENDMARKER, sequence);
      size_type offset = offsets[sequence];
      while(true)
      {
        // Check for a sample.
        size_type sample = gbwt.tryLocate(current);
        if(sample != invalid_sequence())
        {
          samples_found++;
          if(sample != sequence)
          {
            #pragma omp critical
            {
              std::cerr << "build_gbwt: Index verification failed with sequence " << sequence << ", offset "
                        << (offset - offsets[sequence]) << std::endl;
              std::cerr << "build_gbwt: Sample had sequence id " << sample << std::endl;
              failed = true;
            }
            break;
          }
        }

        // Verify LF().
        if(text[offset] == ENDMARKER) { break; }
        edge_type next = gbwt.LF(current);
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
        current = next; offset++;
      }
    }
  }

  if(samples_found != gbwt.samples())
  {
    std::cerr << "build_gbwt: Found " << samples_found << " samples, expected " << gbwt.samples() << std::endl;
    failed = true;
  }

  double seconds = readTimer() - start;

  if(failed) { std::cout << "Index verification failed" << std::endl; }
  else { std::cout << "Index verified in " << seconds << " seconds" << std::endl; }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------

const size_type RANDOM_SEED    = 0xDEADBEEF;
const size_type LOCATE_QUERIES = 10000;
const size_type MAX_LENGTH     = 100;

void
verifyLocate(const std::string& base_name)
{
  GBWT gbwt;
  sdsl::load_from_file(gbwt, base_name + GBWT::EXTENSION);

  std::mt19937_64 rng(RANDOM_SEED);
  std::vector<node_type>  nodes(LOCATE_QUERIES);
  std::vector<range_type> ranges(LOCATE_QUERIES);
  size_type total_length = 0;
  for(size_type i = 0; i < LOCATE_QUERIES; i++)
  {
    size_type node_size = 0;
    while(node_size == 0)
    {
      nodes[i] = rng() % gbwt.effective();
      if(nodes[i] > 0) { nodes[i] += gbwt.header.offset; }
      node_size = gbwt.count(nodes[i]);
    }
    ranges[i].first = rng() % node_size;
    ranges[i].second = ranges[i].first + rng() % std::min(MAX_LENGTH, node_size - ranges[i].first);
    total_length += Range::length(ranges[i]);
  }
  std::cout << "Generated " << LOCATE_QUERIES << " ranges of total length " << total_length << std::endl;

  size_type direct_hash = FNV_OFFSET_BASIS;
  {
    double start = readTimer();
    size_type found = 0;
    for(size_type i = 0; i < LOCATE_QUERIES; i++)
    {
      std::vector<size_type> result;
      for(size_type j = ranges[i].first; j <= ranges[i].second; j++)
      {
        result.push_back(gbwt.locate(nodes[i], j));
      }
      removeDuplicates(result, false);
      for(size_type res : result) { direct_hash = fnv1a_hash(res, direct_hash); }
      found += result.size();
    }
    double seconds = readTimer() - start;
    printTime("Direct locate()", found, seconds);
  }

  size_type fast_hash = FNV_OFFSET_BASIS;
  {
    double start = readTimer();
    size_type found = 0;
    for(size_type i = 0; i < LOCATE_QUERIES; i++)
    {
      std::vector<size_type> result = gbwt.locate(nodes[i], ranges[i]);
      for(size_type res : result) { fast_hash = fnv1a_hash(res, fast_hash); }
      found += result.size();
    }
    double seconds = readTimer() - start;
    printTime("Fast locate()", found, seconds);
  }

  if(direct_hash != fast_hash) { std::cout << "Verification failed" << std::endl; }
  else { std::cout << "Verification successful" << std::endl; }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
