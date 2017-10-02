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

#include <random>
#include <unistd.h>

#include <gbwt/dynamic_gbwt.h>

using namespace gbwt;

//------------------------------------------------------------------------------

void printUsage(int exit_code = EXIT_SUCCESS);

std::vector<size_type> startOffsets(const std::string& base_name);
std::vector<SearchState> queryRanges(const DynamicGBWT& index);

template<class GBWTType>
void verifyExtract(const GBWTType& index, const std::string& base_name, const std::vector<size_type>& offsets);

template<class GBWTType>
void verifySamples(const GBWTType& index);

template<class GBWTType>
void verifyLocate(const GBWTType& index, const std::vector<SearchState>& queries);

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

  DynamicGBWT dynamic_index;
  text_buffer_type input(base_name);
  dynamic_index.insert(input, batch_size * MILLION);

  std::string gbwt_name = base_name + DynamicGBWT::EXTENSION;
  sdsl::store_to_file(dynamic_index, gbwt_name);

  double seconds = readTimer() - start;

  printStatistics(dynamic_index, base_name);

  std::cout << "Indexed " << dynamic_index.size() << " nodes in " << seconds << " seconds (" << (dynamic_index.size() / seconds) << " nodes/second)" << std::endl;
  std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  if(verify_index)
  {
    std::cout << "Verifying the index..." << std::endl;
    std::vector<size_type> offsets = startOffsets(base_name);
    std::vector<SearchState> queries = queryRanges(dynamic_index);
    std::cout << std::endl;

    GBWT compressed_index;
    sdsl::load_from_file(compressed_index, gbwt_name);

    sdsl::util::clear(dynamic_index);
    sdsl::load_from_file(dynamic_index, gbwt_name);

    std::cout << "Verifying extract() in compressed GBWT..." << std::endl;
    verifyExtract(compressed_index, base_name, offsets);

    std::cout << "Verifying extract() in dynamic GBWT..." << std::endl;
    verifyExtract(dynamic_index, base_name, offsets);

    std::cout << "Verifying samples in compressed GBWT..." << std::endl;
    verifySamples(compressed_index);

    std::cout << "Verifying samples in dynamic GBWT..." << std::endl;
    verifySamples(dynamic_index);

    std::cout << "Verifying locate() in compressed GBWT..." << std::endl;
    verifyLocate(compressed_index, queries);

    std::cout << "Verifying locate() in dynamic GBWT..." << std::endl;
    verifyLocate(dynamic_index, queries);
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

std::vector<size_type>
startOffsets(const std::string& base_name)
{
  std::vector<size_type> offsets;
  text_buffer_type text(base_name);
  bool seq_start = true;
  for(size_type i = 0; i < text.size(); i++)
  {
    if(seq_start) { offsets.push_back(i); seq_start = false; }
    if(text[i] == ENDMARKER) { seq_start = true; }
  }

  std::cout << "Found the starting offsets for " << offsets.size() << " sequences" << std::endl;
  return offsets;
}

const size_type RANDOM_SEED    = 0xDEADBEEF;
const size_type LOCATE_QUERIES = 10000;
const size_type MAX_LENGTH     = 100;

std::vector<SearchState>
queryRanges(const DynamicGBWT& index)
{
  std::mt19937_64 rng(RANDOM_SEED);
  std::vector<SearchState> result(LOCATE_QUERIES);
  size_type total_length = 0;
  for(size_type i = 0; i < LOCATE_QUERIES; i++)
  {
    size_type node_size = 0;
    while(node_size == 0)
    {
      result[i].node = rng() % index.effective();
      if(result[i].node > 0) { result[i].node += index.header.offset; }
      node_size = index.nodeSize(result[i].node);
    }
    result[i].range.first = rng() % node_size;
    result[i].range.second = result[i].range.first + rng() % std::min(MAX_LENGTH, node_size - result[i].range.first);
    total_length += Range::length(result[i].range);
  }

  std::cout << "Generated " << LOCATE_QUERIES << " ranges of total length " << total_length << std::endl;
  return result;
}

//------------------------------------------------------------------------------

template<class GBWTType>
void
verifyExtract(const GBWTType& index, const std::string& base_name, const std::vector<size_type>& offsets)
{
  double start = readTimer();

  bool failed = false;
  if(index.sequences() != offsets.size())
  {
    std::cerr << "verifyExtract(): Expected " << offsets.size() << " sequences, got " << index.sequences() << std::endl;
    failed = true;
  }
  std::vector<range_type> blocks = Range::partition(range_type(0, offsets.size() - 1), 4 * omp_get_max_threads());

  #pragma omp parallel for schedule(dynamic, 1)
  for(size_type block = 0; block < blocks.size(); block++)
  {
    text_buffer_type text(base_name);
    for(size_type sequence = blocks[block].first; sequence <= blocks[block].second; sequence++)
    {
      std::vector<node_type> result = index.extract(sequence);
      for(size_type i = 0; i < result.size(); i++)
      {
        if(result[i] != text[offsets[sequence] + i])
        {
          #pragma omp critical
          {
            std::cerr << "verifyExtract(): Verification failed with sequence " << sequence << ", offset " << i << std::endl;
            std::cerr << "verifyExtract(): Expected " << text[offsets[sequence] + i] << ", got " << result[i] << std::endl;
            failed = true;
          }
          break;
        }
      }
      if(text[offsets[sequence] + result.size()] != ENDMARKER)
      {
        #pragma omp critical
        {
          std::cerr << "verifyExtract(): Verification failed with sequence " << sequence << ", offset " << result.size() << std::endl;
          std::cerr << "verifyExtract(): Expected " << text[offsets[sequence] + result.size()] << ", got endmarker" << std::endl;
          failed = true;
        }
      }
    }
  }

  double seconds = readTimer() - start;

  if(failed) { std::cout << "extract() verification failed" << std::endl; }
  else { std::cout << "extract() verified in " << seconds << " seconds" << std::endl; }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------

template<class GBWTType>
void
verifySamples(const GBWTType& index)
{
  double start = readTimer();

  bool failed = false;
  std::atomic<size_type> samples_found(0);
  std::vector<range_type> blocks = Range::partition(range_type(0, index.sequences() - 1), 4 * omp_get_max_threads());

  #pragma omp parallel for schedule(dynamic, 1)
  for(size_type block = 0; block < blocks.size(); block++)
  {
    for(size_type sequence = blocks[block].first; sequence <= blocks[block].second; sequence++)
    {
      edge_type current(ENDMARKER, sequence);
      do
      {
        size_type sample = index.tryLocate(current);
        if(sample != invalid_sequence())
        {
          samples_found++;
          if(sample != sequence)
          {
            #pragma omp critical
            {
              std::cerr << "verifySamples(): Verification failed with sequence " << sequence << ", position " << current << std::endl;
              std::cerr << "verifySamples(): Sample had sequence id " << sample << std::endl;
              failed = true;
            }
            break;
          }
        }
        current = index.LF(current);
      }
      while(current.first != ENDMARKER);
    }
  }

  if(samples_found != index.samples())
  {
    std::cerr << "verifySamples(): Found " << samples_found << " samples, expected " << index.samples() << std::endl;
    failed = true;
  }

  double seconds = readTimer() - start;

  if(failed) { std::cout << "Sample verification failed" << std::endl; }
  else { std::cout << "Samples verified in " << seconds << " seconds" << std::endl; }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------

template<class GBWTType>
void
verifyLocate(const GBWTType& index, const std::vector<SearchState>& queries)
{
  size_type direct_hash = FNV_OFFSET_BASIS;
  {
    double start = readTimer();
    size_type found = 0;
    for(size_type i = 0; i < LOCATE_QUERIES; i++)
    {
      std::vector<size_type> result;
      for(size_type j = queries[i].range.first; j <= queries[i].range.second; j++)
      {
        result.push_back(index.locate(queries[i].node, j));
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
      std::vector<size_type> result = index.locate(queries[i]);
      for(size_type res : result) { fast_hash = fnv1a_hash(res, fast_hash); }
      found += result.size();
    }
    double seconds = readTimer() - start;
    printTime("Fast locate()", found, seconds);
  }

  if(direct_hash != fast_hash) { std::cout << "locate() verification failed" << std::endl; }
  else { std::cout << "locate() verification successful" << std::endl; }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
