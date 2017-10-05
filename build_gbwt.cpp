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

const size_type MAX_ERRORS   = 100; // Do not print more error messages.
size_type errors             = 0;

const size_type RANDOM_SEED  = 0xDEADBEEF;
const size_type QUERIES      = 20000;
const size_type QUERY_LENGTH = 60;

void printUsage(int exit_code = EXIT_SUCCESS);

std::vector<size_type> startOffsets(const std::string& base_name);
std::vector<std::vector<node_type>> generateQueries(const std::string& base_name);

template<class GBWTType>
void verifyExtract(const GBWTType& index, const std::string& base_name, const std::vector<size_type>& offsets);

template<class GBWTType>
void verifySamples(const GBWTType& index);

template<class GBWTType>
void verifyLocate(const GBWTType& index, const std::vector<SearchState>& queries);

std::vector<SearchState> verifyFind(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const std::vector<std::vector<node_type>>& queries);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2) { printUsage(); }

  size_type batch_size = DynamicGBWT::INSERT_BATCH_SIZE / MILLION;
  bool verify_index = false;
  std::string index_base, input_base, output_base;
  int c = 0;
  while((c = getopt(argc, argv, "b:i:o:v")) != -1)
  {
    switch(c)
    {
    case 'b':
      batch_size = std::stoul(optarg); break;
    case 'i':
      index_base = optarg; break;
    case 'o':
      output_base = optarg; break;
    case 'v':
      verify_index = true; break;
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }

  size_type input_files = argc - optind;
  size_type input_size = 0;
  if(index_base.empty() && output_base.empty() && input_files == 1) { output_base = argv[optind]; }
  if(input_files == 0 || output_base.empty()) { printUsage(EXIT_FAILURE); }
  if(verify_index && !(input_files == 1 && index_base.empty()))
  {
    std::cerr << "build_gbwt: Verification only works with indexes for a single file" << std::endl;
    verify_index = false;
  }

  std::cout << "GBWT construction" << std::endl;
  std::cout << std::endl;

  if(!(index_base.empty())) { printHeader("Index name"); std::cout << index_base << std::endl; }
  printHeader("Input files"); std::cout << input_files << std::endl;
  printHeader("Output name"); std::cout << output_base << std::endl;
  if(batch_size != 0) { printHeader("Batch size"); std::cout << batch_size << " million" << std::endl; }
  std::cout << std::endl;

  double start = readTimer();

  DynamicGBWT dynamic_index;
  if(!(index_base.empty()))
  {
    sdsl::load_from_file(dynamic_index, index_base + DynamicGBWT::EXTENSION);
    printStatistics(dynamic_index, index_base);
  }

  while(optind < argc)
  {
    input_base = argv[optind];
    printHeader("Input name"); std::cout << input_base << std::endl;
    text_buffer_type input(input_base);
    input_size += input.size();
    dynamic_index.insert(input, batch_size * MILLION);
    optind++;
  }
  std::cout << std::endl;

  std::string gbwt_name = output_base + DynamicGBWT::EXTENSION;
  sdsl::store_to_file(dynamic_index, gbwt_name);
  printStatistics(dynamic_index, output_base);

  double seconds = readTimer() - start;

  std::cout << "Indexed " << input_size << " nodes in " << seconds << " seconds (" << (input_size / seconds) << " nodes/second)" << std::endl;
  std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  if(verify_index)
  {
    std::cout << "Verifying the index..." << std::endl;
    std::vector<size_type> offsets = startOffsets(input_base);
    std::vector<std::vector<node_type>> queries = generateQueries(input_base);
    std::cout << std::endl;

    GBWT compressed_index;
    sdsl::load_from_file(compressed_index, gbwt_name);

    sdsl::util::clear(dynamic_index);
    sdsl::load_from_file(dynamic_index, gbwt_name);

    std::cout << "Verifying find()..." << std::endl;
    std::vector<SearchState> results = verifyFind(compressed_index, dynamic_index, queries);

    std::cout << "Verifying locate() in compressed GBWT..." << std::endl;
    verifyLocate(compressed_index, results);

    std::cout << "Verifying locate() in dynamic GBWT..." << std::endl;
    verifyLocate(dynamic_index, results);

    std::cout << "Verifying samples in compressed GBWT..." << std::endl;
    verifySamples(compressed_index);

    std::cout << "Verifying samples in dynamic GBWT..." << std::endl;
    verifySamples(dynamic_index);

    std::cout << "Verifying extract() in compressed GBWT..." << std::endl;
    verifyExtract(compressed_index, input_base, offsets);

    std::cout << "Verifying extract() in dynamic GBWT..." << std::endl;
    verifyExtract(dynamic_index, input_base, offsets);

    if(errors > 0) { std::cout << "Index verification failed" << std::endl; }
    else { std::cout << "Index verification successful" << std::endl; }
    std::cout << std::endl;
  }

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  std::cerr << "Usage: build_gbwt [options] input1 [input2 ...]" << std::endl;
  std::cerr << "  -b N  Insert in batches of N million nodes (default: " << (DynamicGBWT::INSERT_BATCH_SIZE / MILLION) << ")" << std::endl;
  std::cerr << "  -i X  Insert the sequences into an existing index with base name X" << std::endl;
  std::cerr << "  -o X  Use base name X for output (default: the only input)" << std::endl;
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

std::vector<std::vector<node_type>>
generateQueries(const std::string& base_name)
{
  std::mt19937_64 rng(RANDOM_SEED);
  std::vector<std::vector<node_type>> result;
  text_buffer_type text(base_name);
  if(text.size() <= QUERY_LENGTH)
  {
    std::cerr << "generateQueries(): Text length " << text.size() << " is too small for query length " << QUERY_LENGTH << std::endl;
    return result;
  }

  size_type attempts = 0;
  while(result.size() < QUERIES && attempts < 2 * QUERIES)
  {
    size_type start_offset = rng() % (text.size() - QUERY_LENGTH);  // We do not want queries containing the final endmarker.
    std::vector<node_type> candidate(text.begin() + start_offset, text.begin() + start_offset + QUERY_LENGTH);
    bool ok = true;
    for(node_type node : candidate)
    {
      if(node == ENDMARKER) { ok = false; break; }
    }
    if(ok) { result.push_back(candidate); }
    attempts++;
  }

  std::cout << "Generated " << result.size() << " queries of total length " << (result.size() * QUERY_LENGTH) << std::endl;
  return result;
}

//------------------------------------------------------------------------------

template<class GBWTType>
void
verifyExtract(const GBWTType& index, const std::string& base_name, const std::vector<size_type>& offsets)
{
  double start = readTimer();

  size_type initial_errors = errors;
  if(index.sequences() != offsets.size())
  {
    errors++;
    if(errors <= MAX_ERRORS)
    {
      std::cerr << "verifyExtract(): Expected " << offsets.size() << " sequences, got " << index.sequences() << std::endl;
    }
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
            errors++;
            if(errors <= MAX_ERRORS)
            {
              std::cerr << "verifyExtract(): Verification failed with sequence " << sequence << ", offset " << i << std::endl;
              std::cerr << "verifyExtract(): Expected " << text[offsets[sequence] + i] << ", got " << result[i] << std::endl;
            }
          }
          break;
        }
      }
      if(text[offsets[sequence] + result.size()] != ENDMARKER)
      {
        #pragma omp critical
        {
          errors++;
          if(errors <= MAX_ERRORS)
          {
            std::cerr << "verifyExtract(): Verification failed with sequence " << sequence << ", offset " << result.size() << std::endl;
            std::cerr << "verifyExtract(): Expected " << text[offsets[sequence] + result.size()] << ", got endmarker" << std::endl;
          }
        }
      }
    }
  }

  double seconds = readTimer() - start;

  if(errors > initial_errors) { std::cout << "extract() verification failed" << std::endl; }
  else { std::cout << "extract() verified in " << seconds << " seconds" << std::endl; }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------

template<class GBWTType>
void
verifySamples(const GBWTType& index)
{
  double start = readTimer();

  size_type initial_errors = errors;
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
              errors++;
              if(errors <= MAX_ERRORS)
              {
                std::cerr << "verifySamples(): Verification failed with sequence " << sequence << ", position " << current << std::endl;
                std::cerr << "verifySamples(): Sample had sequence id " << sample << std::endl;
              }
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
    errors++;
    if(errors <= MAX_ERRORS)
    {
      std::cerr << "verifySamples(): Found " << samples_found << " samples, expected " << index.samples() << std::endl;
    }
  }

  double seconds = readTimer() - start;

  if(errors > initial_errors) { std::cout << "Sample verification failed" << std::endl; }
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
    for(SearchState query : queries)
    {
      if(query.empty()) { continue; }
      std::vector<size_type> result;
      for(size_type i = query.range.first; i <= query.range.second; i++)
      {
        result.push_back(index.locate(query.node, i));
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
    for(SearchState query : queries)
    {
      std::vector<size_type> result = index.locate(query);
      for(size_type res : result) { fast_hash = fnv1a_hash(res, fast_hash); }
      found += result.size();
    }
    double seconds = readTimer() - start;
    printTime("Fast locate()", found, seconds);
  }

  if(direct_hash != fast_hash) { errors++; std::cout << "locate() verification failed" << std::endl; }
  else { std::cout << "locate() verification successful" << std::endl; }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------

std::vector<SearchState>
verifyFind(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const std::vector<std::vector<node_type>>& queries)
{
  std::vector<SearchState> results(queries.size());
  size_type total_length = 0;
  size_type initial_errors = errors;

  {
    double start = readTimer();
    for(size_type i = 0; i < queries.size(); i++)
    {
      results[i] = compressed_index.find(queries[i].begin(), queries[i].end());
      total_length += results[i].size();
    }
    double seconds = readTimer() - start;
    printTime("Compressed GBWT", queries.size(), seconds);
  }

  {
    double start = readTimer();
    for(size_type i = 0; i < queries.size(); i++)
    {
      SearchState result = dynamic_index.find(queries[i].begin(), queries[i].end());
      if(result != results[i])
      {
        errors++;
        if(errors <= MAX_ERRORS)
        {
          std::cerr << "verifyFind(): Mismatching results with query " << i << std::endl;
          std::cerr << "verifyFind(): Compressed: " << results[i] << ", dynamic: " << result << std::endl;
        }
      }
    }
    double seconds = readTimer() - start;
    printTime("Dynamic GBWT", queries.size(), seconds);
  }

  std::cout << "Generated " << results.size() << " ranges of total length " << total_length << std::endl;
  if(errors > initial_errors) { std::cout << "find() verification failed" << std::endl; }
  else { std::cout << "find() verification successful" << std::endl; }
  std::cout << std::endl;

  return results;
}

//------------------------------------------------------------------------------
