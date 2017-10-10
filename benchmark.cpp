/*
  Copyright (c) 2017 Jouni Siren

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

const size_type RANDOM_SEED  = 0xDEADBEEF;
const size_type QUERIES      = 20000;
const size_type QUERY_LENGTH = 60;

void printUsage(int exit_code = EXIT_SUCCESS);

std::vector<std::vector<node_type>> generateQueries(const std::string& base_name);

template<class GBWTType>
std::vector<SearchState> findBenchmark(const GBWTType& index, const std::vector<std::vector<node_type>>& queries);

size_type totalLength(const std::vector<SearchState>& states);

template<class GBWTType>
void locateBenchmark(const GBWTType& index, const std::vector<SearchState>& queries);

void extractBenchmark(const GBWT& compressed_index, const DynamicGBWT& dynamic_index);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 3) { printUsage(); }
  std::string index_base = argv[1], query_base = argv[2];

  std::cout << "GBWT benchmark" << std::endl;
  std::cout << std::endl;

  printHeader("Index name"); std::cout << index_base << std::endl;
  printHeader("Query name"); std::cout << query_base << std::endl;
  std::cout << std::endl;

  double start = readTimer();

  GBWT compressed_index;
  sdsl::load_from_file(compressed_index, index_base + GBWT::EXTENSION);
  printStatistics(compressed_index, index_base);

  DynamicGBWT dynamic_index;
  sdsl::load_from_file(dynamic_index, index_base + DynamicGBWT::EXTENSION);
  printStatistics(dynamic_index, index_base);

  std::cout << "find() benchmarks:" << std::endl;
  std::vector<std::vector<node_type>> queries = generateQueries(query_base);
  std::vector<SearchState> results = findBenchmark(compressed_index, queries);
  findBenchmark(dynamic_index, queries);
  std::cout << results.size() << " ranges of total length " << totalLength(results) << std::endl;
  std::cout << std::endl;

  locateBenchmark(compressed_index, results);
  locateBenchmark(dynamic_index, results);

  extractBenchmark(compressed_index, dynamic_index);

  double seconds = readTimer() - start;
  std::cout << "Benchmarks completed in " << seconds << " seconds" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  std::cerr << "Usage: benchmark index_base query_base" << std::endl;
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

std::string indexType(const GBWT&) { return "Compressed GBWT"; }
std::string indexType(const DynamicGBWT&) { return "Dynamic GBWT"; }

size_type
totalLength(const std::vector<SearchState>& states)
{
  size_type result = 0;
  for(SearchState state : states) { result += state.size(); }
  return result;
}

//------------------------------------------------------------------------------

template<class GBWTType>
std::vector<SearchState>
findBenchmark(const GBWTType& index, const std::vector<std::vector<node_type>>& queries)
{
  std::vector<SearchState> results(queries.size());

  double start = readTimer();
  for(size_type i = 0; i < queries.size(); i++)
  {
    results[i] = index.find(queries[i].begin(), queries[i].end());
  }
  double seconds = readTimer() - start;

  printTime(indexType(index), queries.size(), seconds);
  return results;
}

//------------------------------------------------------------------------------

template<class GBWTType>
void
locateBenchmark(const GBWTType& index, const std::vector<SearchState>& queries)
{
  std::cout << "locate() benchmarks (" << indexType(index) << "):" << std::endl;

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
      found += result.size();
    }
    double seconds = readTimer() - start;
    printTime("Direct", found, seconds);
  }

  {
    double start = readTimer();
    size_type found = 0;
    for(SearchState query : queries)
    {
      std::vector<size_type> result = index.locate(query);
      found += result.size();
    }
    double seconds = readTimer() - start;
    printTime("Fast", found, seconds);
  }

  std::cout << std::endl;
}

//------------------------------------------------------------------------------

template<class GBWTType>
void
extractBenchmark(const GBWTType& index)
{
  double start = readTimer();
  size_type total_length = 0;
  for(size_type i = 0; i < index.sequences(); i++)
  {
    std::vector<node_type> sequence = index.extract(i);
    total_length += sequence.size() + 1;
  }
  double seconds = readTimer() - start;
  printTime(indexType(index), index.sequences(), seconds);
  if(total_length != index.size())
  {
    std::cerr << "extractBenchmark(): " << indexType(index) << ": Total length " << total_length << ", expected " << index.size() << std::endl;
  }
}

void
extractBenchmark(const GBWT& compressed_index, const DynamicGBWT& dynamic_index)
{
  std::cout << "extract() benchmarks:" << std::endl;
  extractBenchmark(compressed_index);
  extractBenchmark(dynamic_index);
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
