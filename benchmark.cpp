/*
  Copyright (c) 2017, 2018 Jouni Siren

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

const std::string tool_name = "GBWT benchmark";

const size_type RANDOM_SEED  = 0xDEADBEEF;

void printUsage(int exit_code = EXIT_SUCCESS);

std::vector<SearchState> findBenchmark(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, size_type find_queries, size_type pattern_length, std::vector<std::vector<node_type>>& queries);

void bidirectionalBenchmark(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const std::vector<std::vector<node_type>>& queries);

template<class GBWTType>
void locateBenchmark(const GBWTType& index, const std::vector<SearchState>& queries);

void extractBenchmark(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, size_type extract_queries);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2) { printUsage(); }

  int c = 0;
  bool find = false, locate = false, extract = false;
  size_type find_queries = 0, pattern_length = 0, extract_queries = 0;
  while((c = getopt(argc, argv, "f:p:le:")) != -1)
  {
    switch(c)
    {
    case 'f':
      find = true;
      find_queries = std::stoul(optarg); break;
    case 'p':
      pattern_length = std::stoul(optarg); break;
    case 'l':
      locate = true; break;
    case 'e':
      extract = true;
      extract_queries = std::stoul(optarg); break;
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }

  if(optind >= argc) { printUsage(EXIT_FAILURE); }
  std::string index_base = argv[optind];
  if(find)
  {
    if(find_queries == 0 || pattern_length == 0)
    {
      std::cerr << "benchmark: Number of queries and pattern length must be non-zero" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  if(locate)
  {
    if(!find)
    {
      std::cerr << "benchmark: Cannot benchmark locate() queries without find() queries" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  if(extract)
  {
    if(extract_queries == 0)
    {
      std::cerr << "benchmark: Number of queries must be non-zero" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  Version::print(std::cout, tool_name);
  printHeader("Index name"); std::cout << index_base << std::endl;
  if(find || locate || extract)
  {
    printHeader("Queries");
    if(find) { std::cout << "find(" << find_queries << ", " << pattern_length << ") "; }
    if(locate) { std::cout << "locate() "; }
    if(extract) { std::cout << "extract(" << extract_queries << ") "; }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  GBWT compressed_index;
  sdsl::load_from_file(compressed_index, index_base + GBWT::EXTENSION);
  printStatistics(compressed_index, index_base);
  if(!(find || locate || extract)) { return 0; }

  DynamicGBWT dynamic_index;
  sdsl::load_from_file(dynamic_index, index_base + DynamicGBWT::EXTENSION);
  printStatistics(dynamic_index, index_base);

  if(find)
  {
    std::vector<std::vector<node_type>> queries;
    std::vector<SearchState> results = findBenchmark(compressed_index, dynamic_index, find_queries, pattern_length, queries);
    if(compressed_index.bidirectional())
    {
      bidirectionalBenchmark(compressed_index, dynamic_index, queries);
    }
    queries.clear();
    if(locate)
    {
      locateBenchmark(compressed_index, results);
      locateBenchmark(dynamic_index, results);
    }
  }
  if(extract)
  {
    extractBenchmark(compressed_index, dynamic_index, extract_queries);    
  }

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: benchmark [options] index_base" << std::endl;
  std::cerr << "  -f N  Benchmark N find() queries (requires -p)" << std::endl;
  std::cerr << "  -p N  Use patterns of length N" << std::endl;
  std::cerr << "  -l    Benchmark locate() queries (requires -f)" << std::endl;
  std::cerr << "  -e N  Benchmark N extract() queries" << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

std::vector<std::vector<node_type>>
generateQueries(const DynamicGBWT& index, size_type find_queries, size_type pattern_length)
{
  std::mt19937_64 rng(RANDOM_SEED);
  std::vector<std::vector<node_type>> result;

  size_type attempts = 0;
  while(result.size() < find_queries && attempts < 2 * find_queries)
  {
    std::vector<node_type> sequence = index.extract(rng() % index.sequences());
    if(sequence.size() >= pattern_length)
    {
      size_type start_offset = rng() % (sequence.size() + 1 - pattern_length);
      result.push_back(std::vector<node_type>(sequence.begin() + start_offset, sequence.begin() + start_offset + pattern_length));
    }
    attempts++;
  }

  std::cout << "Generated " << result.size() << " queries of total length " << (result.size() * pattern_length) << std::endl;
  return result;
}

std::string indexType(const GBWT&) { return "Compressed GBWT"; }
std::string indexType(const DynamicGBWT&) { return "Dynamic GBWT"; }

//------------------------------------------------------------------------------

template<class GBWTType>
size_type
findBenchmark(const GBWTType& index, const std::vector<std::vector<node_type>>& queries, std::vector<SearchState>& results)
{
  double start = readTimer();
  size_type total_length = 0, total_size = 0;
  for(size_type i = 0; i < queries.size(); i++)
  {
    results[i] = index.find(queries[i].begin(), queries[i].end());
    total_length += queries[i].size();
    total_size += results[i].size();
  }
  double seconds = readTimer() - start;
  printTimeLength(indexType(index), queries.size(), total_length, seconds);
  return total_size;
}    

std::vector<SearchState>
findBenchmark(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, size_type find_queries, size_type pattern_length, std::vector<std::vector<node_type>>& queries)
{
  std::cout << "find() benchmarks:" << std::endl;

  queries = generateQueries(dynamic_index, find_queries, pattern_length);
  std::vector<SearchState> results(queries.size());

  size_type compressed_length = findBenchmark(compressed_index, queries, results);
  size_type dynamic_length = findBenchmark(dynamic_index, queries, results);

  if(compressed_length != dynamic_length)
  {
    std::cerr << "findBenchmark(): Total length mismatch: "
              << compressed_length << " (" << indexType(compressed_index) << "), "
              << dynamic_length << " (" << indexType(dynamic_index) << ")" << std::endl;
  }

  std::cout << "Found " << results.size() << " ranges of total length " << compressed_length << std::endl;
  std::cout << std::endl;

  return results;
}

//------------------------------------------------------------------------------

template<class GBWTType>
BidirectionalState
bidirectionalSearch(const GBWTType& index, const std::vector<node_type>& query)
{
  if(query.empty()) { return BidirectionalState(); }

  size_type midpoint = query.size() / 2;
  BidirectionalState state = index.bdFind(query[midpoint]);
  for(size_type i = midpoint + 1; !(state.empty()) && i < query.size(); i++)
  {
    state = index.bdExtendForward(state, query[i]);
  }
  for(size_type i = midpoint; !(state.empty()) && i > 0; i--)
  {
    state = index.bdExtendBackward(state, query[i - 1]);
  }

  return state;
}

template<class GBWTType>
size_type
bidirectionalBenchmark(const GBWTType& index, const std::vector<std::vector<node_type>>& queries)
{
  double start = readTimer();
  size_type total_length = 0, total_size = 0;
  for(size_type i = 0; i < queries.size(); i++)
  {
    BidirectionalState result = bidirectionalSearch(index, queries[i]);
    total_length += queries[i].size();
    total_size += result.size();
  }
  double seconds = readTimer() - start;
  printTimeLength(indexType(index), queries.size(), total_length, seconds);
  return total_size;
}    

void
bidirectionalBenchmark(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const std::vector<std::vector<node_type>>& queries)
{
  std::cout << "Bidirectional benchmarks:" << std::endl;

  size_type compressed_length = bidirectionalBenchmark(compressed_index, queries);
  size_type dynamic_length = bidirectionalBenchmark(dynamic_index, queries);

  if(compressed_length != dynamic_length)
  {
    std::cerr << "bidirectionalBenchmark(): Total length mismatch: "
              << compressed_length << " (" << indexType(compressed_index) << "), "
              << dynamic_length << " (" << indexType(dynamic_index) << ")" << std::endl;
  }

  std::cout << "Found " << queries.size() << " ranges of total length " << compressed_length << std::endl;
  std::cout << std::endl;
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
size_type
extractBenchmark(const GBWTType& index, size_type extract_queries)
{
  double start = readTimer();
  std::mt19937_64 rng(RANDOM_SEED);
  size_type total_length = 0;
  for(size_type i = 0; i < extract_queries; i++)
  {
    std::vector<node_type> sequence = index.extract(rng() % index.sequences());
    total_length += sequence.size();
  }
  double seconds = readTimer() - start;
  printTimeLength(indexType(index), extract_queries, total_length, seconds);
  return total_length;
}

void
extractBenchmark(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, size_type extract_queries)
{
  std::cout << "extract() benchmarks:" << std::endl;
  size_type compressed_length = extractBenchmark(compressed_index, extract_queries);
  size_type dynamic_length = extractBenchmark(dynamic_index, extract_queries);
  if(compressed_length != dynamic_length)
  {
    std::cerr << "extractBenchmark(): Total length " << compressed_length << " (" << indexType(compressed_index) << "), "
      << dynamic_length << " (" << indexType(dynamic_index) << ")" << std::endl;
  }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
