/*
  Copyright (c) 2017, 2018, 2019, 2020 Jouni Siren

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

#include <fstream>
#include <random>
#include <string>
#include <unistd.h>

#include <gbwt/dynamic_gbwt.h>
#include <gbwt/fast_locate.h>

using namespace gbwt;

//------------------------------------------------------------------------------

const std::string tool_name = "GBWT benchmark";

const size_type RANDOM_SEED  = 0xDEADBEEF;

void printUsage(int exit_code = EXIT_SUCCESS);

void compareIndexes(const GBWT& first, const GBWT& second, const std::string& first_name, const std::string& second_name);

void extendedStatistics(const DynamicGBWT& index);

std::vector<SearchState> findBenchmark(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const FastLocate& r_index, size_type find_queries, size_type pattern_length, std::vector<vector_type>& queries, std::vector<size_type>& first_occs);

void bidirectionalBenchmark(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const std::vector<vector_type>& queries);

template<class GBWTType>
void locateBenchmark(const GBWTType& index, const std::vector<SearchState>& queries, const std::vector<size_type>& first_occs);
template<>
void locateBenchmark<FastLocate>(const FastLocate& index, const std::vector<SearchState>& queries, const std::vector<size_type>& first_occs);

void extractBenchmark(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, size_type extract_queries);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2) { printUsage(); }

  int c = 0;
  bool compare = false, find = false, locate = false, extract = false, statistics = false, breakdown = false;
  bool use_r_index = false;
  size_type find_queries = 0, pattern_length = 0, extract_queries = 0;
  std::string compare_base;
  while((c = getopt(argc, argv, "c:f:p:le:rsS")) != -1)
  {
    switch(c)
    {
    case 'c':
      compare = true;
      compare_base = optarg;
      break;
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
    case 'r':
      use_r_index = true; break;
    case 's':
      statistics = true;
      break;
    case 'S':
      breakdown = true;
      break;
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
  sdsl::simple_sds::load_from(compressed_index, index_base + GBWT::EXTENSION);
  printStatistics(compressed_index, index_base);

  if(breakdown)
  {
    std::string filename = index_base + ".html";
    std::cout << "Writing space breakdown to " << filename << std::endl;
    std::ofstream out(filename, std::ios_base::binary);
    if(!out)
    {
      std::cerr << "benchmark: Cannot open output file " << filename << std::endl;
      std::exit(EXIT_FAILURE);
    }
    sdsl::write_structure<sdsl::HTML_FORMAT>(compressed_index, out);
    out.close();
    std::cout << std::endl;
  }

  if(compare)
  {
    GBWT second_index;
    sdsl::simple_sds::load_from(second_index, compare_base + GBWT::EXTENSION);
    printStatistics(second_index, compare_base);
    compareIndexes(compressed_index, second_index, index_base, compare_base);
  }

  if(!(find || locate || extract || statistics)) { return 0; }

  DynamicGBWT dynamic_index(compressed_index);
  printStatistics(dynamic_index, index_base);

  if(statistics) { extendedStatistics(dynamic_index); }

  if(!(find || locate || extract)) { return 0; }

  FastLocate r_index;
  if(use_r_index)
  {
    if(!sdsl::load_from_file(r_index, index_base + FastLocate::EXTENSION))
    {
      std::cerr << "benchmark: Cannot load the r-index from " << (index_base + FastLocate::EXTENSION) << std::endl;
      std::exit(EXIT_FAILURE);
    }
    r_index.setGBWT(compressed_index);
    printStatistics(r_index, index_base);
  }

  if(find)
  {
    std::vector<vector_type> queries;
    std::vector<size_type> first_occs;
    std::vector<SearchState> result = findBenchmark(compressed_index, dynamic_index, r_index, find_queries, pattern_length, queries, first_occs);
    if(compressed_index.bidirectional())
    {
      bidirectionalBenchmark(compressed_index, dynamic_index, queries);
    }
    queries.clear();
    if(locate)
    {
      locateBenchmark(compressed_index, result, first_occs);
      locateBenchmark(dynamic_index, result, first_occs);
      if(use_r_index) { locateBenchmark(r_index, result, first_occs); }
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
  std::cerr << "  -c X  Compare to the index with base name X" << std::endl;
  std::cerr << "  -f N  Benchmark N find() queries (requires -p)" << std::endl;
  std::cerr << "  -p N  Use patterns of length N" << std::endl;
  std::cerr << "  -l    Benchmark locate() queries (requires -f)" << std::endl;
  std::cerr << "  -e N  Benchmark N extract() queries" << std::endl;
  std::cerr << "  -r    Also use the r-index for find()/locate() queries" << std::endl;
  std::cerr << "  -s    Print extended statistics" << std::endl;
  std::cerr << "  -S    Write size breakdown to index_base.html" << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

//------------------------------------------------------------------------------

void
compareIndexes(const GBWT& first, const GBWT& second, const std::string& first_name, const std::string& second_name)
{
  std::cout << "Headers" << std::endl;
  if(first.header != second.header)
  {
    printHeader(first_name); std::cout << first.header << std::endl;
    printHeader(second_name); std::cout << second.header << std::endl;
  }
  std::cout << std::endl;

  // FIXME bwt

  std::cout << "Samples" << std::endl;
  if(first.da_samples.size() != second.da_samples.size() || first.da_samples.records() != second.da_samples.records())
  {
    printHeader(first_name); std::cout << first.da_samples.size() << " samples in " << first.da_samples.records() << " records" << std::endl;
    printHeader(second_name); std::cout << second.da_samples.size() << " samples in " << second.da_samples.records() << " records" << std::endl;
  }
  for(comp_type comp = 0; comp < first.effective() && comp < second.effective(); comp++)
  {
    if(first.da_samples.isSampled(comp) != second.da_samples.isSampled(comp))
    {
      printHeader(first_name); std::cout << "Record " << comp << " has samples: " << first.da_samples.isSampled(comp) << std::endl;
      printHeader(second_name); std::cout << "Record " << comp << " has samples: " << second.da_samples.isSampled(comp) << std::endl;
      break;
    }
  }
  // FIXME other structures in DASamples
  for(size_type i = 0; i < first.da_samples.size() && i < second.da_samples.size(); i++)
  {
    if(first.da_samples.array[i] != second.da_samples.array[i])
    {
      printHeader(first_name); std::cout << "Sample " << i << ": " << first.da_samples.array[i] << std::endl;
      printHeader(second_name); std::cout << "Sample " << i << ": " << second.da_samples.array[i] << std::endl;
      break;
    }
  }
  std::cout << std::endl;

  std::cout << "Metadata" << std::endl;
  if(first.metadata != second.metadata)
  {
    printHeader(first_name); std::cout << first.metadata << std::endl;
    printHeader(second_name); std::cout << second.metadata << std::endl;
  }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------

void
printStatistics(const std::string& header, const std::map<size_type, size_type>& distribution,
                size_type endmarker, size_type total, size_type n)
{
  std::cout << header << std::endl;
  printHeader("Total"); std::cout << total << std::endl;
  printHeader("Average"); std::cout << (total / static_cast<double>(n)) << std::endl;
  printHeader("Endmarker"); std::cout << endmarker << std::endl;

  double multiplier = 0.0;
  auto iter = distribution.begin();
  size_type records_seen = 0, value = 0;
  for(size_type nines = 1; nines <= 6; nines++)
  {
    multiplier = multiplier / 10.0 + 0.9;
    double limit = multiplier * n;
    while(iter != distribution.end() && static_cast<double>(records_seen) < limit)
    {
      value = iter->first;
      records_seen += iter->second;
      ++iter;
    }
    printHeader(std::to_string(multiplier)); std::cout << "at most " << value << std::endl;
  }
  while(iter != distribution.end())
  {
    value = iter->first;
    ++iter;
  }
  printHeader("Maximum"); std::cout << value << std::endl;

  std::cout << std::endl;
}

void
extendedStatistics(const DynamicGBWT& index)
{
  if(index.effective() < 2) { return; }

  std::map<size_type, size_type> run_distribution, outdegree_distribution;
  size_type endmarker_runs = 0, endmarker_outdegree = 0;
  size_type total_runs = 0, total_outdegree = 0;

  {
    const DynamicRecord& endmarker = index.record(ENDMARKER);
    endmarker_runs = endmarker.runs().first;
    endmarker_outdegree = endmarker.outdegree();
  }
  for(comp_type comp = 1; comp < index.effective(); comp++)
  {
    const DynamicRecord& record = index.record(index.toNode(comp));
    size_type record_runs = record.runs().first;
    run_distribution[record_runs]++;
    outdegree_distribution[record.outdegree()]++;
    total_runs += record_runs;
    total_outdegree += record.outdegree();
  }

  printStatistics("Runs", run_distribution, endmarker_runs, total_runs, index.effective() - 1);
  printStatistics("Outdegrees", outdegree_distribution, endmarker_outdegree, total_outdegree, index.effective() - 1);
}

//------------------------------------------------------------------------------

std::vector<vector_type>
generateQueries(const DynamicGBWT& index, size_type find_queries, size_type pattern_length)
{
  std::vector<vector_type> result;
  std::mt19937_64 rng(RANDOM_SEED);

  size_type attempts = 0;
  while(result.size() < find_queries && attempts < 10 * find_queries)
  {
    node_type node = index.toNode(rng() % (index.effective() - 1) + 1);
    if(index.empty(node)) { attempts++; continue; }
    edge_type position(node, rng() % index.nodeSize(node));
    vector_type sequence = index.extract(position, pattern_length);
    if(sequence.size() == pattern_length) { result.push_back(sequence); }
    attempts++;
  }

  std::cout << "Generated " << result.size() << " queries of total length " << (result.size() * pattern_length) << std::endl;
  return result;
}

//------------------------------------------------------------------------------

template<class GBWTType>
size_type
findBenchmark(const GBWTType& index, const std::vector<vector_type>& queries, std::vector<SearchState>& result, std::vector<size_type>&)
{
  double start = readTimer();
  size_type total_length = 0, total_size = 0;
  for(size_type i = 0; i < queries.size(); i++)
  {
    result[i] = index.find(queries[i].begin(), queries[i].end());
    total_length += queries[i].size();
    total_size += result[i].size();
  }
  double seconds = readTimer() - start;
  printTimeLength(indexType(index), queries.size(), total_length, seconds);
  return total_size;
}    

template<>
size_type
findBenchmark<FastLocate>(const FastLocate& index, const std::vector<vector_type>& queries, std::vector<SearchState>& result, std::vector<size_type>& first_occs)
{
  double start = readTimer();
  size_type total_length = 0, total_size = 0;
  for(size_type i = 0; i < queries.size(); i++)
  {
    result[i] = index.find(queries[i].begin(), queries[i].end(), first_occs[i]);
    total_length += queries[i].size();
    total_size += result[i].size();
  }
  double seconds = readTimer() - start;
  printTimeLength(indexType(index), queries.size(), total_length, seconds);
  return total_size;
}    

std::vector<SearchState>
findBenchmark(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const FastLocate& r_index, size_type find_queries, size_type pattern_length, std::vector<vector_type>& queries, std::vector<size_type>& first_occs)
{
  std::cout << "find() benchmarks:" << std::endl;

  bool use_r_index = !(r_index.empty());

  queries = generateQueries(dynamic_index, find_queries, pattern_length);
  std::vector<SearchState> result(queries.size());
  if(use_r_index) { first_occs.resize(queries.size()); }

  size_type compressed_length = findBenchmark(compressed_index, queries, result, first_occs);
  size_type dynamic_length = findBenchmark(dynamic_index, queries, result, first_occs);

  size_type r_length = compressed_length;
  if(use_r_index)
  {
    r_length = findBenchmark(r_index, queries, result, first_occs);
  }

  if(compressed_length != dynamic_length || compressed_length != r_length)
  {
    std::cerr << "findBenchmark(): Total length mismatch: "
              << compressed_length << " (" << indexType(compressed_index) << "), "
              << dynamic_length << " (" << indexType(dynamic_index) << ")";
    if(use_r_index)
    {
      std::cerr << ", " << r_length << "(" << indexType(r_index) << ")";
    }
    std::cerr << std::endl;
  }

  std::cout << "Found " << result.size() << " ranges of total length " << compressed_length << std::endl;
  std::cout << std::endl;

  return result;
}

//------------------------------------------------------------------------------

template<class GBWTType>
BidirectionalState
bidirectionalSearch(const GBWTType& index, const vector_type& query)
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
bidirectionalBenchmark(const GBWTType& index, const std::vector<vector_type>& queries)
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
bidirectionalBenchmark(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const std::vector<vector_type>& queries)
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
locateBenchmark(const GBWTType& index, const std::vector<SearchState>& queries, const std::vector<size_type>&)
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

template<>
void
locateBenchmark<FastLocate>(const FastLocate& index, const std::vector<SearchState>& queries, const std::vector<size_type>& first_occs)
{
  std::cout << "locate() benchmarks (" << indexType(index) << "):" << std::endl;

  {
    double start = readTimer();
    size_type found = 0;
    for(size_type i = 0; i < queries.size(); i++)
    {
      if(queries[i].empty()) { continue; }
      std::vector<size_type> result = index.locate(queries[i]);
      found += result.size();
    }
    double seconds = readTimer() - start;
    printTime("Slow", found, seconds);
  }

  {
    double start = readTimer();
    size_type found = 0;
    for(size_type i = 0; i < queries.size(); i++)
    {
      if(queries[i].empty()) { continue; }
      std::vector<size_type> result = index.locate(queries[i], first_occs[i]);
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
    vector_type sequence = index.extract(rng() % index.sequences());
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
