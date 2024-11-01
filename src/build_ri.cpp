/*
  Copyright (c) 2018, 2020, 2021 Jouni Siren
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

#include <string>
#include <unistd.h>

#include <gbwt/fast_locate.h>
#include <gbwt/test.h>

using namespace gbwt;

//------------------------------------------------------------------------------

const std::string tool_name = "R-index construction";

constexpr size_type MAX_ERRORS   = 100; // Do not print more error messages.
size_type errors                 = 0;

void printUsage(int exit_code = EXIT_SUCCESS);

typedef std::pair<SearchState, size_type> find_result;

std::vector<find_result> verifyFind(const GBWT& gbwt_index, const FastLocate& r_index, const std::string& base_name, std::vector<vector_type>& queries);
void verifyLocate(const GBWT& gbwt_index, const FastLocate& r_index, const std::vector<find_result>& queries);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2) { printUsage(); }
  Verbosity::set(Verbosity::FULL);

  std::string base_name;
  bool verify_index = false;
  int c = 0;
  while((c = getopt(argc, argv, "t:v")) != -1)
  {
    switch(c)
    {
    case 't':
      omp_set_num_threads(std::stoul(optarg)); break;
    case 'v':
      verify_index = true; break;
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }
  if(optind + 1 != argc)
  {
    printUsage(EXIT_FAILURE);
  }
  base_name = argv[optind];

  Version::print(std::cout, tool_name);
  printHeader("Base name"); std::cout << base_name << std::endl;
  std::cout << std::endl;

  double start = readTimer();

  std::string gbwt_name = base_name + GBWT::EXTENSION;
  GBWT gbwt_index;
  sdsl::simple_sds::load_from(gbwt_index, gbwt_name);
  printStatistics(gbwt_index, base_name);

  FastLocate r_index(gbwt_index);
  std::string r_index_name = base_name + FastLocate::EXTENSION;
  if(!sdsl::store_to_file(r_index, r_index_name))
  {
    std::cerr << "build_ri: Cannot write the r-index to " << r_index_name << std::endl;
    std::exit(EXIT_FAILURE);
  }
  printStatistics(r_index, base_name);

  double seconds = readTimer() - start;

  std::cout << "Indexed " << gbwt_index.size() << " nodes in " << seconds << " seconds (" << (gbwt_index.size() / seconds) << " nodes/second)" << std::endl;
  std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  if(verify_index)
  {
    std::cout << "Verifying the index..." << std::endl;
    double verify_start = readTimer();
    std::cout << std::endl;

    std::vector<vector_type> queries;
    std::vector<find_result> result = verifyFind(gbwt_index, r_index, base_name, queries);
    verifyLocate(gbwt_index, r_index, result);

    double verify_seconds = readTimer() - verify_start;
    if(errors > 0) { std::cout << "Index verification failed" << std::endl; }
    else { std::cout << "Index verified in " << verify_seconds << " seconds" << std::endl; }
    std::cout << std::endl;
  }

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: build_ri [options] base_name" << std::endl;
  std::cerr << "  -t N  Extract the samples using N threads" << std::endl;
  std::cerr << "  -v    Verify the index using the sequences in file 'base_name'" << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

//------------------------------------------------------------------------------

/*
  find() queries: Ensure that both index types give the same results.
*/

std::vector<find_result>
verifyFind(const GBWT& gbwt_index, const FastLocate& r_index, const std::string& base_name, std::vector<vector_type>& queries)
{
  std::cout << "Verifying find()..." << std::endl;

  double start = readTimer();
  size_type initial_errors = errors;
  queries = generateQueries(base_name, true);
  std::vector<find_result> result(queries.size());

  for(size_type i = 0; i < queries.size(); i++)
  {
    result[i].first = gbwt_index.find(queries[i].begin(), queries[i].end());
    SearchState r_index_result = r_index.find(queries[i].begin(), queries[i].end(), result[i].second);
    if(r_index_result != result[i].first)
    {
      errors++;
      if(errors <= MAX_ERRORS)
      {
        std::cerr << "verifyFind(): Mismatching results with query " << i << std::endl;
        std::cerr << "verifyFind(): Expected " << result[i].first << ", got " << r_index_result << std::endl;
      }
    }
  }

  double seconds = readTimer() - start;
  if(errors > initial_errors) { std::cout << "find() verification failed" << std::endl; }
  else { std::cout << "find() verified in " << seconds << " seconds" << std::endl; }
  std::cout << std::endl;

  return result;
}

//------------------------------------------------------------------------------

/*
  locate() queries: Ensure that both index types give the same results.
*/

size_type
totalLength(const std::vector<find_result>& states)
{
  size_type result = 0;
  for(find_result state : states) { result += state.first.size(); }
  return result;
}

void
verifyLocate(const GBWT& gbwt_index, const FastLocate& r_index, const std::vector<find_result>& queries)
{
  std::cout << "Verifying locate()..." << std::endl;

  double start = readTimer();
  size_type initial_errors = errors;
  std::cout << queries.size() << " ranges of total length " << totalLength(queries) << std::endl;
  std::vector<range_type> blocks = Range::partition(range_type(0, queries.size() - 1), 4 * omp_get_max_threads());

  #pragma omp parallel for schedule(dynamic, 1)
  for(size_type block = 0; block < blocks.size(); block++)
  {
    for(size_type i = blocks[block].first; i <= blocks[block].second; i++)
    {
      if(queries[i].first.empty()) { continue; }
      size_type gbwt_hash = 0, r_index_slow = 0, r_index_fast = 0;
      {
        std::vector<size_type> result = gbwt_index.locate(queries[i].first);
        for(size_type res : result) { gbwt_hash ^= wang_hash_64(res); }
      }
      {
        std::vector<size_type> result = r_index.locate(queries[i].first);
        for(size_type res : result) { r_index_slow ^= wang_hash_64(res); }
      }
      {
        std::vector<size_type> result = r_index.locate(queries[i].first, queries[i].second);
        for(size_type res : result) { r_index_fast ^= wang_hash_64(res); }
      }

      if(r_index_slow != gbwt_hash || r_index_fast != gbwt_hash)
      {
        #pragma omp critical
        {
          errors++;
          if(errors <= MAX_ERRORS)
          {
            std::cerr << "verifyLocate(): Hash mismatch with query " << i << std::endl;
            std::cerr << "verifyLocate(): Expected " << gbwt_hash << ", got " << r_index_slow << " (slow), " << r_index_fast << " (fast)" << std::endl;
          }
        }
      }
    }
  }

  double seconds = readTimer() - start;
  if(errors > initial_errors) { std::cout << "locate() verification failed" << std::endl; }
  else { std::cout << "locate() verified in " << seconds << " seconds" << std::endl; }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
