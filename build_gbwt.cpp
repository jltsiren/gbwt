/*
  Copyright (c) 2018 Jouni Siren
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

#include <atomic>

#include <gbwt/dynamic_gbwt.h>

using namespace gbwt;

//------------------------------------------------------------------------------

const std::string tool_name = "GBWT construction";

const size_type MAX_ERRORS   = 100; // Do not print more error messages.
size_type errors             = 0;

const size_type RANDOM_SEED  = 0xDEADBEEF;
const size_type QUERIES      = 20000;
const size_type QUERY_LENGTH = 60;

void printUsage(int exit_code = EXIT_SUCCESS);

std::vector<SearchState> verifyFind(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const std::string& query_base, std::vector<vector_type>& queries);
void verifyBidirectional(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const std::vector<vector_type>& queries, const std::vector<SearchState>& find_results);
void verifyLocate(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const std::vector<SearchState>& queries);
void verifyExtract(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const std::string& base_name, bool both_orientations);
void verifySamples(const GBWT& compressed_index, const DynamicGBWT& dynamic_index);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2) { printUsage(); }

  size_type batch_size = DynamicGBWT::INSERT_BATCH_SIZE / MILLION, sample_interval = DynamicGBWT::SAMPLE_INTERVAL;
  bool verify_index = false, both_orientations = false, build_index = true;
  std::string index_base, input_base, output_base;
  int c = 0;
  while((c = getopt(argc, argv, "b:fi:lo:rs:v")) != -1)
  {
    switch(c)
    {
    case 'b':
      batch_size = std::stoul(optarg); break;
    case 'f':
      both_orientations = false; break;
    case 'i':
      index_base = optarg; break;
    case 'l':
      build_index = false; break;
    case 'o':
      output_base = optarg; break;
    case 'r':
      both_orientations = true; break;
    case 's':
      sample_interval = std::stoul(optarg); break;
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
  if(!build_index && !verify_index)
  {
    std::cerr << "build_gbwt: Index can only be loaded for verification" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string gbwt_name = output_base + DynamicGBWT::EXTENSION;

  Version::print(std::cout, tool_name);

  if(!(index_base.empty())) { printHeader("Index name"); std::cout << index_base << std::endl; }
  printHeader("Input files"); std::cout << input_files << std::endl;
  printHeader("Output name"); std::cout << output_base << std::endl;
  if(batch_size != 0) { printHeader("Batch size"); std::cout << batch_size << " million" << std::endl; }
  printHeader("Orientation"); std::cout << (both_orientations ? "both" : "forward only") << std::endl;
  printHeader("Sample interval"); std::cout << sample_interval << std::endl;
  std::cout << std::endl;

  if(build_index)
  {
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
      input_size += input.size() * (both_orientations ? 2 : 1);
      dynamic_index.insert(input, batch_size * MILLION, both_orientations, sample_interval);
      optind++;
    }
    std::cout << std::endl;

    sdsl::store_to_file(dynamic_index, gbwt_name);
    printStatistics(dynamic_index, output_base);

    double seconds = readTimer() - start;

    std::cout << "Indexed " << input_size << " nodes in " << seconds << " seconds (" << (input_size / seconds) << " nodes/second)" << std::endl;
    std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
    std::cout << std::endl;
  }
  else { input_base = argv[optind]; }

  if(verify_index)
  {
    std::cout << "Verifying the index..." << std::endl;
    double verify_start = readTimer();
    std::cout << std::endl;

    GBWT compressed_index;
    sdsl::load_from_file(compressed_index, gbwt_name);
    DynamicGBWT dynamic_index;
    sdsl::load_from_file(dynamic_index, gbwt_name);

    std::vector<vector_type> queries;
    std::vector<SearchState> results = verifyFind(compressed_index, dynamic_index, input_base, queries);
    if(both_orientations)
    {
      verifyBidirectional(compressed_index, dynamic_index, queries, results);
    }
    verifyLocate(compressed_index, dynamic_index, results);
    verifyExtract(compressed_index, dynamic_index, input_base, both_orientations);
    verifySamples(compressed_index, dynamic_index);

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

  std::cerr << "Usage: build_gbwt [options] input1 [input2 ...]" << std::endl;
  std::cerr << "  -b N  Insert in batches of N million nodes (default: " << (DynamicGBWT::INSERT_BATCH_SIZE / MILLION) << ")" << std::endl;
  std::cerr << "  -f    Index the sequences only in forward orientation (default)" << std::endl;
  std::cerr << "  -i X  Insert the sequences into an existing index with base name X" << std::endl;
  std::cerr << "  -l    Load an existing index instead of building it" << std::endl;
  std::cerr << "  -o X  Use base name X for output (default: the only input)" << std::endl;
  std::cerr << "  -r    Index the sequences also in reverse orientation" << std::endl;
  std::cerr << "  -s N  Sample sequence ids at one out of N positions (default: " << DynamicGBWT::SAMPLE_INTERVAL << "; use 0 for no samples)" << std::endl;
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

std::vector<vector_type>
generateQueries(const std::string& base_name)
{
  std::mt19937_64 rng(RANDOM_SEED);
  std::vector<vector_type> result;
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
    vector_type candidate(text.begin() + start_offset, text.begin() + start_offset + QUERY_LENGTH);
    bool ok = true;
    for(auto node : candidate)
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

/*
  find() queries: Ensure that both index types give the same results.
*/

std::vector<SearchState>
verifyFind(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const std::string& query_base, std::vector<vector_type>& queries)
{
  std::cout << "Verifying find()..." << std::endl;

  double start = readTimer();
  size_type initial_errors = errors;
  queries = generateQueries(query_base);
  std::vector<SearchState> results(queries.size());

  for(size_type i = 0; i < queries.size(); i++)
  {
    results[i] = compressed_index.find(queries[i].begin(), queries[i].end());
    SearchState dynamic_result = dynamic_index.find(queries[i].begin(), queries[i].end());
    if(dynamic_result != results[i])
    {
      errors++;
      if(errors <= MAX_ERRORS)
      {
        std::cerr << "verifyFind(): Mismatching results with query " << i << std::endl;
        std::cerr << "verifyFind(): " << indexType(compressed_index) << ": " << results[i] << ", " << indexType(dynamic_index) << ": " << dynamic_result << std::endl;
      }
    }
  }

  double seconds = readTimer() - start;
  if(errors > initial_errors) { std::cout << "find() verification failed" << std::endl; }
  else { std::cout << "find() verified in " << seconds << " seconds" << std::endl; }
  std::cout << std::endl;

  return results;
}

//------------------------------------------------------------------------------

/*
  Bidirectional queries: Ensure that both index types give the same results and that the
  results match those returned by find() queries.
*/

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

void
verifyBidirectional(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const std::vector<vector_type>& queries, const std::vector<SearchState>& find_results)
{
  std::cout << "Verifying bidirectional search..." << std::endl;

  double start = readTimer();
  size_type initial_errors = errors;

  for(size_type i = 0; i < queries.size(); i++)
  {
    BidirectionalState compressed_state = bidirectionalSearch(compressed_index, queries[i]);
    BidirectionalState dynamic_state = bidirectionalSearch(dynamic_index, queries[i]);
    vector_type reverse_query;
    reversePath(queries[i], reverse_query);
    BidirectionalState find_state(find_results[i], dynamic_index.find(reverse_query.begin(), reverse_query.end()));
    if(compressed_state != dynamic_state || compressed_state != find_state)
    {
      errors++;
      if(errors <= MAX_ERRORS)
      {
        std::cerr << "verifyBidirectional(): Mismatching results with query " << i << std::endl;
        std::cerr << "verifyBidirectional(): " << indexType(compressed_index) << ": " << compressed_state << ", " << indexType(dynamic_index) << ": " << dynamic_state << ", find(): " << find_state << std::endl;
      }
    }
  }

  double seconds = readTimer() - start;
  if(errors > initial_errors) { std::cout << "Bidirectional search verification failed" << std::endl; }
  else { std::cout << "Bidirectional search verified in " << seconds << " seconds" << std::endl; }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------

/*
  locate() queries: Ensure that both index types and both algorithms give the same results.
*/

template<class GBWTType>
size_type
directLocate(const GBWTType& index, SearchState query)
{
  size_type hash = 0;
  std::vector<size_type> results;
  for(size_type i = query.range.first; i <= query.range.second; i++)
  {
    results.push_back(index.locate(query.node, i));
  }
  removeDuplicates(results, false);
  for(size_type res : results) { hash ^= wang_hash_64(res); }
  return hash;
}

template<class GBWTType>
size_type
fastLocate(const GBWTType& index, SearchState query)
{
  size_type hash = 0;
  std::vector<size_type> results = index.locate(query);
  for(size_type res : results) { hash ^= wang_hash_64(res); }
  return hash;
}

void
verifyLocate(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const std::vector<SearchState>& queries)
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
      SearchState query = queries[i];
      if(query.empty()) { continue; }

      size_type compressed_direct = directLocate(compressed_index, query);
      size_type dynamic_direct = directLocate(dynamic_index, query);
      size_type compressed_fast = fastLocate(compressed_index, query);
      size_type dynamic_fast = fastLocate(dynamic_index, query);

      if(compressed_direct != dynamic_direct || dynamic_direct != compressed_fast || compressed_fast != dynamic_fast)
      {
        #pragma omp critical
        {
          errors++;
          if(errors <= MAX_ERRORS)
          {
            std::cerr << "verifyLocate(): Hash mismatch with query " << i << std::endl;
            std::cerr << "verifyLocate(): " << compressed_direct << " (direct), " << compressed_fast << " (fast) in " << indexType(compressed_index) << std::endl;
            std::cerr << "verifyLocate(): " << dynamic_direct << " (direct), " << dynamic_fast << " (fast) in " << indexType(dynamic_index) << std::endl;
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

/*
  extract() queries: Ensure that the index contains the correct sequences.
*/

void
tryExtract(const GBWT& compressed_index, const DynamicGBWT& dynamic_index,
           text_buffer_type& text, const std::vector<size_type>& offsets,
           size_type sequence, bool both_orientations, bool is_reverse)
{
  if(is_reverse && !both_orientations) { return; }
  size_type seq_id = (both_orientations ? Path::encode(sequence, is_reverse) : sequence);

  // Extract the sequences.
  vector_type compressed_result = compressed_index.extract(seq_id);
  vector_type dynamic_result = dynamic_index.extract(seq_id);
  vector_type correct_sequence; correct_sequence.reserve(compressed_result.size());
  for(size_type i = offsets[sequence]; text[i] != ENDMARKER; i++) { correct_sequence.push_back(text[i]); }
  if(is_reverse) { reversePath(correct_sequence); }

  // Compare the lengths.
  if(compressed_result.size() != correct_sequence.size() || compressed_result.size() != dynamic_result.size())
  {
    #pragma omp critical
    {
      errors++;
      if(errors <= MAX_ERRORS)
      {
        std::cerr << "verifyExtract(): Length mismatch with sequence " << sequence << (is_reverse ? " (reverse)" : " (forward)") << std::endl;
        std::cerr << "verifyExtract(): Text: " << correct_sequence.size() << ", "
          << indexType(compressed_index) << ": " << compressed_result.size() << ", "
            << indexType(dynamic_index) << ": " << dynamic_result.size() << std::endl;
      }
    }
    return;
  }

  // Compare the sequences.
  for(size_type i = 0; i < compressed_result.size(); i++)
  {
    if(compressed_result[i] != correct_sequence[i] || compressed_result[i] != dynamic_result[i])
    {
      #pragma omp critical
      {
        errors++;
        if(errors <= MAX_ERRORS)
        {
          std::cerr << "verifyExtract(): Mismatch at sequence " << sequence << ", offset " << i << (is_reverse ? " (reverse)" : " (forward)") << std::endl;
          std::cerr << "verifyExtract(): Text: " << correct_sequence[i] << ", "
            << indexType(compressed_index) << ": " << compressed_result[i] << ", "
            << indexType(dynamic_index) << ": " << dynamic_result[i] << std::endl;
        }
      }
      return;
    }
  }
}

void
verifyExtract(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const std::string& base_name, bool both_orientations)
{
  std::cout << "Verifying extract()..." << std::endl;

  double start = readTimer();
  size_type initial_errors = errors;
  std::vector<size_type> offsets = startOffsets(base_name);
  size_type expected_sequences = offsets.size() * (both_orientations ? 2 : 1);
  if(compressed_index.sequences() != expected_sequences || compressed_index.sequences() != dynamic_index.sequences())
  {
    errors++;
    if(errors <= MAX_ERRORS)
    {
      std::cerr << "verifyExtract(): Mismatching number of sequences" << std::endl;
      std::cerr << "verifyExtract(): Input: " << expected_sequences << ", "
                << indexType(compressed_index) << ": " << compressed_index.sequences() << ", "
                << indexType(dynamic_index) << ": " << dynamic_index.sequences() << std::endl;
    }
    std::cout << "extract() verification failed" << std::endl;
    return;
  }
  std::vector<range_type> blocks = Range::partition(range_type(0, offsets.size() - 1), 4 * omp_get_max_threads());

  #pragma omp parallel for schedule(dynamic, 1)
  for(size_type block = 0; block < blocks.size(); block++)
  {
    text_buffer_type text(base_name);
    for(size_type sequence = blocks[block].first; sequence <= blocks[block].second; sequence++)
    {
      tryExtract(compressed_index, dynamic_index, text, offsets, sequence, both_orientations, false);
      tryExtract(compressed_index, dynamic_index, text, offsets, sequence, both_orientations, true);
    }
  }

  double seconds = readTimer() - start;
  if(errors > initial_errors) { std::cout << "extract() verification failed" << std::endl; }
  else { std::cout << "extract() verified in " << seconds << " seconds" << std::endl; }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------

/*
  Ensure that all samples have the correct sequence identifiers.
*/

template<class GBWTType>
bool
trySample(const GBWTType& index, size_type sequence, edge_type& current, std::atomic<size_type>& samples_found)
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
          std::cerr << "verifySamples(): " << indexType(index) << ": Verification failed with sequence " << sequence << ", position " << current << std::endl;
          std::cerr << "verifySamples(): Sample had sequence id " << sample << std::endl;
        }
      }
      return false;
    }
  }
  current = index.LF(current);
  return true;
}

void
verifySamples(const GBWT& compressed_index, const DynamicGBWT& dynamic_index)
{
  std::cout << "Verifying samples..." << std::endl;

  double start = readTimer();
  size_type initial_errors = errors;
  std::atomic<size_type> found_compressed(0), found_dynamic(0);
  std::vector<range_type> blocks = Range::partition(range_type(0, compressed_index.sequences() - 1), 4 * omp_get_max_threads());

  #pragma omp parallel for schedule(dynamic, 1)
  for(size_type block = 0; block < blocks.size(); block++)
  {
    for(size_type sequence = blocks[block].first; sequence <= blocks[block].second; sequence++)
    {
      edge_type curr_compressed(ENDMARKER, sequence), curr_dynamic(ENDMARKER, sequence);
      do
      {
        if(!trySample(compressed_index, sequence, curr_compressed, found_compressed)) { break; }
        if(!trySample(dynamic_index, sequence, curr_dynamic, found_dynamic)) { break; }
        if(curr_compressed != curr_dynamic)
        {
          #pragma omp critical
          {
            errors++;
            if(errors <= MAX_ERRORS)
            {
              std::cerr << "verifySamples(): Position mismatch between indexes" << std::endl;
              std::cerr << "verifySamples(): " << indexType(compressed_index) << ": " << curr_compressed << ", "
                << indexType(dynamic_index) << ": " << curr_dynamic << std::endl;
            }
          }
          break;
        }
      }
      while(curr_compressed.first != ENDMARKER);
    }
  }

  if(found_compressed != compressed_index.samples() || found_dynamic != dynamic_index.samples() || found_compressed != found_dynamic)
  {
    errors++;
    if(errors <= MAX_ERRORS)
    {
      std::cerr << "verifySamples(): Mismatch in the number of samples" << std::endl;
      std::cerr << "verifySamples(): " << indexType(compressed_index) << ": " << found_compressed << ", "
        << indexType(dynamic_index) << ": " << found_dynamic << std::endl;
    }
  }

  double seconds = readTimer() - start;
  if(errors > initial_errors) { std::cout << "Sample verification failed" << std::endl; }
  else { std::cout << "Samples verified in " << seconds << " seconds" << std::endl; }
  std::cout << std::endl;
}

//------------------------------------------------------------------------------
