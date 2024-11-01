/*
  Copyright (c) 2018, 2019, 2020, 2021, 2023 Jouni Siren
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

#include <set>
#include <unistd.h>

#include <atomic>

#include <gbwt/dynamic_gbwt.h>
#include <gbwt/test.h>
#include <gbwt/variants.h>

using namespace gbwt;

//------------------------------------------------------------------------------

const std::string tool_name = "GBWT construction";

constexpr size_type MAX_ERRORS   = 100; // Do not print more error messages.
size_type errors                 = 0;

std::vector<SearchState> verifyFind(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const std::string& query_base, std::vector<vector_type>& queries);
void verifyBidirectional(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const std::vector<vector_type>& queries, const std::vector<SearchState>& find_results);
void verifyLocate(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const std::vector<SearchState>& queries);
void verifyExtract(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const std::string& base_name, bool both_orientations);
void verifyInverseLF(const GBWT& compressed_index, const DynamicGBWT& dynamic_index);
void verifySamples(const GBWT& compressed_index, const DynamicGBWT& dynamic_index);

//------------------------------------------------------------------------------

struct Config
{
  Config(int argc, char** argv);

  // Returns `true` if ok.
  bool validate() const;

  static void usage();

  // GBWT construction parameters.
  size_type batch_size = DynamicGBWT::INSERT_BATCH_SIZE / MILLION;
  size_type sample_interval = DynamicGBWT::SAMPLE_INTERVAL;

  // Special options.
  bool verify_index = false;
  bool both_orientations = false;
  bool build_index = true;
  bool build_empty = false;
  bool resample = false;

  // Input.
  enum input_type { input_sdsl, input_parse, input_text };
  input_type input = input_sdsl;
  std::vector<std::string> input_files;
  std::string index_base;

  // Build from parse.
  bool skip_overlaps = false;
  bool check_overlaps = false;
  std::set<std::string> phasing_files;

  // Output.
  std::string output_base;
  bool sdsl_format = false;
};

//------------------------------------------------------------------------------

struct MetadataConfig
{
  // Take config from the GBWT index.
  void from(const DynamicGBWT& index);

  // Use this metadata in the index.
  void applyTo(DynamicGBWT& index) const;

  // Add new contigs to the index.
  void addContigs(DynamicGBWT& index) const;

  bool need_sample_names = true;
  bool have_sample_names = false;
  bool use_contig_names = true;
  bool use_path_names = true;
  size_type contig_id = 0;

  std::set<size_type> samples;
  std::set<range_type> haplotypes;
  std::vector<std::string> sample_names, contig_names;
};

//------------------------------------------------------------------------------

size_type insertParse(DynamicGBWT& index, const std::string& filename, const Config& config, MetadataConfig& metadata);
size_type insertTextFile(DynamicGBWT& index, const std::string& filename, const Config& config);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  Verbosity::set(Verbosity::FULL);
  Config config(argc, argv);
  if(!(config.validate())) { std::exit(EXIT_FAILURE); }
  std::string gbwt_name = config.output_base + DynamicGBWT::EXTENSION;

  Version::print(std::cout, tool_name);

  if(!(config.index_base.empty())) { printHeader("Index name"); std::cout << config.index_base << std::endl; }
  printHeader("Input files"); std::cout << config.input_files.size();
  if(config.input == Config::input_parse)
  {
    std::cout << " (VCF parses";
    if(config.check_overlaps) { std::cout << "; checking overlaps"; }
    if(config.skip_overlaps) { std::cout << "; skipping overlaps"; }
    if(!(config.phasing_files.empty())) { std::cout << "; using " << config.phasing_files.size() << " phasing files"; }
    std::cout << ")";
  }
  std::cout << std::endl;
  printHeader("Output name"); std::cout << config.output_base << (config.sdsl_format ? " (SDSL format)" : " (simple-sds format)") << std::endl;
  if(config.batch_size != 0) { printHeader("Batch size"); std::cout << config.batch_size << " million" << std::endl; }
  printHeader("Orientation"); std::cout << (config.both_orientations ? "both" : "forward only") << std::endl;
  printHeader("Sample interval"); std::cout << config.sample_interval << std::endl;
  std::cout << std::endl;

  if(config.build_index)
  {
    double start = readTimer();

    // Load the index and determine whether we should use sample/contig/path names.
    // We assume that each VCF parse adds new contigs for the same samples.
    // We take the sample names either from the index we load or from the first parse.
    DynamicGBWT dynamic_index;
    MetadataConfig metadata;
    if(config.index_base.empty())
    {
      if(config.input == Config::input_parse) { dynamic_index.addMetadata(); }
    }
    else
    {
      sdsl::simple_sds::load_from(dynamic_index, config.index_base + DynamicGBWT::EXTENSION);
      printStatistics(dynamic_index, config.index_base);
      metadata.from(dynamic_index);
    }

    size_type input_size = 0;
    if(config.build_empty)
    {
      std::cout << "Building an empty GBWT" << std::endl;
    }
    for(const std::string& input_base : config.input_files)
    {
      if(config.build_empty) { continue; }
      printHeader("Input name"); std::cout << input_base << std::endl;
      if(config.input == Config::input_parse)
      {
        input_size += insertParse(dynamic_index, input_base, config, metadata);
      }
      else if(config.input == Config::input_sdsl)
      {
        text_buffer_type input(input_base);
        input_size += input.size() * (config.both_orientations ? 2 : 1);
        dynamic_index.insert(input, config.batch_size * MILLION, config.both_orientations, config.sample_interval);
      }
      else if(config.input == Config::input_text)
      {
        input_size += insertTextFile(dynamic_index, input_base, config);
      }
    }
    std::cout << std::endl;

    // Set metadata if we built from parse.
    if(config.input == Config::input_parse)
    {
      if(config.index_base.empty())
      {
        metadata.applyTo(dynamic_index);
      }
      else if(dynamic_index.hasMetadata())
      {
        metadata.addContigs(dynamic_index);
      }
    }

    // Serialize and print statistics.
    if(config.sdsl_format)
    {
      if(!sdsl::store_to_file(dynamic_index, gbwt_name))
      {
        std::cerr << "build_gbwt: Cannot write the index to " << gbwt_name << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    else { sdsl::simple_sds::serialize_to(dynamic_index, gbwt_name); }
    printStatistics(dynamic_index, config.output_base);

    double seconds = readTimer() - start;

    std::cout << "Indexed " << input_size << " nodes in " << seconds << " seconds (" << (input_size / seconds) << " nodes/second)" << std::endl;
    std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
    std::cout << std::endl;
  }

  if(config.resample)
  {
    std::string input_base = config.input_files.front();
    std::cout << "Resampling the index..." << std::endl;
    double resample_start = readTimer();
    std::cout << std::endl;

    GBWT compressed_index;
    sdsl::simple_sds::load_from(compressed_index, input_base + GBWT::EXTENSION);
    compressed_index.resample(config.sample_interval);
    printStatistics(compressed_index, config.output_base);

    if(config.sdsl_format)
    {
      if(!sdsl::store_to_file(compressed_index, gbwt_name))
      {
        std::cerr << "build_gbwt: Cannot write the index to " << gbwt_name << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    else { sdsl::simple_sds::serialize_to(compressed_index, gbwt_name); }

    double resample_seconds = readTimer() - resample_start;
    std::cout << "Resampled the index in " << resample_seconds << " seconds" << std::endl;
    std::cout << std::endl;
  }

  if(config.verify_index)
  {
    std::string input_base = config.input_files.front();
    std::cout << "Verifying the index..." << std::endl;
    double verify_start = readTimer();
    std::cout << std::endl;

    GBWT compressed_index;
    sdsl::simple_sds::load_from(compressed_index, gbwt_name);
    DynamicGBWT dynamic_index(compressed_index);

    std::vector<vector_type> queries;
    std::vector<SearchState> result = verifyFind(compressed_index, dynamic_index, input_base, queries);
    if(config.both_orientations)
    {
      verifyBidirectional(compressed_index, dynamic_index, queries, result);
    }
    verifyLocate(compressed_index, dynamic_index, result);
    verifyExtract(compressed_index, dynamic_index, input_base, config.both_orientations);
    if(config.both_orientations)
    {
      verifyInverseLF(compressed_index, dynamic_index);
    }
    verifySamples(compressed_index, dynamic_index);

    double verify_seconds = readTimer() - verify_start;
    if(errors > 0) { std::cout << "Index verification failed" << std::endl; }
    else { std::cout << "Index verified in " << verify_seconds << " seconds" << std::endl; }
    std::cout << std::endl;
  }

  return 0;
}

//------------------------------------------------------------------------------

Config::Config(int argc, char** argv)
{
  if(argc < 2) { usage(); std::exit(EXIT_FAILURE); }

  int c = 0;
  while((c = getopt(argc, argv, "b:cefF:i:lL:o:OpP:rRs:Stv")) != -1)
  {
    switch(c)
    {
    case 'b':
      this->batch_size = std::stoul(optarg); break;
    case 'c':
      this->check_overlaps = true; break;
    case 'e':
      this->build_empty = true; break;
    case 'f':
      this->both_orientations = false; break;
    case 'F':
      readRows(optarg, this->input_files, true); break;
    case 'i':
      this->index_base = optarg; break;
    case 'l':
      this->build_index = false; break;
    case 'L':
      {
        std::vector<std::string> rows;
        readRows(optarg, rows, true);
        this->phasing_files.insert(rows.begin(), rows.end());
      }
      break;
    case 'o':
      this->output_base = optarg; break;
    case 'O':
      this->sdsl_format = true; break;
    case 'p':
      this->input = input_parse; break;
    case 'P':
      this->phasing_files.insert(optarg); break;
    case 'r':
      this->both_orientations = true; break;
    case 'R':
      this->resample = true; this->build_index = false; break;
    case 's':
      this->sample_interval = std::stoul(optarg); break;
    case 'S':
      this->skip_overlaps = true; break;
    case 't':
      this->input = input_text; break;
    case 'v':
      this->verify_index = true; break;
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }

  // Append the remaining arguments to the list of input files and check the parameters.
  while(optind < argc) { this->input_files.push_back(argv[optind]); optind++; }

  // If we have a single input and no specified index / output, use that as the
  // base name for output.
  if(this->index_base.empty() && this->output_base.empty() && this->input_files.size() == 1)
  {
    this->output_base = this->input_files.front();
  }
}

bool
Config::validate() const
{
  if(this->input_files.empty() || this->output_base.empty()) { return false; }
  if(this->resample && this->input_files.size() != 1)
  {
    std::cerr << "build_gbwt: Resampling requires a single input" << std::endl;
    return false;
  }
  if(this->verify_index && !(this->input_files.size() == 1 && this->index_base.empty() && this->input == input_sdsl))
  {
    std::cerr << "build_gbwt: Verification only works with indexes for a single SDSL input" << std::endl;
    return false;
  }
  if(this->build_empty)
  {
    if(!(this->index_base.empty()))
    {
      std::cerr << "build_gbwt: Cannot load an index when building an empty GBWT" << std::endl;
      return false;
    }
    if(this->output_base.empty())
    {
      std::cerr << "build_gbwt: No output file specified for empty GBWT" << std::endl;
      return false;
    }
    if(this->input == Config::input_parse)
    {
      std::cerr << "build_gbwt: Cannot build empty GBWT from a parsed VCF file" << std::endl;
      return false;
    }
    if(this->verify_index)
    {
      std::cerr << "build_gbwt: Cannot verify an empty GBWT" << std::endl;
      return false;
    }
  }

  return true;
}

void
Config::usage()
{
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: build_gbwt [options] input1 [input2 ...]" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Input / output:" << std::endl;
  std::cerr << "  -e    Build an empty GBWT (the input will not be read)" << std::endl;
  std::cerr << "  -F X  Read a list of input files from X, one file per line (no input1 needed; may repeat)" << std::endl;
  std::cerr << "  -i X  Insert the sequences into an existing index with base name X" << std::endl;
  std::cerr << "  -l    Load an existing index instead of building it" << std::endl;
  std::cerr << "  -R    Resample sequence ids in the loaded index (implies -l)" << std::endl;
  std::cerr << "  -o X  Use base name X for output (default: the only input)" << std::endl;
  std::cerr << "  -O    Output SDSL format instead of simple-sds format" << std::endl;
  std::cerr << "  -t    The input is a text file (each line is a comma-separated sequence of nodes)" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Construction parameters:" << std::endl;
  std::cerr << "  -b N  Insert in batches of N million nodes (default: " << (DynamicGBWT::INSERT_BATCH_SIZE / MILLION) << ")" << std::endl;
  std::cerr << "  -f    Index the sequences only in forward orientation (default)" << std::endl;
  std::cerr << "  -r    Index the sequences also in reverse orientation" << std::endl;
  std::cerr << "  -s N  Sample sequence ids at one out of N positions (default: " << DynamicGBWT::SAMPLE_INTERVAL << "; use 0 for no samples)" << std::endl;
  std::cerr << "  -v    Verify the index after construction" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Construction from a parsed VCF file:" << std::endl;
  std::cerr << "  -p    The input is a parsed VCF file" << std::endl;
  std::cerr << "  -P X  Only use the phasing information in file X (may repeat)" << std::endl;
  std::cerr << "  -L X  Read a list of phasing files from X, one file per line (may repeat)" << std::endl;
  std::cerr << "  -c    Check for overlapping variants in haplotypes" << std::endl;
  std::cerr << "  -S    Skip overlapping variants" << std::endl;
  std::cerr << std::endl;
}

//------------------------------------------------------------------------------

void
MetadataConfig::from(const DynamicGBWT& index)
{
  this->need_sample_names = false;
  this->use_contig_names = (index.hasMetadata() && index.metadata.hasContigNames());
  if(index.hasMetadata()) { this->contig_id = index.metadata.contigs(); }
  this->use_path_names = (index.hasMetadata() && index.metadata.hasPathNames());
}

void
MetadataConfig::applyTo(DynamicGBWT& index) const
{
  if(this->have_sample_names) { index.metadata.setSamples(this->sample_names); }
  else { index.metadata.setSamples(this->samples.size()); }

  index.metadata.setHaplotypes(this->haplotypes.size());

  if(this->use_contig_names) { index.metadata.setContigs(this->contig_names); }
  else { index.metadata.setContigs(this->contig_id); }
}

void
MetadataConfig::addContigs(DynamicGBWT& index) const
{
  if(this->use_contig_names)
  {
    index.metadata.addContigs(this->contig_names);
  }
  else
  {
    index.metadata.clearContigNames();
    index.metadata.setContigs(this->contig_id);
  }
}

//------------------------------------------------------------------------------

size_type
insertParse(DynamicGBWT& index, const std::string& filename, const Config& config, MetadataConfig& metadata)
{
  // Load the parse and determine if we still want to use sample/contig names.
  VariantPaths variants;
  if(!sdsl::load_from_file(variants, filename))
  {
    std::cerr << "build_gbwt: Cannot load variants from " << filename << std::endl;
    std::exit(EXIT_FAILURE);
  }
  metadata.need_sample_names &= variants.hasSampleNames();
  if(metadata.need_sample_names && !(metadata.have_sample_names))
  {
    metadata.sample_names = variants.getSampleNames();
    metadata.have_sample_names = true;
  }
  metadata.use_contig_names &= variants.hasContigName();
  if(metadata.use_contig_names)
  {
    metadata.contig_names.emplace_back(variants.getContigName());
  }

  // Build GBWT from the parse.
  if(config.check_overlaps) { checkOverlaps(variants, std::cerr, true); }
  std::set<range_type> overlaps;
  size_type node_width = variants.nodeWidth(config.both_orientations);
  size_type old_size = index.size();
  GBWTBuilder builder(node_width, config.batch_size * MILLION, config.sample_interval);
  builder.swapIndex(index);
  generateHaplotypes(variants, config.phasing_files,
    [](size_type) -> bool { return true; },
    [&](const Haplotype& haplotype)
    {
      builder.insert(haplotype.path, config.both_orientations);
      metadata.samples.insert(haplotype.sample);
      metadata.haplotypes.insert(range_type(haplotype.sample, haplotype.phase));
      if(metadata.use_path_names)
      {
        builder.index.metadata.addPath({
          static_cast<PathName::path_name_type>(haplotype.sample),
          static_cast<PathName::path_name_type>(metadata.contig_id),
          static_cast<PathName::path_name_type>(haplotype.phase),
          static_cast<PathName::path_name_type>(haplotype.count)
        });
      }
    },
    [&](size_type site, size_type allele) -> bool
    {
      if(config.check_overlaps) { overlaps.insert(range_type(site, allele)); }
      return config.skip_overlaps;
    });
  builder.finish();
  builder.swapIndex(index);
  if(config.check_overlaps && !overlaps.empty())
  {
    std::cerr << overlaps.size() << " unresolved overlaps:" << std::endl;
    for(range_type overlap : overlaps)
    {
      std::cerr << "- site " << overlap.first << ", allele " << overlap.second << std::endl;
    }
  }

  metadata.contig_id++;
  return index.size() - old_size;
}

size_type
insertTextFile(DynamicGBWT& index, const std::string& filename, const Config& config)
{
  // TODO: How to determine width?
  GBWTBuilder builder(64, config.batch_size * MILLION, config.sample_interval);
  builder.swapIndex(index);

  size_type input_size = 0;
  std::ifstream infile(filename, std::ios_base::binary);
  if(!infile)
  {
    std::cerr << "build_gbwt: Cannot open " << filename << " for reading" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string buffer;
  while(std::getline(infile, buffer))
  {
    vector_type sequence;
    size_t start = 0;
    while(start < buffer.length())
    {
      size_t limit = buffer.find(',', start);
      if(limit == std::string::npos) { limit = buffer.length(); }
      sequence.push_back(std::stoul(buffer.substr(start, limit - start)));
      start = limit + 1;
    }
    input_size += (sequence.size() + 1) * (config.both_orientations ? 2 : 1);
    builder.insert(sequence, config.both_orientations);
  }
  infile.close();

  builder.finish();
  builder.swapIndex(index);
  return input_size;
}

//------------------------------------------------------------------------------

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
  TODO We could validate the actual results in DynamicGBWT by extracting backwards
  using Psi(). Psi() can be implemented using the incoming edges.
*/

std::vector<SearchState>
verifyFind(const GBWT& compressed_index, const DynamicGBWT& dynamic_index, const std::string& query_base, std::vector<vector_type>& queries)
{
  std::cout << "Verifying find()..." << std::endl;

  double start = readTimer();
  size_type initial_errors = errors;
  queries = generateQueries(query_base, true);
  std::vector<SearchState> result(queries.size());

  for(size_type i = 0; i < queries.size(); i++)
  {
    result[i] = compressed_index.find(queries[i].begin(), queries[i].end());
    SearchState dynamic_result = dynamic_index.find(queries[i].begin(), queries[i].end());
    if(dynamic_result != result[i])
    {
      errors++;
      if(errors <= MAX_ERRORS)
      {
        std::cerr << "verifyFind(): Mismatching results with query " << i << std::endl;
        std::cerr << "verifyFind(): " << indexType(compressed_index) << ": " << result[i] << ", " << indexType(dynamic_index) << ": " << dynamic_result << std::endl;
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
  std::vector<size_type> result;
  for(size_type i = query.range.first; i <= query.range.second; i++)
  {
    result.push_back(index.locate(query.node, i));
  }
  removeDuplicates(result, false);
  for(size_type res : result) { hash ^= wang_hash_64(res); }
  return hash;
}

template<class GBWTType>
size_type
fastLocate(const GBWTType& index, SearchState query)
{
  size_type hash = 0;
  std::vector<size_type> result = index.locate(query);
  for(size_type res : result) { hash ^= wang_hash_64(res); }
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
  std::vector<size_type> offsets = startOffsets(base_name, true);
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
  Inverse LF(): Try inverting LF() at every position.
*/

void
verifyInverseLF(const GBWT& compressed_index, const DynamicGBWT& dynamic_index)
{
  std::cout << "Verifying inverse LF()..." << std::endl;

  double start = readTimer();
  size_type initial_errors = errors;

  #pragma omp parallel for schedule(dynamic, 1)
  for(size_type sequence = 0; sequence < compressed_index.sequences(); sequence++)
  {
    edge_type prev(ENDMARKER, sequence);
    edge_type curr = compressed_index.start(sequence);
    while(curr.first != ENDMARKER)
    {
      edge_type compressed_pred = compressed_index.inverseLF(curr);
      edge_type dynamic_pred = dynamic_index.inverseLF(curr);
      if(compressed_pred != prev || dynamic_pred != prev)
      {
        errors++;
        if(errors <= MAX_ERRORS)
        {
          std::cerr << "verifyInverseLF(): Mismatching results for sequence " << sequence << ", position " << curr << std::endl;
          std::cerr << "verifyInverseLF(): " << indexType(compressed_index) << ": " << compressed_pred << ", " << indexType(dynamic_index) << ": " << dynamic_pred << ", truth: " << prev << std::endl;
        }
      }
      prev = curr;
      curr = compressed_index.LF(curr);
    }
  }

  double seconds = readTimer() - start;
  if(errors > initial_errors) { std::cout << "Inverse LF() verification failed" << std::endl; }
  else { std::cout << "Inverse LF() verified in " << seconds << " seconds" << std::endl; }
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
