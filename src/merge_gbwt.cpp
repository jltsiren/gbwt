/*
  Copyright (c) 2017, 2018, 2019, 2021 Jouni Siren
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

#include <gbwt/dynamic_gbwt.h>

using namespace gbwt;

//------------------------------------------------------------------------------

const std::string tool_name = "GBWT merging";

void printUsage(int exit_code = EXIT_SUCCESS);

enum MergingAlgorithm { ma_insert, ma_fast, ma_parallel };

std::string algorithmName(MergingAlgorithm algorithm);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 5) { printUsage(); }
  Verbosity::set(Verbosity::FULL);

  size_type batch_size = DynamicGBWT::MERGE_BATCH_SIZE, sample_interval = DynamicGBWT::SAMPLE_INTERVAL;
  MergingAlgorithm algorithm = ma_insert;
  MergeParameters parameters;
  std::string output;
  bool sdsl_format = false;
  int c = 0;
  while((c = getopt(argc, argv, "b:C:fiJ:M:o:OpP:s:S:t:T:")) != -1)
  {
    switch(c)
    {
    case 'b':
      batch_size = std::stoul(optarg); break;
    case 'C':
      parameters.setChunkSize(std::stoul(optarg)); break;
    case 'f':
      algorithm = ma_fast; break;
    case 'i':
      algorithm = ma_insert; break;
    case 'J':
      parameters.setMergeJobs(std::stoul(optarg)); break;
    case 'M':
      parameters.setMergeBuffers(std::stoul(optarg)); break;
    case 'o':
      output = optarg; break;
    case 'O':
      sdsl_format = true; break;
    case 'p':
      algorithm = ma_parallel; break;
    case 'P':
      parameters.setPosBufferSize(std::stoul(optarg)); break;
    case 's':
      sample_interval = std::stoul(optarg); break;
    case 'S':
      omp_set_num_threads(std::stoul(optarg)); break;
    case 't':
      TempFile::setDirectory(optarg); break;
    case 'T':
      parameters.setThreadBufferSize(std::stoul(optarg)); break;
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }

  size_type input_files = argc - optind;
  size_type total_inserted = 0;
  if(input_files <= 1 || output.empty()) { printUsage(EXIT_FAILURE); }

  Version::print(std::cout, tool_name);

  printHeader("Algorithm"); std::cout << algorithmName(algorithm) << std::endl;
  printHeader("Input files"); std::cout << input_files << std::endl;
  printHeader("Output name"); std::cout << output << (sdsl_format ? " (SDSL format)" : " (simple-sds format)") << std::endl;
  if(algorithm == ma_insert)
  {
    printHeader("Batch size"); std::cout << batch_size << std::endl;
    printHeader("Sample interval"); std::cout << sample_interval << std::endl;
  }
  else if(algorithm == ma_parallel)
  {
    printHeader("Temp directory"); std::cout << TempFile::temp_dir << std::endl;
    printHeader("Search threads"); std::cout << omp_get_max_threads() << std::endl;
    printHeader("Pos buffers"); std::cout << parameters.pos_buffer_size << " MB" << std::endl;
    printHeader("Thread buffers"); std::cout << parameters.thread_buffer_size << " MB" << std::endl;
    printHeader("Merge buffers"); std::cout << parameters.merge_buffers << std::endl;
    printHeader("Chunk size"); std::cout << parameters.chunk_size << std::endl;
    printHeader("Merge jobs"); std::cout << parameters.merge_jobs << std::endl;
  }
  std::cout << std::endl;

  double start = readTimer();

  if(algorithm == ma_fast)
  {
    std::vector<GBWT> indexes(argc - optind);
    for(int i = optind; i < argc; i++)
    {
      std::string input_name = argv[i];
      sdsl::simple_sds::load_from(indexes[i - optind], input_name + GBWT::EXTENSION);
      printStatistics(indexes[i - optind], input_name);
      total_inserted += indexes[i - optind].size();
    }
    GBWT merged(indexes);
    if(sdsl_format)
    {
      if(!sdsl::store_to_file(merged, output + GBWT::EXTENSION))
      {
        std::cerr << "merge_gbwt: Cannot write the index to " << (output + GBWT::EXTENSION) << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    else { sdsl::simple_sds::serialize_to(merged, output + GBWT::EXTENSION); }
    printStatistics(merged, output);
  }
  else
  {
    DynamicGBWT index;
    {
      std::string input_name = argv[optind];
      sdsl::simple_sds::load_from(index, input_name + DynamicGBWT::EXTENSION);
      printStatistics(index, input_name);
      optind++;
    }
    while(optind < argc)
    {
      std::string input_name = argv[optind];
      if(algorithm == ma_insert)
      {
        GBWT next;
        sdsl::simple_sds::load_from(next, input_name + GBWT::EXTENSION);
        printStatistics(next, input_name);
        index.merge(next, batch_size, sample_interval);
        total_inserted += next.size();
      }
      else if(algorithm == ma_parallel)
      {
        DynamicGBWT next;
        sdsl::simple_sds::load_from(next, input_name + DynamicGBWT::EXTENSION);
        printStatistics(next, input_name);
        index.merge(next, parameters);
        total_inserted += next.size();
      }
      optind++;
    }
    if(sdsl_format)
    {
      if(!sdsl::store_to_file(index, output + DynamicGBWT::EXTENSION))
      {
        std::cerr << "merge_gbwt: Cannot write the index to " << (output + DynamicGBWT::EXTENSION) << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    else { sdsl::simple_sds::serialize_to(index, output + DynamicGBWT::EXTENSION); }
    printStatistics(index, output);
  }

  double seconds = readTimer() - start;

  std::cout << "Inserted " << total_inserted << " nodes in " << seconds << " seconds ("
            << (total_inserted / seconds) << " nodes/second)" << std::endl;
  std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: merge_gbwt [options] -o output input1 input2 [input3 ...]" << std::endl;
  std::cerr << std::endl;
  std::cerr << "General options:" << std::endl;
  std::cerr << "  -o X  Use X as the base name for output (required)" << std::endl;
  std::cerr << "  -O    Output SDSL format instead of simple-sds format" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Algorithm choice:" << std::endl;
  std::cerr << "  -f    Fast algorithm (node ids must not overlap)" << std::endl;
  std::cerr << "  -i    Insertion algorithm (default)" << std::endl;
  std::cerr << "  -p    Parallel algorithm" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Insertion algorithm (-i):" << std::endl;
  std::cerr << "  -b N  Insert in batches of N sequences (default: " << DynamicGBWT::MERGE_BATCH_SIZE << ")" << std::endl;
  std::cerr << "  -s N  Sample sequence ids at one out of N positions (default: " << DynamicGBWT::SAMPLE_INTERVAL << "; use 0 for no samples)" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Parallel algorithm (-p):" << std::endl;
  std::cerr << "  -C N  Parallelize search in chunks of N sequences (default: " << MergeParameters::CHUNK_SIZE << ")" << std::endl;
  std::cerr << "  -J N  Run N parallel merge jobs (default: " << MergeParameters::MERGE_JOBS << ")" << std::endl;
  std::cerr << "  -M N  Use N merge buffers (default: " << MergeParameters::MERGE_BUFFERS << ")" << std::endl;
  std::cerr << "  -P N  Use N-megabyte position buffers (default: " << MergeParameters::POS_BUFFER_SIZE << ")" << std::endl;
  std::cerr << "  -S N  Use N search threads (default: " << omp_get_max_threads() << ")" << std::endl;
  std::cerr << "  -t X  Use directory X for temporary files (default: " << TempFile::DEFAULT_TEMP_DIR << ")" << std::endl;
  std::cerr << "  -T N  Use N-megabyte thread buffers (default: " << MergeParameters::THREAD_BUFFER_SIZE << ")" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Use base names for the inputs and the output." << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

std::string algorithmName(MergingAlgorithm algorithm)
{
  switch(algorithm)
  {
  case ma_insert:
    return "insert";
  case ma_fast:
    return "fast";
  case ma_parallel:
    return "parallel";
  default:
    return "unknown";
  }
}

//------------------------------------------------------------------------------
