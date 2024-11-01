/*
  Copyright (c) 2018, 2019, 2021 Jouni Siren

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

#include <gbwt/dynamic_gbwt.h>

using namespace gbwt;

//------------------------------------------------------------------------------

const std::string tool_name = "GBWT sequence removal";

void printUsage(int exit_code = EXIT_SUCCESS);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 3) { printUsage(); }
  Verbosity::set(Verbosity::FULL);

  // Parse command line options.
  int c = 0;
  std::string output;
  size_type chunk_size = DynamicGBWT::REMOVE_CHUNK_SIZE;
  bool range = false, sample = false, contig = false;
  bool sdsl_format = false;
  while((c = getopt(argc, argv, "c:o:OrSC")) != -1)
  {
    switch(c)
    {
    case 'c':
      chunk_size = std::max(1ul, std::stoul(optarg));
      break;
    case 'o':
      output = optarg; break;
    case 'O':
      sdsl_format = true; break;
    case 'r':
      range = true; break;
    case 'S':
      sample = true; break;
    case 'C':
      contig = true; break;
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }

  // Check command line options.
  if(optind + 1 >= argc) { printUsage(EXIT_FAILURE); }
  if((range & sample) || (range & contig) || (sample & contig))
  {
    std::cerr << "remove_seq: Options -r, -S, and -C are mutually exclusive" << std::endl;
    std::exit(EXIT_FAILURE);
  }
  std::string base_name = argv[optind]; optind++;
  std::string key;
  if(output.empty()) { output = base_name; }
  std::vector<size_type> seq_ids;
  if(range)
  {
    if(argc != optind + 2) { printUsage(EXIT_FAILURE); }
    size_type start = std::stoul(argv[optind]); optind++;
    size_type stop = std::stoul(argv[optind]); optind++;
    if(stop < start) { printUsage(EXIT_FAILURE); }
    for(size_type seq_id = start; seq_id <= stop; seq_id++)
    {
      seq_ids.push_back(seq_id);
    }
  }
  else if(sample || contig)
  {
    if(argc != optind + 1) { printUsage(EXIT_FAILURE); }
    key = argv[optind];
  }
  else
  {
    while(optind < argc)
    {
      seq_ids.push_back(std::stoul(argv[optind]));
      optind++;
    }
  }

  // Initial output.
  Version::print(std::cout, tool_name);
  printHeader("Input"); std::cout << base_name << std::endl;
  printHeader("Output"); std::cout << output << std::endl;
  if(range)
  {
    printHeader("Range"); std::cout << range_type(seq_ids.front(), seq_ids.back()) << std::endl;
  }
  else if(sample)
  {
    printHeader("Sample"); std::cout << key << std::endl;
  }
  else if(contig)
  {
    printHeader("Contig"); std::cout << key << std::endl;
  }
  else
  {
    printHeader("Sequences"); std::cout << seq_ids.size() << std::endl;
  }
  printHeader("Chunk size"); std::cout << chunk_size << std::endl;
  std::cout << std::endl;

  double start = readTimer();

  // Load index.
  DynamicGBWT index;
  sdsl::simple_sds::load_from(index, base_name + DynamicGBWT::EXTENSION);
  printStatistics(index, base_name);

  // Handle metadata.
  if(sample)
  {
    if(!(index.hasMetadata()) || !(index.metadata.hasSampleNames()) || !(index.metadata.hasPathNames()))
    {
      std::cerr << "remove_seq: Option -S requires sample and path names" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    size_type sample_id = index.metadata.sample(key);
    seq_ids = index.metadata.removeSample(sample_id);
    if(seq_ids.empty())
    {
      std::cerr << "remove_seq: No sequences for sample " << key << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  else if(contig)
  {
    if(!(index.hasMetadata()) || !(index.metadata.hasContigNames()) || !(index.metadata.hasPathNames()))
    {
      std::cerr << "remove_seq: Option -C requires contig and path names" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    size_type contig_id = index.metadata.contig(key);
    seq_ids = index.metadata.removeContig(contig_id);
    if(seq_ids.empty())
    {
      std::cerr << "remove_seq: No sequences for contig " << key << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  else
  {
    if(index.hasMetadata() && index.metadata.hasPathNames())
    {
      std::cerr << "remove_seq: Removing arbitrary sequences would leave the metadata inconsistent" << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // Remove the sequences.
  size_type total_length = index.remove(seq_ids, chunk_size);
  if(total_length > 0)
  {
    if(sdsl_format)
    {
      if(!sdsl::store_to_file(index, output + DynamicGBWT::EXTENSION))
      {
        std::cerr << "remove_seq: Cannot write the index to " << (output + DynamicGBWT::EXTENSION) << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    else { sdsl::simple_sds::serialize_to(index, output + DynamicGBWT::EXTENSION); }
    printStatistics(index, output);
  }

  double seconds = readTimer() - start;

  std::cout << "Removed " << total_length << " nodes in " << seconds << " seconds ("
            << (total_length / seconds) << " nodes/second)" << std::endl;
  std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: remove_seq [options] base_name seq1 [seq2 ...]" << std::endl;
  std::cerr << std::endl;
  std::cerr << "  -c N  Build the RA in chunks of N sequences per thread (default: " << DynamicGBWT::REMOVE_CHUNK_SIZE << ")" << std::endl;
  std::cerr << "  -o X  Use X as the base name for output" << std::endl;
  std::cerr << "  -O    Output SDSL format instead of simple-sds format" << std::endl;
  std::cerr << "  -r    Remove a range of sequences (inclusive; requires 2 sequence ids)" << std::endl;
  std::cerr << "  -S    Remove all sequences for the sample with name seq1 (cannot have seq2)" << std::endl;
  std::cerr << "  -C    Remove all sequences for the contig with name seq1 (cannot have seq2)" << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

//------------------------------------------------------------------------------
