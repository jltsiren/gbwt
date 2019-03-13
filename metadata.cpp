/*
  Copyright (c) 2018, 2019 Jouni Siren

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

#include <gbwt/gbwt.h>

using namespace gbwt;

//------------------------------------------------------------------------------

const std::string tool_name = "GBWT metadata tool";

void printUsage(int exit_code = EXIT_SUCCESS);

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2) { printUsage(); }

  int c = 0;
  bool set_samples = false, set_haplotypes = false, set_contigs = false;
  size_type new_samples = 0, new_haplotypes = 0, new_contigs = 0;
  bool clear_paths = false, clear_samples = false, clear_contigs = false;
  bool remove_metadata = false;
  while((c = getopt(argc, argv, "s:h:c:PSRr")) != -1)
  {
    switch(c)
    {
    case 's':
      set_samples = true;
      new_samples = std::stoul(optarg);
      break;
    case 'h':
      set_haplotypes = true;
      new_haplotypes = std::stoul(optarg);
      break;
    case 'c':
      set_contigs = true;
      new_contigs = std::stoul(optarg);
      break;
    case 'P':
      clear_paths = true;
      break;
    case 'S':
      clear_samples = true;
      break;
    case 'C':
      clear_contigs = true;
      break;
    case 'r':
      remove_metadata = true;
      break;
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }

  if(optind >= argc) { printUsage(EXIT_FAILURE); }
  std::string index_base = argv[optind];

  Version::print(std::cout, tool_name);
  GBWT gbwt;
  if(!sdsl::load_from_file(gbwt, index_base + GBWT::EXTENSION))
  {
    std::cerr << "metadata: Cannot load the index from " << (index_base + GBWT::EXTENSION) << std::endl;
    std::exit(EXIT_FAILURE);
  }

  bool modified = false;
  if(clear_paths)
  {
    gbwt.metadata.clearPaths();
    modified = true;
  }
  if(clear_samples)
  {
    gbwt.metadata.clearSamples();
    modified = true;
  }
  if(set_samples)
  {
    if(gbwt.metadata.get(Metadata::FLAG_SAMPLE_NAMES))
    {
      std::cerr << "metadata: Changing sample count without changing sample names" << std::endl;
    }
    gbwt.metadata.setSamples(new_samples);
    modified = true;
  }
  if(set_haplotypes)
  {
    gbwt.metadata.setHaplotypes(new_haplotypes);
    modified = true;
  }
  if(clear_contigs)
  {
    gbwt.metadata.clearContigs();
    modified = true;
  }
  if(set_contigs)
  {
    if(gbwt.metadata.get(Metadata::FLAG_SAMPLE_NAMES))
    {
      std::cerr << "metadata: Changing contig count without changing contig names" << std::endl;
    }
    gbwt.metadata.setContigs(new_contigs);
    modified = true;
  }
  if(remove_metadata)
  {
    gbwt.clearMetadata();
    modified = true;
  }
  if(modified && !remove_metadata) { gbwt.addMetadata(); }

  printStatistics(gbwt, index_base);
  if(modified)
  {
    if(!sdsl::store_to_file(gbwt, index_base + GBWT::EXTENSION))
    {
      std::cerr << "metadata: Cannot write the index to " << (index_base + GBWT::EXTENSION) << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: metadata [options] basename" << std::endl;
  std::cerr << "  -s N  Set the number of samples to N" << std::endl;
  std::cerr << "  -h N  Set the number of haplotypes to N" << std::endl;
  std::cerr << "  -c N  Set the number of contigs to N" << std::endl;
  std::cerr << "  -P    Remove path names" << std::endl;
  std::cerr << "  -S    Remove sample names" << std::endl;
  std::cerr << "  -C    Remove contig names" << std::endl;
  std::cerr << "  -r    Remove all metadata" << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

//------------------------------------------------------------------------------
