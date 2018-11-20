/*
  Copyright (c) 2018 Jouni Siren

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
  bool remove_metadata = false;
  while((c = getopt(argc, argv, "c:h:rs:")) != -1)
  {
    switch(c)
    {
    case 'c':
      set_contigs = true;
      new_contigs = std::stoul(optarg);
      break;
    case 'h':
      set_haplotypes = true;
      new_haplotypes = std::stoul(optarg);
      break;
    case 'r':
      remove_metadata = true;
      break;
    case 's':
      set_samples = true;
      new_samples = std::stoul(optarg);
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
  sdsl::load_from_file(gbwt, index_base + GBWT::EXTENSION);

  bool modified = false;
  if(set_haplotypes)
  {
    gbwt.metadata.setHaplotypes(new_haplotypes);
    modified = true;
  }
  if(set_samples)
  {
    gbwt.metadata.setSamples(new_samples);
    modified = true;
  }
  if(set_contigs)
  {
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
  if(modified) { sdsl::store_to_file(gbwt, index_base + GBWT::EXTENSION); }

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: metadata [options] basename" << std::endl;
  std::cerr << "  -c N  Set the number of contigs to N" << std::endl;
  std::cerr << "  -h N  Set the number of haplotypes to N" << std::endl;
  std::cerr << "  -r    Remove all metadata" << std::endl;
  std::cerr << "  -s N  Set the number of samples to N" << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

//------------------------------------------------------------------------------
