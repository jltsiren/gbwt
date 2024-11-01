/*
  Copyright (c) 2018, 2019, 2020, 2021 Jouni Siren

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
  Verbosity::set(Verbosity::FULL);

  int c = 0;
  bool print_metadata = true, need_metadata = false;
  bool list_samples = false, list_contigs = false, list_paths = false, list_tags = false;
  bool remove_metadata = false;
  bool sdsl_format = false;
  while((c = getopt(argc, argv, "scptrO")) != -1)
  {
    switch(c)
    {
    case 's':
      list_samples = true; print_metadata = false; need_metadata = true;
      break;
    case 'c':
      list_contigs = true; print_metadata = false; need_metadata = true;
      break;
    case 'p':
      list_paths = true; print_metadata = false; need_metadata = true;
      break;
    case 't':
      list_tags = true; print_metadata = false;
      break;
    case 'r':
      remove_metadata = true;
      break;
    case 'O':
      sdsl_format = true;
      break;
    case '?':
      std::exit(EXIT_FAILURE);
    default:
      std::exit(EXIT_FAILURE);
    }
  }

  if(optind >= argc) { printUsage(EXIT_FAILURE); }
  std::string index_base = argv[optind];

  GBWT index;
  sdsl::simple_sds::load_from(index, index_base + GBWT::EXTENSION);
  if(!(index.hasMetadata()) && need_metadata)
  {
    std::cerr << "metadata_tool: No metadata in the GBWT index" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  bool modified = false;
  if(print_metadata)
  {
    std::cout << index.metadata << std::endl;
  }
  if(list_samples)
  {
    if(!(index.metadata.hasSampleNames()))
    {
      std::cerr << "metadata_tool: No sample names in the GBWT index" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    for(size_type i = 0; i < index.metadata.samples(); i++)
    {
      std::cout << index.metadata.sample(i) << std::endl;
    }
  }
  if(list_contigs)
  {
    if(!(index.metadata.hasContigNames()))
    {
      std::cerr << "metadata_tool: No contig names in the GBWT index" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    for(size_type i = 0; i < index.metadata.contigs(); i++)
    {
      std::cout << index.metadata.contig(i) << std::endl;
    }
  }
  if(list_paths)
  {
    if(!(index.metadata.hasPathNames()))
    {
      std::cerr << "metadata_tool: No path names in the GBWT index" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    for(size_type i = 0; i < index.metadata.paths(); i++)
    {
      PathName path = index.metadata.path(i);
      std::cout << "(";
      if(index.metadata.hasSampleNames()) { std::cout << index.metadata.sample(path.sample); }
      else { std::cout << path.sample; }
      std::cout << ", ";
      if(index.metadata.hasContigNames()) { std::cout << index.metadata.contig(path.contig); }
      else { std::cout << path.contig; }
      std::cout << ", " << path.phase << ", " << path.count << ")" << std::endl;
    }
  }
  if(list_tags)
  {
    // TODO: The `Tags` structure could support iteration directly.
    for(auto iter = index.tags.tags.begin(); iter != index.tags.tags.end(); ++iter)
    {
      std::cout << iter->first << " = " << iter->second << std::endl;
    }
  }
  if(remove_metadata)
  {
    index.clearMetadata();
    modified = true;
  }

  if(modified)
  {
    if(sdsl_format)
    {
      if(!sdsl::store_to_file(index, index_base + GBWT::EXTENSION))
      {
        std::cerr << "metadata: Cannot write the index to " << (index_base + GBWT::EXTENSION) << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    else { sdsl::simple_sds::serialize_to(index, index_base + GBWT::EXTENSION); }
  }

  return 0;
}

//------------------------------------------------------------------------------

void
printUsage(int exit_code)
{
  Version::print(std::cerr, tool_name);

  std::cerr << "Usage: metadata [options] basename" << std::endl;
  std::cerr << "  -s    List sample names" << std::endl;
  std::cerr << "  -c    List contig names" << std::endl;
  std::cerr << "  -p    List path names" << std::endl;
  std::cerr << "  -t    List tags" << std::endl;
  std::cerr << "  -r    Remove all metadata" << std::endl;
  std::cerr << "  -O    Output SDSL format instead of simple-sds format" << std::endl;
  std::cerr << std::endl;

  std::exit(exit_code);
}

//------------------------------------------------------------------------------
