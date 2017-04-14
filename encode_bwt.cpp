/*
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

#include <map>

#include "files.h"
#include "support.h"

using namespace gbwt;

//------------------------------------------------------------------------------

int
main(int argc, char** argv)
{
  if(argc < 2)
  {
    std::cerr << "Usage: encode_bwt base_name" << std::endl;
    std::cerr << std::endl;
    std::cerr << "Builds the suffix array and the BWT." << std::endl;
    std::cerr << std::endl;
    std::exit(EXIT_SUCCESS);
  }

  std::string base_name = argv[1];
  std::string header_name = base_name + HEADER_EXTENSION;
  std::string alphabet_name = base_name + ALPHABET_EXTENSION;
  std::string bwt_name = base_name + BWT_EXTENSION;
  std::string gbwt_name = base_name + GBWT_EXTENSION;
  std::string index_name = base_name + INDEX_EXTENSION;

  std::cout << "Encoding the BWT" << std::endl;
  std::cout << std::endl;

  printHeader("Base name"); std::cout << base_name << std::endl;
  std::cout << std::endl;

  double start = readTimer();

  GBWTHeader header;
  sdsl::load_from_file(header, header_name);
  std::cout << header << std::endl;
  std::cout << std::endl;

  // Open input.
  sdsl::sd_vector<> alphabet;
  sdsl::load_from_file(alphabet, alphabet_name);
  sdsl::sd_vector<>::select_1_type alphabet_select;
  sdsl::util::init_support(alphabet_select, &alphabet);
  text_buffer_type bwt(bwt_name);

  // Prepare output and statistics.
  sdsl::int_vector_buffer<8> gbwt(gbwt_name, std::ios::out);
  std::vector<size_type> node_borders;
  std::map<size_type, size_type> alphabet_distribution; // effective alphabet size -> number of nodes
  std::map<size_type, size_type> run_distribution;  // run count -> number_of_nodes
  std::vector<size_type> counts(header.alphabet_size, 0);

  for(size_type node = 0; node < header.alphabet_size; node++)
  {
    // Determine the effective alphabet.
    size_type start = alphabet_select(node + 1) - node;
    size_type stop = alphabet_select(node + 2) - node - 1;
    std::map<value_type, value_type> effective_alphabet;
    for(size_type i = start; i < stop; i++)
    {
      effective_alphabet[bwt[i]] = 0;
    }
    alphabet_distribution[effective_alphabet.size()]++;
    if(effective_alphabet.size() == 0) { continue; }

    // Write the alphabet.
    node_borders.push_back(gbwt.size());
    ByteCode::write(gbwt, effective_alphabet.size());
    size_type c = 0;
    for(auto iter = effective_alphabet.begin(); iter != effective_alphabet.end(); ++iter)
    {
      iter->second = c; c++;
      ByteCode::write(gbwt, iter->first);
      ByteCode::write(gbwt, counts[iter->first]);
    }

    // Write the runs.
    size_type run_count = 0;
    Run encoder(effective_alphabet.size());
    Run::run_type run(0, 0);
    for(size_type i = start; i < stop; i++)
    {
      if(bwt[i] == run.first) { run.second++; }
      else
      {
        if(run.second > 0)
        {
          counts[run.first] += run.second;
          run.first = effective_alphabet[run.first];
          encoder.write(gbwt, run); run_count++;
        }
        run.first = bwt[i]; run.second = 1;
      }
    }
    if(run.second > 0)
    {
      counts[run.first] += run.second;
      run.first = effective_alphabet[run.first];
      encoder.write(gbwt, run); run_count++;
    }
    run_distribution[run_count]++;
  }
  node_borders.push_back(gbwt.size());
  bwt.close(); gbwt.close();

  // Write the node index.
  sdsl::sd_vector<> node_index(node_borders.begin(), node_borders.end());
  sdsl::store_to_file(node_index, index_name);

  double seconds = readTimer() - start;
  size_type data_size = header.sequences + header.total_length;

  std::cout << "Effective alphabet size distribution:" << std::endl;
  for(auto value : alphabet_distribution)
  {
    std::cout << value << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Run count distribution:" << std::endl;
  for(auto value : run_distribution)
  {
    std::cout << value << std::endl;
  }
  std::cout << std::endl;

  std::cout << "Indexed " << data_size << " nodes in " << seconds << " seconds (" << (data_size / seconds) << " nodes/second)" << std::endl;
  std::cout << "Memory usage " << inGigabytes(memoryUsage()) << " GB" << std::endl;
  std::cout << std::endl;

  return 0;
}

//------------------------------------------------------------------------------
