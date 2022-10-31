/*
  Copyright (c) 2017, 2018, 2019, 2020, 2021 Jouni Siren
  Copyright (c) 2015, 2016, 2017 Genome Research Ltd.

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

#include <gbwt/utils.h>

#include <cstdio>
#include <cstdlib>
#include <mutex>
#include <set>
#include <string>

#include <sys/resource.h>
#include <unistd.h>

namespace gbwt
{

//------------------------------------------------------------------------------

// Numerical class constants.

constexpr size_type Verbosity::SILENT;
constexpr size_type Verbosity::BASIC;
constexpr size_type Verbosity::EXTENDED;
constexpr size_type Verbosity::DEFAULT;
constexpr size_type Verbosity::FULL;

constexpr size_type Version::MAJOR_VERSION;
constexpr size_type Version::MINOR_VERSION;
constexpr size_type Version::PATCH_VERSION;
constexpr size_type Version::GBWT_VERSION;
constexpr size_type Version::METADATA_VERSION;
constexpr size_type Version::VARIANT_VERSION;
constexpr size_type Version::R_INDEX_VERSION;

//------------------------------------------------------------------------------

// Other class variables.

const std::string Version::SOURCE_KEY = "source";
const std::string Version::SOURCE_VALUE = "jltsiren/gbwt";

//------------------------------------------------------------------------------

std::vector<range_type>
Range::partition(range_type range, size_type blocks)
{
  if(empty(range)) { return std::vector<range_type>(); }
  blocks = bound(blocks, 1, length(range));

  std::vector<range_type> result(blocks);
  for(size_type block = 0, start = range.first; block < blocks; block++)
  {
    result[block].first = start;
    if(start <= range.second)
    {
      start += std::max((size_type)1, (range.second + 1 - start) / (blocks - block));
    }
    result[block].second = start - 1;
  }

  return result;
}

//------------------------------------------------------------------------------

size_type Verbosity::level = Verbosity::DEFAULT;

void
Verbosity::set(size_type new_level)
{
  level = Range::bound(new_level, SILENT, FULL);
}

std::string
Verbosity::levelName()
{
  switch(level)
  {
    case SILENT:
      return "silent"; break;
    case BASIC:
      return "basic"; break;
    case EXTENDED:
      return "extended"; break;
    case FULL:
      return "full"; break;
  }
  return "unknown";
}

//------------------------------------------------------------------------------

std::string
Version::str(bool verbose)
{
  std::ostringstream ss;
  if(verbose) { ss << "GBWT version "; }
  else { ss << "v"; }
  ss << MAJOR_VERSION << "." << MINOR_VERSION << "." << PATCH_VERSION;
  if(verbose) { ss << " (file format version " << GBWT_VERSION << ")"; }
  return ss.str();
}

void
Version::print(std::ostream& out, const std::string& tool_name, bool verbose, size_type new_lines)
{
  out << tool_name;
  if(verbose) { out << std::endl; }
  else { out << " "; }
  out << str(verbose);
  for(size_type i = 0; i < new_lines; i++) { out << std::endl; }
}

//------------------------------------------------------------------------------

std::ostream&
printHeader(const std::string& header, std::ostream& out)
{
  size_type indent = DEFAULT_INDENT;
  std::string padding;
  if(header.length() + 1 < indent) { padding = std::string(indent - 1 - header.length(), ' '); }
  out << header << ":" << padding;
  return out;
}

void
printTime(const std::string& header, size_type queries, double seconds, std::ostream& out)
{
  printHeader(header, out);
  out << queries << " queries in " << seconds << " seconds ("
      << inMicroseconds(seconds / queries) << " µs/query)" << std::endl;
}

void
printTimeLength(const std::string& header, size_type queries, size_type total_length, double seconds, std::ostream& out)
{
  printHeader(header, out);
  out << queries << " queries of total length " << total_length << " in " << seconds << " seconds ("
      << inMicroseconds(seconds / total_length) << " µs/character)" << std::endl;
}

//------------------------------------------------------------------------------

double
readTimer()
{
  return omp_get_wtime();
}

size_type
memoryUsage()
{
  rusage usage;
  getrusage(RUSAGE_SELF, &usage);
#if defined(__APPLE__) && defined(__MACH__)
  return usage.ru_maxrss;
#else
  return KILOBYTE * usage.ru_maxrss;
#endif
}

//------------------------------------------------------------------------------

namespace TempFile
{
  const std::string DEFAULT_TEMP_DIR = ".";
  std::string temp_dir = DEFAULT_TEMP_DIR;
  std::mutex tempfile_lock;

  // By storing the filenames in a static object, we can delete the remaining
  // temporary files when std::exit() is called.
  struct Handler
  {
    size_type counter;
    std::set<std::string> filenames;

    Handler() :
      counter(0)
    {
    }

    ~Handler()
    {
      std::lock_guard<std::mutex> lock(tempfile_lock);
      for(auto& filename : this->filenames)
      {
        std::remove(filename.c_str());
      }
      this->filenames.clear();
    }
  } handler;

  void
  setDirectory(const std::string& directory)
  {
    std::lock_guard<std::mutex> lock(tempfile_lock);
    if(directory.empty()) { temp_dir = DEFAULT_TEMP_DIR; }
    else if(directory[directory.length() - 1] != '/') { temp_dir = directory; }
    else { temp_dir = directory.substr(0, directory.length() - 1); }
  }

  std::string
  getName(const std::string& name_part)
  {
    char hostname[32];
    gethostname(hostname, 32); hostname[31] = 0;

    std::string filename;
    {
      std::lock_guard<std::mutex> lock(tempfile_lock);
      filename = temp_dir + '/' + name_part + '_'
        + std::string(hostname) + '_'
        + std::to_string(sdsl::util::pid()) + '_'
        + std::to_string(handler.counter);
      handler.filenames.insert(filename);
      handler.counter++;
    }

    return filename;
  }

  void
  remove(std::string& filename)
  {
    if(!(filename.empty()))
    {
      std::remove(filename.c_str());
      {
        std::lock_guard<std::mutex> lock(tempfile_lock);
        handler.filenames.erase(filename);
      }
      filename.clear();
    }
  }

  void
  forget() {
    std::lock_guard<std::mutex> lock(tempfile_lock);
    handler.filenames.clear();
    handler.counter = 0;
  }
} // namespace TempFile

size_type
readRows(const std::string& filename, std::vector<std::string>& rows, bool skip_empty_rows)
{
  std::ifstream input(filename.c_str(), std::ios_base::binary);
  if(!input)
  {
    std::cerr << "readRows(): Cannot open input file " << filename << std::endl;
    return 0;
  }

  size_type chars = 0;
  while(input)
  {
    std::string buf;
    std::getline(input, buf);
    if(skip_empty_rows && buf.empty()) { continue; }
    rows.push_back(buf);
    chars += buf.length();
  }

  input.close();
  return chars;
}

size_type
fileSize(std::ifstream& file)
{
  std::streamoff curr = file.tellg();

  file.seekg(0, std::ios::end);
  std::streamoff size = file.tellg();
  file.seekg(0, std::ios::beg);
  size -= file.tellg();

  file.seekg(curr, std::ios::beg);
  return size;
}

size_type
fileSize(std::ofstream& file)
{
  std::streamoff curr = file.tellp();

  file.seekp(0, std::ios::end);
  std::streamoff size = file.tellp();
  file.seekp(0, std::ios::beg);
  size -= file.tellp();

  file.seekp(curr, std::ios::beg);
  return size;
}

//------------------------------------------------------------------------------

} // namespace gbwt
