#ifndef GBWT_UTILS_H
#define GBWT_UTILS_H

#include <gbwt/config.h>

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <concepts>
#include <limits>
#include <string_view>
#include <type_traits>
#include <vector>

#include <sdsl/int_vector.hpp>
#include <sdsl/sd_vector.hpp>
#include <sdsl/simple_sds.hpp>

// When we build without OpenMP, use these stubs instead of the real runtime
// calls.
#ifdef GBWT_USE_OPENMP
#include <omp.h>
#else
#include <chrono>
inline int omp_get_max_threads() { return 1; }
inline int omp_get_thread_num() { return 0; }
inline void omp_set_num_threads(int) { }
inline double omp_get_wtime()
{
  return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
#endif

#include <zstd.h>

// Boost is only a dependency when StringArray's shared-memory-backed storage
// is enabled.
#if defined(GBWT_ENABLE_SHARED_MEMORY)
#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/shared_memory_object.hpp>
#include <boost/interprocess/sync/named_mutex.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#endif

namespace gbwt
{

#if defined(GBWT_ENABLE_SHARED_MEMORY)
typedef boost::interprocess::allocator<char, boost::interprocess::managed_shared_memory::segment_manager> SharedMemCharAllocatorType;

// Concept to determine if an allocator is the Boost shared-memory one or not.
template<typename Allocator>
concept SharedMemory = std::same_as<Allocator, SharedMemCharAllocatorType>;
#else
// Nothing is in shared memory when we build without it.
template<typename Allocator>
concept SharedMemory = false;
#endif

}

// Parallel sorting is only available with libstdc++ parallel mode.
#if defined(__GLIBCXX__) && defined(GBWT_USE_OPENMP)
#include <parallel/algorithm>
#endif

namespace gbwt
{

/*
  utils.h: Common utility methods.
*/

//------------------------------------------------------------------------------

/*
  We can save a lot of memory during construction by using 32-bit integers instead of
  64-bit integers in the dynamic records and internal buffers. This limits the number
  of nodes, the number of paths, the length of the paths, and the number of
  occurrences of each node to less than 2^32.
*/

#define GBWT_SAVE_MEMORY

//------------------------------------------------------------------------------

typedef std::uint64_t size_type;
typedef std::uint32_t short_type;
typedef std::uint8_t  byte_type;

typedef size_type node_type;  // Node identifier.
typedef size_type comp_type;  // Record identifier, compacted node identifier.
typedef size_type rank_type;  // Rank of incoming / outgoing edge.

#ifdef GBWT_SAVE_MEMORY
typedef std::pair<short_type, short_type> edge_type;
typedef std::pair<short_type, short_type> run_type;
typedef std::pair<short_type, short_type> sample_type;  // (i, DA[i]) within a record
#else
typedef std::pair<node_type, size_type>   edge_type;
typedef std::pair<rank_type, size_type>   run_type;
typedef std::pair<size_type, size_type>   sample_type;  // (i, DA[i]) within a record
#endif

//------------------------------------------------------------------------------

constexpr size_type BYTE_BITS    = 8;
constexpr size_type S_WORD_BITS  = 32;
constexpr size_type WORD_BITS    = 64;

constexpr size_type LOW_MASK     = 0xFFFFFFFFULL;

constexpr size_type KILOBYTE     = 1024;
constexpr size_type MEGABYTE     = KILOBYTE * KILOBYTE;
constexpr size_type GIGABYTE     = KILOBYTE * MEGABYTE;

constexpr double KILOBYTE_DOUBLE = 1024.0;
constexpr double MILLION_DOUBLE  = 1000000.0;
constexpr double MEGABYTE_DOUBLE = KILOBYTE_DOUBLE * KILOBYTE_DOUBLE;
constexpr double GIGABYTE_DOUBLE = KILOBYTE_DOUBLE * MEGABYTE_DOUBLE;

constexpr size_type MILLION      = 1000000;
constexpr size_type BILLION      = 1000 * MILLION;

constexpr node_type ENDMARKER    = 0;

constexpr node_type invalid_node() { return std::numeric_limits<node_type>::max(); }
constexpr comp_type invalid_comp() { return std::numeric_limits<comp_type>::max(); }
constexpr size_type invalid_sequence() { return std::numeric_limits<size_type>::max(); }
constexpr size_type invalid_offset() { return std::numeric_limits<size_type>::max(); }

inline constexpr edge_type invalid_edge() { return edge_type(invalid_node(), invalid_offset()); }
inline constexpr sample_type invalid_sample() { return sample_type(invalid_offset(), invalid_sequence()); }

//------------------------------------------------------------------------------

typedef sdsl::int_vector<0>        text_type;
typedef sdsl::int_vector_buffer<0> text_buffer_type;

#ifdef GBWT_SAVE_MEMORY
typedef std::vector<short_type>    vector_type;
#else
typedef std::vector<node_type>     vector_type;
#endif

//------------------------------------------------------------------------------

/*
  range_type stores a closed range [first, second]. Empty ranges are indicated by
  first > second. The emptiness check uses +1 to handle the common special case
  [0, -1].
*/

typedef std::pair<size_type, size_type> range_type;

struct Range
{
  static size_type length(range_type range)
  {
    return range.second + 1 - range.first;
  }

  static bool empty(range_type range)
  {
    return (range.first + 1 > range.second + 1);
  }

  static bool empty(size_type sp, size_type ep)
  {
    return (sp + 1 > ep + 1);
  }

  static size_type bound(size_type value, range_type bounds)
  {
    return bound(value, bounds.first, bounds.second);
  }

  static size_type bound(size_type value, size_type low, size_type high)
  {
    return std::max(std::min(value, high), low);
  }

  static constexpr range_type empty_range()
  {
    return range_type(1, 0);
  }

  /*
    Partition the range approximately evenly between the blocks. The actual number of
    blocks will not be greater than the length of the range.
  */
  static std::vector<range_type> partition(range_type range, size_type blocks);
};

template<class A, class B>
std::ostream& operator<<(std::ostream& out, const std::pair<A, B>& data)
{
  return out << "(" << data.first << ", " << data.second << ")";
}

template<class A>
std::ostream& operator<<(std::ostream& out, const std::vector<A>& data)
{
  out << "{ ";
  for(const A& element : data) { out << element << " "; }
  out << "}";
  return out;
}

//------------------------------------------------------------------------------

/*
  Global verbosity setting for index construction. Used in conditions of type
  if(Verbosity::level >= Verbosity::THRESHOLD). While the level can be set directly,
  Verbosity::set() does a few sanity checks.

  SILENT    no status information
  BASIC     basic statistics on the input and the final index
  EXTENDED  intermediate statistics for each batch
  FULL      further details of each batch
*/

struct Verbosity
{
  static size_type level;

  static void set(size_type new_level);
  static std::string levelName();

  constexpr static size_type SILENT   = 0;
  constexpr static size_type BASIC    = 1;
  constexpr static size_type EXTENDED = 2;
  constexpr static size_type FULL     = 3;

  // The default is now SILENT, as it's more convenient for users.
  constexpr static size_type DEFAULT = SILENT;
};

//------------------------------------------------------------------------------

struct Version
{
  static std::string str(bool verbose = false);
  static void print(std::ostream& out, const std::string& tool_name, bool verbose = false, size_type new_lines = 2);

  constexpr static size_type MAJOR_VERSION    = 1;
  constexpr static size_type MINOR_VERSION    = 5;
  constexpr static size_type PATCH_VERSION    = 0;

  constexpr static size_type GBWT_VERSION     = 6;
  constexpr static size_type METADATA_VERSION = 2;
  constexpr static size_type VARIANT_VERSION  = 1;
  constexpr static size_type R_INDEX_VERSION  = 1;

  const static std::string SOURCE_KEY; // source
  const static std::string SOURCE_VALUE; // jltsiren/gbwt
  const static std::string SOURCE_GBWT_RS; // jltsiren/gbwt-rs
};

//------------------------------------------------------------------------------

// NOTE: This is deprecated. Use `sdsl::bits::length` instead.
template<class IntegerType>
size_type
bit_length(IntegerType val)
{
  return sdsl::bits::length(val);
}

//------------------------------------------------------------------------------

/*
  Thomas Wang's integer hash function. In many implementations, std::hash
  is identity function for integers, which leads to performance issues.
*/

inline size_type
wang_hash_64(size_type key)
{
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return key;
}

//------------------------------------------------------------------------------

inline double
inMegabytes(size_type bytes)
{
  return bytes / MEGABYTE_DOUBLE;
}

inline double
inGigabytes(size_type bytes)
{
  return bytes / GIGABYTE_DOUBLE;
}

inline double
inBPC(size_type bytes, size_type size)
{
  return (8.0 * bytes) / size;
}

inline double
inMicroseconds(double seconds)
{
  return seconds * MILLION_DOUBLE;
}

constexpr size_type DEFAULT_INDENT = 18;

std::ostream& printHeader(const std::string& header, std::ostream& out = std::cout);
void printTime(const std::string& header, size_type queries, double seconds, std::ostream& out = std::cout);
void printTimeLength(const std::string& header, size_type queries, size_type total_length, double seconds, std::ostream& out = std::cout);

//------------------------------------------------------------------------------

double readTimer();       // Seconds from an arbitrary time point.
size_type memoryUsage();  // Peak memory usage in bytes.

//------------------------------------------------------------------------------

/*
  Temporary file names have the pattern "prefix_hostname_pid_counter", where
  - prefix is given as an argument to getName();
  - hostname is the name of the host;
  - pid is the process id; and
  - counter is a running counter starting from 0.

  The generated names are stored until the file is deleted with remove(). All
  remaining temporary files are deleted when the program exits (normally or
  with std::exit()).

  TempFile is thread-safe.
*/

namespace TempFile
{
  extern const std::string DEFAULT_TEMP_DIR;
  extern std::string temp_dir;

  void setDirectory(const std::string& directory);
  std::string getName(const std::string& name_part);
  void remove(std::string& filename);  // Also clears the filename.
  // Forget about current temporary files so that they aren't deleted.
  void forget();
}

// Returns the total length of the rows, excluding line ends.
// If 'rows' is nonempty, appends to it.
size_type readRows(const std::string& filename, std::vector<std::string>& rows, bool skip_empty_rows);

size_type fileSize(std::ifstream& file);
size_type fileSize(std::ofstream& file);

//------------------------------------------------------------------------------

/*
  parallelQuickSort() uses less working space than parallelMergeSort(). Calling omp_set_nested(1)
  improves the speed of parallelQuickSort().
*/

template<class Iterator, class Comparator>
void
parallelQuickSort(Iterator first, Iterator last, const Comparator& comp)
{
#if defined(__GLIBCXX__) && defined(GBWT_USE_OPENMP)
  int nested = omp_get_nested();
  omp_set_nested(1);
  __gnu_parallel::sort(first, last, comp, __gnu_parallel::balanced_quicksort_tag());
  omp_set_nested(nested);
#else
  std::sort(first, last, comp);
#endif
}

template<class Iterator>
void
parallelQuickSort(Iterator first, Iterator last)
{
#if defined(__GLIBCXX__) && defined(GBWT_USE_OPENMP)
  int nested = omp_get_nested();
  omp_set_nested(1);
  __gnu_parallel::sort(first, last, __gnu_parallel::balanced_quicksort_tag());
  omp_set_nested(nested);
#else
  std::sort(first, last);
#endif
}

template<class Iterator, class Comparator>
void
parallelMergeSort(Iterator first, Iterator last, const Comparator& comp)
{
#if defined(__GLIBCXX__) && defined(GBWT_USE_OPENMP)
  __gnu_parallel::sort(first, last, comp, __gnu_parallel::multiway_mergesort_tag());
#else
  std::sort(first, last, comp);
#endif
}

template<class Iterator>
void
parallelMergeSort(Iterator first, Iterator last)
{
#if defined(__GLIBCXX__) && defined(GBWT_USE_OPENMP)
  __gnu_parallel::sort(first, last, __gnu_parallel::multiway_mergesort_tag());
#else
  std::sort(first, last);
#endif
}

template<class Iterator, class Comparator>
void
sequentialSort(Iterator first, Iterator last, const Comparator& comp)
{
#if defined(__GLIBCXX__) && defined(GBWT_USE_OPENMP)
  __gnu_parallel::sort(first, last, comp, __gnu_parallel::sequential_tag());
#else
  std::sort(first, last, comp);
#endif
}

template<class Iterator>
void
sequentialSort(Iterator first, Iterator last)
{
#if defined(__GLIBCXX__) && defined(GBWT_USE_OPENMP)
  __gnu_parallel::sort(first, last, __gnu_parallel::sequential_tag());
#else
  std::sort(first, last);
#endif
}

template<class Element>
void
removeDuplicates(std::vector<Element>& vec, bool parallel)
{
  if(parallel) { parallelQuickSort(vec.begin(), vec.end()); }
  else         { sequentialSort(vec.begin(), vec.end()); }
  vec.resize(std::unique(vec.begin(), vec.end()) - vec.begin());
}

//------------------------------------------------------------------------------

// Zstandard wrapper that stores the compressed output in a vector.
// Throws `std::runtime_error` on failure.
class ZstdCompressor
{
public:
  constexpr static int DEFAULT_COMPRESSION_LEVEL = 3;

  explicit ZstdCompressor(int compression_level = DEFAULT_COMPRESSION_LEVEL);
  ~ZstdCompressor();
  ZstdCompressor(const ZstdCompressor&) = delete;
  ZstdCompressor& operator=(const ZstdCompressor&) = delete;

  // Compresses the given data, buffering it internally.
  void compress(std::string_view data);

  // Compresses the given data directly.
  // Flushes the internal input buffer if necessary.
  void compressDirect(std::string_view data);

  // Finalizes the compression.
  // The compressor cannot be used after this call.
  void finish();

  // Returns the compressed data.
  // Valid after `finish()` has been called.
  const std::vector<char>& outputData() const { return this->output; }

private:
  // Compresses and clears the internal input buffer.
  void flushInput();

  // Compresses the given input buffer.
  void compress(ZSTD_inBuffer& buffer);

  // Flushes the internal output buffer to the output vector.
  void flushOutput();

  ZSTD_CCtx* context;

  std::vector<char> input_buffer;
  size_t input_buffer_capacity;
  std::vector<char> output_buffer;
  ZSTD_outBuffer out_buffer;

  std::vector<char> output;
};

// Zstandard decompression wrapper that decompresses data from an internal vector.
// Throws `sdsl::simple_sds::InvalidData` on failure.
class ZstdDecompressor
{
public:
  explicit ZstdDecompressor(std::vector<char>&& input);
  ~ZstdDecompressor();
  ZstdDecompressor(const ZstdDecompressor&) = delete;
  ZstdDecompressor& operator=(const ZstdDecompressor&) = delete;

  // Decompresses the given number of bytes and appends them to the output
  // vector, whatever byte type and allocator it uses.
  template<typename Element, typename Allocator>
  void decompress(size_t bytes, std::vector<Element, Allocator>& output);

  // Returns `true` if all input data has been consumed.
  bool finished();

private:
  // Decompress more data into the output buffer.
  // Assumes that the output buffer has been consumed.
  void fillOutputBuffer();

  ZSTD_DCtx* context;

  std::vector<char> input;
  ZSTD_inBuffer in_buffer;
  std::vector<char> output_buffer;
  ZSTD_outBuffer out_buffer;
  size_t cursor; // Next unread byte in `out_buffer`.
};

template<typename Element, typename Allocator>
void
ZstdDecompressor::decompress(size_t bytes, std::vector<Element, Allocator>& output)
{
  static_assert(sizeof(Element) == 1, "ZstdDecompressor decompresses into a vector of bytes");

  size_t decompressed = 0;
  while(decompressed < bytes)
  {
    if(this->cursor < this->out_buffer.pos)
    {
      size_t to_copy = std::min(bytes - decompressed, this->out_buffer.pos - this->cursor);
      const Element* start_ptr = static_cast<const Element*>(this->out_buffer.dst) + this->cursor;
      output.insert(output.end(), start_ptr, start_ptr + to_copy);
      this->cursor += to_copy;
      decompressed += to_copy;
    }
    else if(this->in_buffer.pos < this->in_buffer.size)
    {
      this->fillOutputBuffer();
    }
    else
    {
      std::string msg = "ZstdDecompressor: Unexpected end of input data";
      throw sdsl::simple_sds::InvalidData(msg);
    }
  }
}

//------------------------------------------------------------------------------

} // namespace gbwt

#endif // GBWT_UTILS_H
