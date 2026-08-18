#include <gtest/gtest.h>

#include <random>

#include <gbwt/support.h>

using namespace gbwt;

namespace
{

//------------------------------------------------------------------------------

class StringArrayTest : public ::testing::Test
{
public:
  void check_array(const StringArray<>& array, const std::vector<std::string>& truth) const
  {
    ASSERT_EQ(array.size(), truth.size()) << "Incorrect array size";
    ASSERT_EQ(array.empty(), truth.empty()) << "Incorrect emptiness";
    size_t total_length = 0;
    for(const std::string& s : truth) { total_length += s.length(); }
    ASSERT_EQ(array.length(), total_length) << "Incorrect total length of the strings";

    for(size_type i = 0; i < array.size(); i++)
    {
      std::string correct = truth[i];
      ASSERT_EQ(array.str(i), correct) << "Incorrect string " << i;
      EXPECT_EQ(array.length(i), correct.length()) << "Incorrect length for string " << i;
      std::string_view view = array.view(i);
      std::string from_view(view.data(), view.size());
      EXPECT_EQ(from_view, correct) << "Incorrect view of string " << i;
    }
  }

  void try_remove(const std::vector<std::string>& original, size_type i) const
  {
    std::vector<std::string> copy = original;
    if(i < original.size()) { copy.erase(copy.begin() + i); }
    StringArray<> truth(copy);

    StringArray<> removed(original);
    removed.remove(i);
    EXPECT_EQ(removed, truth) << "Remove failed for " << i << " / " << original.size();
  }

  void check_file_size(const StringArray<>& original, std::ifstream& in) const
  {
    size_type expected_size = original.simple_sds_size() * sizeof(sdsl::simple_sds::element_type);
    size_type file_size = fileSize(in);
    ASSERT_EQ(expected_size, file_size) << "Incorrect file size";
  }

  void simple_sds_duplicate_from(StringArray<>& array, const std::string& filename, const std::function<std::string(std::string_view)>& transform) const
  {
    std::ifstream in(filename, std::ios_base::binary);
    ASSERT_TRUE(in) << "Cannot open input file " << filename;
    in.exceptions(std::ios::badbit | std::ios::failbit);
    ASSERT_NO_THROW(array.simple_sds_load_duplicate(in, transform)) << "simple_sds_load_duplicate() failed";
    in.close();
  }

  void zstd_compress_to(const StringArray<>& array, const std::string& filename) const
  {
    std::ofstream out(filename, std::ios_base::binary);
    ASSERT_TRUE(out) << "Cannot open output file " << filename;
    out.exceptions(std::ios::badbit | std::ios::failbit);
    ASSERT_NO_THROW(array.simple_sds_compress(out)) << "simple_sds_compress() failed";
    out.close();
  }

  void zstd_decompress_from(StringArray<>& array, const std::string& filename) const
  {
    std::ifstream in(filename, std::ios_base::binary);
    ASSERT_TRUE(in) << "Cannot open input file " << filename;
    in.exceptions(std::ios::badbit | std::ios::failbit);
    ASSERT_NO_THROW(array.simple_sds_decompress(in)) << "simple_sds_decompress() failed";
    in.close();
  }

  void zstd_compress_even_to(const StringArray<>& array, const std::string& filename) const
  {
    std::ofstream out(filename, std::ios_base::binary);
    ASSERT_TRUE(out) << "Cannot open output file " << filename;
    out.exceptions(std::ios::badbit | std::ios::failbit);
    ASSERT_NO_THROW(array.simple_sds_compress_even(out)) << "simple_sds_compress_even() failed";
    out.close();
  }

  void zstd_decompress_duplicate_from(StringArray<>& array, const std::string& filename, const std::function<std::string(std::string_view)>& transform) const
  {
    std::ifstream in(filename, std::ios_base::binary);
    ASSERT_TRUE(in) << "Cannot open input file " << filename;
    in.exceptions(std::ios::badbit | std::ios::failbit);
    ASSERT_NO_THROW(array.simple_sds_decompress_duplicate(in, transform)) << "simple_sds_decompress_duplicate() failed";
    in.close();
  }
};

TEST_F(StringArrayTest, DefaultEmptyArray)
{
  std::vector<std::string> truth;
  StringArray<> array;
  this->check_array(array, truth);
}

TEST_F(StringArrayTest, FromEmptyArray)
{
  std::vector<std::string> truth;
  StringArray<> array(truth);
  this->check_array(array, truth);
}

TEST_F(StringArrayTest, NonEmptyArray)
{
  std::vector<std::string> truth
  {
    "first", "second", "third", "fourth"
  };
  StringArray<> array(truth);
  this->check_array(array, truth);
}

TEST_F(StringArrayTest, FromMap)
{
  std::map<std::string, std::string> source
  {
    { "A-key", "A-value" },
    { "C-key", "C-value" },
    { "B-key", "B-value" },
  };
  std::vector<std::string> truth
  {
    "A-key", "A-value", "B-key", "B-value", "C-key", "C-value"
  };
  StringArray<> array(source);
  this->check_array(array, truth);
}

TEST_F(StringArrayTest, Choose)
{
  std::vector<std::string> original
  {
    "first", "second", "third", "fourth"
  };
  std::vector<std::string> truth
  {
    "second", "fourth"
  };

  // Choose odd positions (that contain even numbers).
  StringArray<> array(original.size(),
  [&original](size_type i) -> std::string_view
  {
    return std::string_view(original[i]);
  },
  [](size_type i) -> bool
  {
    return ((i & 1) != 0);
  });
  this->check_array(array, truth);
}

TEST_F(StringArrayTest, Remove)
{
  // Try removing the string at each position of the original array.
  std::vector<std::string> original
  {
    "first", "second", "third", "fourth"
  };
  for(size_type i = 0; i <= original.size(); i++)
  {
    this->try_remove(original, i);
  }

  // Try removing the only string.
  std::vector<std::string> one
  {
    "one"
  };
  this->try_remove(one, 0);
}

TEST_F(StringArrayTest, SDSLEmpty)
{
  std::vector<std::string> truth;
  StringArray<> original(truth);

  std::string filename = TempFile::getName("string-array");
  sdsl::store_to_file(original, filename);

  StringArray<> copy; sdsl::load_from_file(copy, filename);
  ASSERT_EQ(copy, original) << "SDSL serialization changed the empty array";

  TempFile::remove(filename);
}

TEST_F(StringArrayTest, SimpleSDSEmpty)
{
  std::vector<std::string> truth;
  StringArray<> original(truth);

  std::string filename = TempFile::getName("string-array");
  sdsl::simple_sds::serialize_to(original, filename);

  StringArray<> copy;
  std::ifstream in(filename, std::ios_base::binary);
  this->check_file_size(original, in);
  copy.simple_sds_load(in);
  in.close();
  ASSERT_EQ(copy, original) << "Simple-SDS serialization changed the empty array";

  TempFile::remove(filename);
}

TEST_F(StringArrayTest, ZstdEmpty)
{
  std::vector<std::string> truth;
  StringArray<> original(truth);

  std::string filename = TempFile::getName("string-array");
  this->zstd_compress_to(original, filename);

  StringArray<> copy;
  this->zstd_decompress_from(copy, filename);
  ASSERT_EQ(copy, original) << "Zstd compression changed the empty array";

  TempFile::remove(filename);
}

std::string
reverse_string(std::string_view s)
{
  std::string str(s);
  std::reverse(str.begin(), str.end());
  return str;
}

StringArray<>
duplicate_array(const std::vector<std::string>& source)
{
  return StringArray<>(2 * source.size(),
  [&](size_type i) -> size_type
  {
    return source[i / 2].length();
  },
  [&](size_type i) -> std::string
  {
    std::string value = source[i / 2];
    if(i & 1) { value = reverse_string(value); }
    return value;
  });
}

TEST_F(StringArrayTest, DuplicateEmpty)
{
  std::vector<std::string> source;
  StringArray<> original(source);
  StringArray<> truth = duplicate_array(source);

  std::string filename = TempFile::getName("string-array");
  sdsl::simple_sds::serialize_to(original, filename);

  StringArray<> copy;
  this->simple_sds_duplicate_from(copy, filename, reverse_string);
  ASSERT_EQ(copy, truth) << "Simple-SDS serialization changed the empty array";

  TempFile::remove(filename);
}

TEST_F(StringArrayTest, ZstdDuplicateEmpty)
{
  std::vector<std::string> source;
  StringArray<> original = duplicate_array(source);

  std::string filename = TempFile::getName("string-array");
  this->zstd_compress_even_to(original, filename);

  StringArray<> copy;
  this->zstd_decompress_duplicate_from(copy, filename, reverse_string);
  ASSERT_EQ(copy, original) << "Zstd compression changed the empty array";

  TempFile::remove(filename);
}

TEST_F(StringArrayTest, ZstdDuplicateEmptyStrings)
{
  std::vector<std::string> source
  {
    "", "", ""
  };
  StringArray<> original = duplicate_array(source);

  std::string filename = TempFile::getName("string-array");
  this->zstd_compress_even_to(original, filename);

  StringArray<> copy;
  this->zstd_decompress_duplicate_from(copy, filename, reverse_string);
  ASSERT_EQ(copy, original) << "Zstd compression changed the array with empty strings";

  TempFile::remove(filename);
}

TEST_F(StringArrayTest, SDSLNonEmpty)
{
  std::vector<std::string> truth
  {
    "first",
    "second",
    "third",
    "fourth"
  };
  StringArray<> original(truth);

  std::string filename = TempFile::getName("string-array");
  sdsl::store_to_file(original, filename);

  StringArray<> copy; sdsl::load_from_file(copy, filename);
  ASSERT_EQ(copy, original) << "SDSL serialization changed the non-empty array";

  TempFile::remove(filename);
}

TEST_F(StringArrayTest, SimpleSDSNonEmpty)
{
  std::vector<std::string> truth
  {
    "first",
    "second",
    "third",
    "fourth"
  };
  StringArray<> original(truth);

  std::string filename = TempFile::getName("string-array");
  sdsl::simple_sds::serialize_to(original, filename);

  StringArray<> copy;
  std::ifstream in(filename, std::ios_base::binary);
  this->check_file_size(original, in);
  copy.simple_sds_load(in);
  in.close();
  ASSERT_EQ(copy, original) << "Simple-SDS serialization changed the non-empty array";

  TempFile::remove(filename);
}

TEST_F(StringArrayTest, ZstdNonEmpty)
{
  std::vector<std::string> truth
  {
    "first",
    "second",
    "third",
    "fourth"
  };
  StringArray<> original(truth);

  std::string filename = TempFile::getName("string-array");
  this->zstd_compress_to(original, filename);

  StringArray<> copy;
  this->zstd_decompress_from(copy, filename);
  ASSERT_EQ(copy, original) << "Zstd compression changed the non-empty array";

  TempFile::remove(filename);
}

TEST_F(StringArrayTest, DuplicateNonEmpty)
{
  std::vector<std::string> source
  {
    "first",
    "second",
    "third",
    "fourth"
  };
  StringArray<> original(source);
  StringArray<> truth = duplicate_array(source);

  std::string filename = TempFile::getName("string-array");
  sdsl::simple_sds::serialize_to(original, filename);

  StringArray<> copy;
  this->simple_sds_duplicate_from(copy, filename, reverse_string);
  ASSERT_EQ(copy, truth) << "Simple-SDS serialization changed the non-empty array";

  TempFile::remove(filename);
}

TEST_F(StringArrayTest, ZstdDuplicateNonEmpty)
{
  std::vector<std::string> source
  {
    "first",
    "second",
    "third",
    "fourth"
  };
  StringArray<> original = duplicate_array(source);

  std::string filename = TempFile::getName("string-array");
  this->zstd_compress_even_to(original, filename);

  StringArray<> copy;
  this->zstd_decompress_duplicate_from(copy, filename, reverse_string);
  ASSERT_EQ(copy, original) << "Zstd compression changed the non-empty array";

  TempFile::remove(filename);
}

TEST_F(StringArrayTest, SimpleSDSWithEmptyStrings)
{
  // Here we test that the compression still works when there is an empty
  // string in the middle and the sd_vector used for the offsets contains
  // duplicate values. There is also an empty string at the end, as this
  // used to fail in the past.
  std::vector<std::string> truth
  {
    "first",
    "second",
    "",
    "fourth",
    ""
  };
  StringArray<> original(truth);

  std::string filename = TempFile::getName("string-array");
  sdsl::simple_sds::serialize_to(original, filename);

  StringArray<> copy;
  std::ifstream in(filename, std::ios_base::binary);
  this->check_file_size(original, in);
  copy.simple_sds_load(in);
  in.close();
  ASSERT_EQ(copy, original) << "Simple-SDS serialization changed the array with empty strings";

  TempFile::remove(filename);
}

TEST_F(StringArrayTest, ZstdWithEmptyStrings)
{
  std::vector<std::string> truth
  {
    "first",
    "second",
    "",
    "fourth",
    ""
  };
  StringArray<> original(truth);

  std::string filename = TempFile::getName("string-array");
  this->zstd_compress_to(original, filename);

  StringArray<> copy;
  this->zstd_decompress_from(copy, filename);
  ASSERT_EQ(copy, original) << "Zstd compression changed the array with empty strings";

  TempFile::remove(filename);
}

TEST_F(StringArrayTest, DuplicateWithEmptyStrings)
{
  std::vector<std::string> source
  {
    "first",
    "second",
    "",
    "fourth",
    ""
  };
  StringArray<> original(source);
  StringArray<> truth = duplicate_array(source);

  std::string filename = TempFile::getName("string-array");
  sdsl::simple_sds::serialize_to(original, filename);

  StringArray<> copy;
  this->simple_sds_duplicate_from(copy, filename, reverse_string);
  ASSERT_EQ(copy, truth) << "Simple-SDS serialization changed the array with empty strings";

  TempFile::remove(filename);
}

TEST_F(StringArrayTest, ZstdDuplicateWithEmptyStrings)
{
  std::vector<std::string> source
  {
    "first",
    "second",
    "",
    "fourth",
    ""
  };
  StringArray<> original = duplicate_array(source);

  std::string filename = TempFile::getName("string-array");
  this->zstd_compress_even_to(original, filename);

  StringArray<> copy;
  this->zstd_decompress_duplicate_from(copy, filename, reverse_string);
  ASSERT_EQ(copy, original) << "Zstd compression changed the array with empty strings";

  TempFile::remove(filename);
}

TEST_F(StringArrayTest, ZstdDuplicateWithShortOrEmptyStrings)
{
  std::vector<std::string> source
  {
    "f",
    "s",
    "",
    "f",
    ""
  };
  StringArray<> original = duplicate_array(source);

  std::string filename = TempFile::getName("string-array");
  this->zstd_compress_even_to(original, filename);

  StringArray<> copy;
  this->zstd_decompress_duplicate_from(copy, filename, reverse_string);
  ASSERT_EQ(copy, original) << "Zstd compression changed the array with empty strings";

  TempFile::remove(filename);
}

//------------------------------------------------------------------------------

#if defined(GBWT_ENABLE_SHARED_MEMORY)

class StringArraySharedMemoryTest : public ::testing::Test
{
public:
  std::string segment_name;

  void SetUp() override
  {
    this->segment_name = "gbwt_test_string_array_" + std::to_string(sdsl::util::pid());
    bi::shared_memory_object::remove(this->segment_name.c_str());
  }

  void TearDown() override
  {
    bi::shared_memory_object::remove(this->segment_name.c_str());
  }
};

/*
  Which members an instantiation has is part of StringArray's interface, so it
  is checked here rather than left to whatever the rest of the tests happen to
  name. A segment is only meaningful for the allocator that puts characters in
  one, and building content from a source is only meaningful for an allocator
  that does not need to be told where to put it.
*/

static_assert(std::constructible_from<StringArray<SharedMemCharAllocatorType>, bi::managed_shared_memory*, std::string>,
              "The shared-memory allocator should accept a segment and a prefix");
static_assert(!std::constructible_from<StringArray<>, bi::managed_shared_memory*, std::string>,
              "The plain allocator should not accept a segment and a prefix");

static_assert(std::constructible_from<StringArray<>, std::vector<std::string>>,
              "The plain allocator should build an array from a vector of strings");
static_assert(!std::constructible_from<StringArray<SharedMemCharAllocatorType>, std::vector<std::string>>,
              "The shared-memory allocator should not build an array without being given a segment");

// A requires-expression only turns a failed requirement into `false` while
// substituting template arguments, so these are written against a parameter
// rather than naming the instantiations directly.
template<typename ArrayType>
concept HasAttach = requires(bi::managed_shared_memory* segment) { ArrayType::attach(segment, "arr"); };
template<typename ArrayType>
concept HasSimpleSdsDecompress = requires(ArrayType& array, std::istream& in) { array.simple_sds_decompress(in); };
template<typename ArrayType, typename OtherArrayType>
concept HasSwapWith = requires(ArrayType& array, OtherArrayType& another) { array.swap(another); };

static_assert(HasAttach<StringArray<SharedMemCharAllocatorType>>,
              "The shared-memory allocator should have attach()");
static_assert(!HasAttach<StringArray<>>,
              "The plain allocator should not have attach()");

static_assert(HasSimpleSdsDecompress<StringArray<>>,
              "The plain allocator should have simple_sds_decompress()");
static_assert(!HasSimpleSdsDecompress<StringArray<SharedMemCharAllocatorType>>,
              "The shared-memory allocator should not have simple_sds_decompress()");

// Exchanging representations only makes sense between arrays that free their
// characters the same way.
static_assert(HasSwapWith<StringArray<>, StringArray<>>,
              "swap() should work between arrays of the same type");
static_assert(!HasSwapWith<StringArray<>, StringArray<SharedMemCharAllocatorType>>,
              "swap() should not cross allocators");

// A second StringArray attaching to a name a prior one already published
// under, from the same process, should find and reuse that data instead of
// trying to construct it again under the same name.
TEST_F(StringArraySharedMemoryTest, AttachAfterPublish)
{
  std::vector<std::string> source { "first", "second", "third" };
  bi::managed_shared_memory segment(bi::create_only, this->segment_name.c_str(), 65536);

  StringArray<SharedMemCharAllocatorType> writer(source, &segment, "arr");
  ASSERT_EQ(writer.size(), source.size()) << "Writer has the wrong size";

  StringArray<SharedMemCharAllocatorType> reader = StringArray<SharedMemCharAllocatorType>::attach(&segment, "arr");
  ASSERT_EQ(reader.size(), source.size()) << "Reader did not attach to the published data";
  for(size_type i = 0; i < source.size(); i++)
  {
    EXPECT_EQ(reader.str(i), source[i]) << "Wrong string " << i << " after attaching";
  }
}

// attach() must fail instead of silently returning an empty array when
// nothing has been published under the given name yet.
TEST_F(StringArraySharedMemoryTest, AttachToMissingNameFails)
{
  bi::managed_shared_memory segment(bi::create_only, this->segment_name.c_str(), 65536);
  ASSERT_THROW(StringArray<SharedMemCharAllocatorType>::attach(&segment, "arr"), std::runtime_error)
    << "Attaching to a name nothing was published under should fail instead of silently succeeding";
}

// Unlike attach(), the lazy (non-attach) shared-memory constructor must
// never throw, regardless of `shared_memory` or `object_prefix_in_shared_memory`.
TEST_F(StringArraySharedMemoryTest, LazyConstructorNeverThrows)
{
  StringArray<SharedMemCharAllocatorType> no_segment;
  EXPECT_TRUE(no_segment.empty()) << "A StringArray with no segment should be empty";

  bi::managed_shared_memory segment(bi::create_only, this->segment_name.c_str(), 65536);
  StringArray<SharedMemCharAllocatorType> unpublished_name(&segment, "never_published");
  EXPECT_TRUE(unpublished_name.empty()) << "A StringArray naming an unpublished object should be empty, not throw";
}

// The lazy constructor never checks for existing data, even when a real
// segment already has something published under the given name: only
// attach() (or a load call) actually looks.
TEST_F(StringArraySharedMemoryTest, LazyConstructorIgnoresExistingPublication)
{
  std::vector<std::string> source { "first", "second" };
  bi::managed_shared_memory segment(bi::create_only, this->segment_name.c_str(), 65536);
  StringArray<SharedMemCharAllocatorType> writer(source, &segment, "arr");
  ASSERT_EQ(writer.size(), source.size()) << "Writer has the wrong size";

  StringArray<SharedMemCharAllocatorType> lazy(&segment, "arr");
  EXPECT_TRUE(lazy.empty()) << "The plain constructor should stay empty even though \"arr\" is already published";
}

// Two independent handles to the same segment stand in for two separate
// processes: the second must attach to the first's published data instead
// of decompressing its own copy under the same name, since republishing
// under a name that already exists would throw (see managed_shared_memory's
// construct<>() semantics).
TEST_F(StringArraySharedMemoryTest, SimpleSDSLoadDuplicateThenAttach)
{
  std::vector<std::string> source { "first", "second", "third", "fourth" };
  StringArray<> truth = duplicate_array(source);

  std::string filename = TempFile::getName("string-array");
  sdsl::simple_sds::serialize_to(StringArray<>(source), filename);

  bi::managed_shared_memory writer_segment(bi::create_only, this->segment_name.c_str(), 1024 * 1024);
  StringArray<SharedMemCharAllocatorType> writer(&writer_segment, "arr");
  {
    std::ifstream in(filename, std::ios_base::binary);
    writer.simple_sds_load_duplicate(in, reverse_string);
  }
  ASSERT_EQ(writer.size(), truth.size()) << "Writer has the wrong size after loading";
  for(size_type i = 0; i < truth.size(); i++)
  {
    EXPECT_EQ(writer.str(i), truth.str(i)) << "Writer has the wrong string " << i << " after loading";
  }

  bi::managed_shared_memory reader_segment(bi::open_only, this->segment_name.c_str());
  StringArray<SharedMemCharAllocatorType> reader(&reader_segment, "arr");
  {
    std::ifstream in(filename, std::ios_base::binary);
    reader.simple_sds_load_duplicate(in, reverse_string);
  }
  ASSERT_EQ(reader.size(), truth.size()) << "Reader did not attach to the published data";
  for(size_type i = 0; i < truth.size(); i++)
  {
    EXPECT_EQ(reader.str(i), truth.str(i)) << "Reader has the wrong string " << i << " after attaching";
  }

  TempFile::remove(filename);
}

// As SimpleSDSLoadDuplicateThenAttach, but for the zstd-compressed format
// and simple_sds_decompress_duplicate().
TEST_F(StringArraySharedMemoryTest, ZstdDecompressDuplicateThenAttach)
{
  std::vector<std::string> source { "first", "second", "third", "fourth" };
  StringArray<> truth = duplicate_array(source);

  std::string filename = TempFile::getName("string-array");
  {
    std::ofstream out(filename, std::ios_base::binary);
    truth.simple_sds_compress_even(out);
  }

  bi::managed_shared_memory writer_segment(bi::create_only, this->segment_name.c_str(), 1024 * 1024);
  StringArray<SharedMemCharAllocatorType> writer(&writer_segment, "arr");
  {
    std::ifstream in(filename, std::ios_base::binary);
    writer.simple_sds_decompress_duplicate(in, reverse_string);
  }
  ASSERT_EQ(writer.size(), truth.size()) << "Writer has the wrong size after decompressing";
  for(size_type i = 0; i < truth.size(); i++)
  {
    EXPECT_EQ(writer.str(i), truth.str(i)) << "Writer has the wrong string " << i << " after decompressing";
  }

  bi::managed_shared_memory reader_segment(bi::open_only, this->segment_name.c_str());
  StringArray<SharedMemCharAllocatorType> reader(&reader_segment, "arr");
  {
    std::ifstream in(filename, std::ios_base::binary);
    reader.simple_sds_decompress_duplicate(in, reverse_string);
  }
  ASSERT_EQ(reader.size(), truth.size()) << "Reader did not attach to the published data";
  for(size_type i = 0; i < truth.size(); i++)
  {
    EXPECT_EQ(reader.str(i), truth.str(i)) << "Reader has the wrong string " << i << " after attaching";
  }

  TempFile::remove(filename);
}

//------------------------------------------------------------------------------

// Copying between the two allocators is the only way to get characters from
// one kind of storage to the other, so the tests below cover each direction
// and each way it can fail. They are here rather than with the other
// StringArray tests because two allocators to copy between only exist when
// shared memory is compiled in.

// Reading a shared-memory array into an ordinary one has to produce an
// independent heap copy, not a second reference to the segment.
TEST_F(StringArraySharedMemoryTest, CopyFromSharedToPlain)
{
  std::vector<std::string> source { "first", "second", "third" };
  bi::managed_shared_memory segment(bi::create_only, this->segment_name.c_str(), 65536);
  StringArray<SharedMemCharAllocatorType> shared(source, &segment, "arr");

  StringArray<> constructed(shared);
  ASSERT_EQ(constructed.size(), source.size()) << "Copy-constructed array has the wrong size";
  for(size_type i = 0; i < source.size(); i++)
  {
    EXPECT_EQ(constructed.str(i), source[i]) << "Copy-constructed array has the wrong string " << i;
  }

  StringArray<> assigned(std::vector<std::string> { "replaced" });
  assigned = shared;
  ASSERT_EQ(assigned.size(), source.size()) << "Assigned array kept its old size";
  for(size_type i = 0; i < source.size(); i++)
  {
    EXPECT_EQ(assigned.str(i), source[i]) << "Assigned array has the wrong string " << i;
  }
}

// Assigning into an array that names a prefix nothing is published under yet
// publishes the characters there, so that another handle can attach to them.
TEST_F(StringArraySharedMemoryTest, AssignFromPlainToShared)
{
  std::vector<std::string> source { "first", "second", "third" };
  StringArray<> plain(source);
  bi::managed_shared_memory segment(bi::create_only, this->segment_name.c_str(), 65536);

  StringArray<SharedMemCharAllocatorType> shared(&segment, "arr");
  shared = plain;
  ASSERT_EQ(shared.size(), source.size()) << "Assigned shared-memory array has the wrong size";
  for(size_type i = 0; i < source.size(); i++)
  {
    EXPECT_EQ(shared.str(i), source[i]) << "Assigned shared-memory array has the wrong string " << i;
  }

  StringArray<SharedMemCharAllocatorType> reader = StringArray<SharedMemCharAllocatorType>::attach(&segment, "arr");
  EXPECT_EQ(reader, plain) << "Assigning into a shared-memory array did not publish the strings";
}

// Whatever is published under a prefix may be in use by another process, so
// assignment must refuse to replace it rather than doing so silently.
TEST_F(StringArraySharedMemoryTest, AssignOverPublishedSharedFails)
{
  std::vector<std::string> source { "first", "second" };
  bi::managed_shared_memory segment(bi::create_only, this->segment_name.c_str(), 65536);
  StringArray<SharedMemCharAllocatorType> writer(source, &segment, "arr");

  StringArray<SharedMemCharAllocatorType> second_handle(&segment, "arr");
  StringArray<> replacement(std::vector<std::string> { "different" });
  ASSERT_THROW(second_handle = replacement, std::runtime_error)
    << "Assigning over already published strings should fail instead of replacing them";

  StringArray<SharedMemCharAllocatorType> reader = StringArray<SharedMemCharAllocatorType>::attach(&segment, "arr");
  EXPECT_EQ(reader, writer) << "The published strings should be untouched after the failed assignment";
}

// An array with no segment has nowhere to put the characters.
TEST_F(StringArraySharedMemoryTest, AssignToSegmentlessSharedFails)
{
  StringArray<SharedMemCharAllocatorType> no_segment;
  StringArray<> plain(std::vector<std::string> { "first", "second" });
  ASSERT_THROW(no_segment = plain, std::runtime_error)
    << "Assigning into a shared-memory array with no segment should fail";
}

// Copying between two shared-memory arrays hands over the published objects
// instead of duplicating them in the segment.
TEST_F(StringArraySharedMemoryTest, CopyBetweenSharedArraysShares)
{
  std::vector<std::string> source { "first", "second" };
  bi::managed_shared_memory segment(bi::create_only, this->segment_name.c_str(), 65536);
  StringArray<SharedMemCharAllocatorType> writer(source, &segment, "arr");

  StringArray<SharedMemCharAllocatorType> copy(writer);
  EXPECT_EQ(copy.strings, writer.strings) << "Copying a shared-memory array should not copy the characters";
  EXPECT_EQ(copy, writer) << "The copy should hold the same strings";
}

// An rvalue source cannot hand its characters over to a different allocator,
// so it is copied and stays readable afterwards.
TEST_F(StringArraySharedMemoryTest, MoveAcrossAllocatorsCopies)
{
  std::vector<std::string> source { "first", "second" };
  bi::managed_shared_memory segment(bi::create_only, this->segment_name.c_str(), 65536);
  StringArray<SharedMemCharAllocatorType> shared(source, &segment, "arr");

  StringArray<> plain;
  plain = std::move(shared);
  ASSERT_EQ(plain.size(), source.size()) << "Moving from a shared-memory array gave the wrong size";
  EXPECT_EQ(plain, shared) << "The moved-from array should still hold the same strings";
}

// Comparison works whichever allocator either side keeps its characters in.
TEST_F(StringArraySharedMemoryTest, ComparisonAcrossAllocators)
{
  std::vector<std::string> source { "first", "second", "third" };
  bi::managed_shared_memory segment(bi::create_only, this->segment_name.c_str(), 65536);
  StringArray<SharedMemCharAllocatorType> shared(source, &segment, "arr");
  StringArray<> same(source);
  StringArray<> shorter(std::vector<std::string> { "first", "second" });
  StringArray<> different(std::vector<std::string> { "first", "second", "fourth" });

  EXPECT_TRUE(shared == same) << "Equal arrays with different allocators should compare equal";
  EXPECT_TRUE(same == shared) << "Comparison should not depend on which side is which";
  EXPECT_TRUE(shared != shorter) << "Arrays of different sizes should not compare equal";
  EXPECT_TRUE(shorter != shared) << "Arrays of different sizes should not compare equal either way";
  EXPECT_TRUE(shared != different) << "Arrays with different characters should not compare equal";
  EXPECT_TRUE(different != shared) << "Arrays with different characters should not compare equal either way";
}

// Two arrays with the same characters but different string boundaries hold
// different strings, so the index has to be part of the comparison.
TEST_F(StringArraySharedMemoryTest, ComparisonAcrossAllocatorsUsesIndex)
{
  bi::managed_shared_memory segment(bi::create_only, this->segment_name.c_str(), 65536);
  StringArray<SharedMemCharAllocatorType> shared(std::vector<std::string> { "ab", "c" }, &segment, "arr");
  StringArray<> regrouped(std::vector<std::string> { "a", "bc" });

  EXPECT_TRUE(shared != regrouped) << "Arrays splitting the same characters differently should not compare equal";
}

#endif

//------------------------------------------------------------------------------

TEST(DictionaryTest, Empty)
{
  Dictionary empty;

  EXPECT_EQ(empty.size(), static_cast<size_type>(0)) << "Empty dictionary contains " << empty.size() << " keys";
  EXPECT_TRUE(empty.empty()) << "Empty dictionary is not empty";

  size_t offset = empty.find("key");
  EXPECT_EQ(offset, empty.size()) << "Missing keys are not reported as missing";
}

TEST(DictionaryTest, Keys)
{
  std::vector<std::string> keys
  {
    "first", "second", "third", "fourth", "fifth"
  };

  Dictionary dict(keys);
  ASSERT_EQ(dict.size(), keys.size()) << "Expected " << keys.size() << " keys, got " << dict.size();
  EXPECT_FALSE(dict.empty()) << "The dictionary is empty";

  bool ok = true;
  for(size_type i = 0; i < keys.size(); i++)
  {
    ok &= (dict[i] == keys[i]);
    ok &= (dict.find(keys[i]) == i);
  }
  EXPECT_TRUE(ok) << "The dictionary does not contain the correct keys";

  size_t offset = dict.find("key");
  EXPECT_EQ(offset, dict.size()) << "Missing keys are not reported as missing";

  dict.remove(keys.size());
  ASSERT_EQ(dict.size(), keys.size()) << "Removing an invalid key changed Dictionary size";

  constexpr size_type REMOVED_KEY = 2;
  dict.remove(REMOVED_KEY);
  ASSERT_EQ(dict.size(), keys.size() - 1) << "Expected " << (keys.size() - 1) << " keys after removal, got " << dict.size();

  ok = true;
  for(size_type i = 0; i < keys.size(); i++)
  {
    if(i < REMOVED_KEY)
    {
      ok &= (dict[i] == keys[i]);
      ok &= (dict.find(keys[i]) == i);
    }
    else if(i == REMOVED_KEY)
    {
      ok &= (dict.find(keys[i]) == dict.size());
    }
    else
    {
      ok &= (dict[i - 1] == keys[i]);
      ok &= (dict.find(keys[i]) == i - 1);
    }
  }
  EXPECT_TRUE(ok) << "The dictionary does not contain the correct keys after removal";
}

TEST(DictionaryTest, Comparisons)
{
  std::vector<std::string> keys
  {
    "first", "second", "third", "fourth", "fifth"
  };
  std::vector<std::string> first_keys
  {
    "first", "second", "third"
  };
  std::vector<std::string> second_keys
  {
    "fourth", "fifth"
  };
  Dictionary empty, all(keys), first(first_keys), second(second_keys);

  EXPECT_NE(empty, all) << "Empty dictionary is equal to the full dictionary";
  EXPECT_NE(empty, first) << "Empty dictionary is equal to the first dictionary";
  EXPECT_NE(empty, second) << "Empty dictionary is equal to the second dictionary";
  EXPECT_NE(all, first) << "Full dictionary is equal to the first dictionary";
  EXPECT_NE(all, second) << "Full dictionary is equal to the second dictionary";
  EXPECT_NE(first, second) << "The first and second dictionaries are equal";

  empty.append(first);
  EXPECT_EQ(empty, first) << "Appending to an empty dictionary does not work";

  first.append(second);
  EXPECT_EQ(first, all) << "Appending to a non-empty dictionary does not work";
}

TEST(DictionaryTest, Merging)
{
  std::vector<std::string> keys
  {
    "first", "second", "third", "fourth", "fifth"
  };
  std::vector<std::string> first_keys
  {
    "first", "second", "third"
  };
  std::vector<std::string> second_keys
  {
    "fifth", "first", "fourth"
  };

  Dictionary first(first_keys), second(second_keys);
  Dictionary merged(first, second);

  EXPECT_EQ(merged.size(), keys.size()) << "Expected " << keys.size() << " keys, got " << merged.size();
  for(const std::string& key : keys)
  {
    EXPECT_LT(merged.find(key), merged.size()) << "The dictionary does not contain " << key;
  }
}

TEST(DictionaryTest, Serialization)
{
  std::vector<std::string> keys
  {
    "first", "second", "third", "fourth", "fifth"
  };
  Dictionary original(keys);
  size_t expected_size = original.simple_sds_size() * sizeof(sdsl::simple_sds::element_type);

  std::string sdsl_filename = TempFile::getName("Dictionary");
  sdsl::store_to_file(original, sdsl_filename);
  Dictionary sdsl_copy; sdsl::load_from_file(sdsl_copy, sdsl_filename);
  TempFile::remove(sdsl_filename);
  EXPECT_EQ(original, sdsl_copy) << "SDSL serialization failed";

  std::string simple_sds_filename = TempFile::getName("Dictionary");
  sdsl::simple_sds::serialize_to(original, simple_sds_filename);
  std::ifstream in(simple_sds_filename, std::ios_base::binary);
  size_t bytes = fileSize(in);
  ASSERT_EQ(bytes, expected_size) << "Invalid Simple-SDS file size";
  Dictionary simple_sds_copy; simple_sds_copy.simple_sds_load(in);
  in.close();
  TempFile::remove(simple_sds_filename);
  EXPECT_EQ(original, simple_sds_copy) << "Simple-SDS serialization failed";
}

//------------------------------------------------------------------------------

class TagsTest : public ::testing::Test
{
public:
  void check_tags(const Tags& tags, const std::map<std::string, std::string>& truth) const
  {
    ASSERT_EQ(tags.size(), truth.size()) << "Incorrect tags size";
    ASSERT_EQ(tags.empty(), truth.empty()) << "Incorrect emptiness";

    for(auto iter = truth.begin(); iter != truth.end(); ++iter)
    {
      EXPECT_TRUE(tags.contains(iter->first)) << "Key " << iter->first << " is missing";
      EXPECT_EQ(tags.get(iter->first), iter->second) << "Invalid value for key " << iter->first;
    }
  }

  void check_file_size(const Tags& original, std::ifstream& in) const
  {
    size_type expected_size = original.simple_sds_size() * sizeof(sdsl::simple_sds::element_type);
    size_type file_size = fileSize(in);
    ASSERT_EQ(expected_size, file_size) << "Incorrect file size";
  }
};

TEST_F(TagsTest, Empty)
{
  std::map<std::string, std::string> truth;
  Tags tags;
  this->check_tags(tags, truth);
}

TEST_F(TagsTest, NonEmpty)
{
  std::map<std::string, std::string> truth
  {
    { "first-key", "first-value" },
    { "second-key", "second-value" },
    { "third-key", "third-value" },
  };
  Tags tags;
  for(auto iter = truth.begin(); iter != truth.end(); ++iter) { tags.set(iter->first, iter->second); }
  this->check_tags(tags, truth);
}

TEST_F(TagsTest, MissingKeys)
{
  std::map<std::string, std::string> truth
  {
    { "first-key", "first-value" },
    { "second-key", "second-value" },
    { "third-key", "third-value" },
  };
  Tags tags;
  for(auto iter = truth.begin(); iter != truth.end(); ++iter) { tags.set(iter->first, iter->second); }

  ASSERT_FALSE(tags.contains("key")) << "Tags contains an invalid key";
  ASSERT_TRUE(tags.get("key").empty()) << "Non-empty value for an invalid key";
}

TEST_F(TagsTest, RemoveKey)
{
  std::map<std::string, std::string> truth
  {
    { "first-key", "first-value" },
    { "second-key", "second-value" },
    { "third-key", "third-value" },
  };
  Tags tags;
  for(auto iter = truth.begin(); iter != truth.end(); ++iter) { tags.set(iter->first, iter->second); }

  tags.unset("second-key");
  truth.erase("second-key");
  this->check_tags(tags, truth);

  ASSERT_FALSE(tags.contains("second-key")) << "Tags still contains the removed key";
  ASSERT_TRUE(tags.get("second-key").empty()) << "Non-empty value for the removed key";
}

TEST_F(TagsTest, NormalizedKeys)
{
  std::map<std::string, std::string> source
  {
    { "First-Key", "first-value" },
    { "Second-Key", "second-value" },
    { "Third-Key", "third-value" },
  };
  Tags tags;
  for(auto iter = source.begin(); iter != source.end(); ++iter) { tags.set(iter->first, iter->second); }
  this->check_tags(tags, source); // Check with the original keys.

  std::map<std::string, std::string> truth
  {
    { "first-key", "first-value" },
    { "second-key", "second-value" },
    { "third-key", "third-value" },
  };
  this->check_tags(tags, truth); // Check with normalized keys.
}

TEST_F(TagsTest, SerializeEmpty)
{
  Tags original;

  std::string filename = TempFile::getName("tags");
  sdsl::store_to_file(original, filename);

  Tags copy; sdsl::load_from_file(copy, filename);
  ASSERT_EQ(copy, original) << "Serialization changed empty tags";

  TempFile::remove(filename);
}

TEST_F(TagsTest, CompressEmpty)
{
  Tags original;

  std::string filename = TempFile::getName("tags");
  sdsl::simple_sds::serialize_to(original, filename);

  Tags copy;
  std::ifstream in(filename, std::ios_base::binary);
  this->check_file_size(original, in);
  copy.simple_sds_load(in);
  in.close();
  ASSERT_EQ(copy, original) << "Compression changed empty tags";

  TempFile::remove(filename);
}

TEST_F(TagsTest, SerializeNonEmpty)
{
  std::map<std::string, std::string> truth
  {
    { "first-key", "first-value" },
    { "second-key", "second-value" },
    { "third-key", "third-value" },
  };
  Tags original;
  for(auto iter = truth.begin(); iter != truth.end(); ++iter) { original.set(iter->first, iter->second); }

  std::string filename = TempFile::getName("tags");
  sdsl::store_to_file(original, filename);

  Tags copy; sdsl::load_from_file(copy, filename);
  ASSERT_EQ(copy, original) << "Serialization changed non-empty tags";

  TempFile::remove(filename);
}

TEST_F(TagsTest, CompressNonEmpty)
{
  std::map<std::string, std::string> truth
  {
    { "first-key", "first-value" },
    { "second-key", "second-value" },
    { "third-key", "third-value" },
  };
  Tags original;
  for(auto iter = truth.begin(); iter != truth.end(); ++iter) { original.set(iter->first, iter->second); }

  std::string filename = TempFile::getName("tags");
  sdsl::simple_sds::serialize_to(original, filename);

  Tags copy;
  std::ifstream in(filename, std::ios_base::binary);
  this->check_file_size(original, in);
  copy.simple_sds_load(in);
  in.close();
  ASSERT_EQ(copy, original) << "Compression changed non-empty tags";

  TempFile::remove(filename);
}

//------------------------------------------------------------------------------

} // namespace
