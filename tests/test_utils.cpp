#include <gtest/gtest.h>

#include <cstring>
#include <random>

#include <gbwt/utils.h>

using namespace gbwt;

namespace
{

//------------------------------------------------------------------------------

class ZstdTest : public ::testing::Test
{
public:
  std::vector<char> random_bytes(size_t size, unsigned seed) const
  {
    std::mt19937_64 rng(seed);
    std::vector<char> data(size);
    for(size_t i = 0; i < size; i += sizeof(std::uint64_t))
    {
      std::uint64_t rand_val = rng();
      size_t chunk_size = std::min(sizeof(std::uint64_t), size - i);
      std::memcpy(data.data() + i, &rand_val, chunk_size);
    }
    return data;
  }

  void compress_direct(const std::vector<char>& input, std::vector<char>& output, const std::string& test_case) const
  {
    ZstdCompressor compressor;
    ASSERT_NO_THROW(compressor.compressDirect(std::string_view(input.data(), input.size())))
      << "compressDirect() failed for " << test_case;
    ASSERT_NO_THROW(compressor.finish()) << "finish() failed for " << test_case;
    output = compressor.outputData();
  }

  void compress_incremental(const std::vector<char>& input, std::vector<char>& output, size_t chunk_size, const std::string& test_case) const
  {
    ZstdCompressor compressor;
    size_t offset = 0;
    while(offset < input.size())
    {
      size_t to_compress = std::min(chunk_size, input.size() - offset);
      ASSERT_NO_THROW(compressor.compress(std::string_view(input.data() + offset, to_compress)))
        << "compress() failed for " << test_case << " at offset " << offset;
      offset += to_compress;
    }
    ASSERT_NO_THROW(compressor.finish()) << "finish() failed for " << test_case;
    output = compressor.outputData();
  }

  void try_compress(const std::vector<char>& input, std::vector<char>& output, size_t chunk_size, const std::string& test_case) const
  {
    std::vector<char> direct, incremental;
    this->compress_direct(input, direct, test_case);
    this->compress_incremental(input, incremental, chunk_size, test_case);
    ASSERT_EQ(direct, incremental) << "Mismatch between direct and incremental compression for " << test_case;
    output = direct;
    ASSERT_TRUE(output.size() > 0) << "No data produced for " << test_case;
  }

  void decompress_direct(const std::vector<char>& input, size_t original_size, std::vector<char>& output, const std::string& test_case) const
  {
    std::vector<char> copy = input;
    ZstdDecompressor decompressor(std::move(copy));
    ASSERT_NO_THROW(decompressor.decompress(original_size, output))
      << "decompress() failed for " << test_case;
    ASSERT_TRUE(decompressor.finished()) << "Decompressor not finished after direct decompression for " << test_case;
    ASSERT_THROW(decompressor.decompress(1, output), sdsl::simple_sds::InvalidData)
      << "Decompressor did not throw on decompressing more data than available for " << test_case;
  }

  void decompress_incremental(const std::vector<char>& input, size_t original_size, std::vector<char>& output, size_t chunk_size, const std::string& test_case) const
  {
    std::vector<char> copy = input;
    ZstdDecompressor decompressor(std::move(copy));
    size_t decompressed = 0;
    do
    {
      // Even with an empty input, there are still some bytes to decompress.
      size_t to_decompress = std::min(chunk_size, original_size - decompressed);
      ASSERT_NO_THROW(decompressor.decompress(to_decompress, output))
        << "decompress() failed for " << test_case << " at offset " << decompressed;
      decompressed += to_decompress;
    }
    while(decompressed < original_size);
    ASSERT_TRUE(decompressor.finished()) << "Decompressor not finished after incremental decompression for " << test_case;
    ASSERT_THROW(decompressor.decompress(1, output), sdsl::simple_sds::InvalidData)
      << "Decompressor did not throw on decompressing more data than available for " << test_case;
  }

  void try_decompress(const std::vector<char>& input, const std::vector<char>& original, size_t chunk_size, const std::string& test_case) const
  {
    std::vector<char> direct, incremental;
    this->decompress_direct(input, original.size(), direct, test_case);
    this->decompress_incremental(input, original.size(), incremental, chunk_size, test_case);
    ASSERT_EQ(direct, incremental) << "Mismatch between direct and incremental decompression for " << test_case;
    ASSERT_EQ(direct, original) << "Decompressed data does not match original for " << test_case;
  }
};

TEST_F(ZstdTest, NoData)
{
  std::vector<char> empty, compressed;
  this->try_compress(empty, compressed, 128, "no data");
  this->try_decompress(compressed, empty, 128, "no data");
}

TEST_F(ZstdTest, SmallRandomData)
{
  std::vector<char> original = this->random_bytes(1024, 42);
  std::vector<char> compressed;
  this->try_compress(original, compressed, 128, "small random data");
  this->try_decompress(compressed, original, 128, "small random data");
}

TEST_F(ZstdTest, SmallTexts)
{
  std::vector<std::string> files
  {
    "Makefile",
    "test_bwtmerge.cpp",
    "test_metadata.cpp",
    "test_queries.cpp",
    "test_support.cpp",
    "test_utils.cpp",
    "test_variants.cpp"
  };

  for(const std::string& filename : files)
  {
    std::ifstream in(filename, std::ios_base::binary);
    ASSERT_TRUE(in) << "Cannot open " << filename;
    std::vector<char> original((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    in.close();

    std::vector<char> compressed;
    this->try_compress(original, compressed, 128, filename);
    this->try_decompress(compressed, original, 128, filename);
  }
}

TEST_F(ZstdTest, LargeRandomData)
{
  std::vector<char> original = this->random_bytes(10 * 1024 * 1024, 4242); // 10 MB
  std::vector<char> compressed;
  this->try_compress(original, compressed, 8192, "large random data");
  this->try_decompress(compressed, original, 8192, "large random data");
}

//------------------------------------------------------------------------------

} // namespace
