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

#include <gtest/gtest.h>

#include <random>

#include <gbwt/bwtmerge.h>

using namespace gbwt;

namespace
{

//------------------------------------------------------------------------------

void
initLargeArray(std::vector<edge_type>& large_array, size_type seed = 0xDEADBEEF)
{
  constexpr size_type TOTAL_VALUES = 2 * BlockArray::BLOCK_SIZE;
  constexpr size_type NODES = MILLION;
  constexpr size_type OFFSETS = 1000;

  large_array.clear();
  large_array.reserve(TOTAL_VALUES);
  std::mt19937_64 rng(seed);
  for(size_type i = 0; i < TOTAL_VALUES; i++)
  {
    large_array.emplace_back(rng() % NODES, rng() % OFFSETS);
  }
}

//------------------------------------------------------------------------------

class BlockArrayTest : public ::testing::Test
{
public:
  std::vector<size_type> block_offsets;
  size_type last_block_size, last_block_offset;

  BlockArrayTest()
  {
  }

  void SetUp() override
  {
    this->block_offsets = { 0, 42, 51 };
    this->last_block_size = 37;
  }

  size_type expectedSize() { return (this->block_offsets.size() - 1) * BlockArray::BLOCK_SIZE + this->last_block_size; }
  size_type expectedBlocks() { return this->block_offsets.size(); }

  // Fill the initial array.
  void fillArray(BlockArray& array)
  {
    // Fill blocks with repeats of 0..255 to ensure the correct byte order.
    // Use various offsets to distinguish the blocks.
    for(size_type block = 0; block < this->expectedBlocks(); block++)
    {
      size_type limit = (block + 1 == this->expectedBlocks() ? this->last_block_size : BlockArray::BLOCK_SIZE);
      for(size_type i = 0; i < limit; i++) { array.push_back((i + this->block_offsets[block]) & 0xFF); }
    }
  }

  // Fill the last block.
  void fillLastBlock(BlockArray& array)
  {
    size_type limit = array.blocks() * BlockArray::BLOCK_SIZE;
    for(size_type i = array.size(); i < limit; i++) { array.push_back((i + this->block_offsets.back()) & 0xFF); }
  }

  size_type wrongValues(const BlockArray& array)
  {
    size_type wrong_values = 0;
    for(size_type i = 0; i < array.size(); i++)
    {
      if(array[i] != ((array.offset(i) + this->block_offsets[array.block(i)]) & 0xFF)) { wrong_values++; }
    }
    return wrong_values;
  }

  size_type wrongValues(const BlockArray& array, sdsl::int_vector_buffer<8>& buffer)
  {
    size_type wrong_values = 0, limit = std::min(array.size(), buffer.size());
    for(size_type i = 0; i < limit; i++)
    {
      if(array[i] != buffer[i]) { wrong_values++; }
    }
    return wrong_values;
  }

  void testSerialization(const BlockArray& array, const std::string& filename, const std::string& test_name)
  {
    sdsl::store_to_file(array, filename);

    BlockArray new_array;
    sdsl::load_from_file(new_array, filename);
    EXPECT_EQ(new_array.size(), array.size()) << test_name << " array size changed after loading";
    size_type wrong_values = wrongValues(new_array);
    EXPECT_EQ(wrong_values, 0u) << test_name << " array contains " << wrong_values << " wrong values";

    sdsl::int_vector_buffer<8> buffer(filename);
    EXPECT_EQ(new_array.size(), buffer.size()) << test_name << " array size is different as sdsl::int_vector_buffer<8>";
    wrong_values = wrongValues(new_array, buffer);
    EXPECT_EQ(wrong_values, 0u) << test_name << " _array differs at " << wrong_values << " positions from sdsl::int_vector_buffer<8>";
  }
};

TEST_F(BlockArrayTest, BasicTests)
{
  // Empty array.
  BlockArray array;
  EXPECT_EQ(array.size(), 0u) << "Initial array size was non-zero";
  EXPECT_EQ(array.blocks(), 0u) << "Initial array block count was non-zero";
  EXPECT_TRUE(array.empty()) << "The initial array was not empty";

  // Array with values.
  std::vector<byte_type> values { 42, 51, 33 };
  for(auto value : values) { array.push_back(value); }
  ASSERT_EQ(array.size(), 3u) << "Array size was wrong after insertions";
  EXPECT_EQ(array.blocks(), 1u) << "Block count was wrong after insertions";
  EXPECT_FALSE(array.empty()) << "The array was empty after insertions";
  for(size_type i = 0; i < array.size(); i++)
  {
    EXPECT_EQ(array[i], values[i]) << "Wrong value at offset " << i;
  }

  // Array with modified values.
  std::vector<byte_type> new_values { 123, 22, 222 };
  for(size_type i = 0; i < array.size(); i++) { array[i] = new_values[i]; }
  ASSERT_EQ(array.size(), 3u) << "Array size was wrong after updates";
  EXPECT_EQ(array.blocks(), 1u) << "Block count was wrong after updates";
  EXPECT_FALSE(array.empty()) << "The array was empty after updates";
  for(size_type i = 0; i < array.size(); i++)
  {
    EXPECT_EQ(array[i], new_values[i]) << "Wrong value at offset " << i;
  }

  // Cleared array.
  array.clear();
  EXPECT_EQ(array.size(), 0u) << "Array size was non-zero after clear()";
  EXPECT_EQ(array.blocks(), 0u) << "Block count was non-zero after clear()";
  EXPECT_TRUE(array.empty()) << "The array was non-empty after clear()";
}

TEST_F(BlockArrayTest, LargeArray)
{
  // Large array.
  BlockArray array;
  fillArray(array);
  ASSERT_EQ(array.size(), expectedSize()) << "Wrong size for the large array";
  EXPECT_EQ(array.blocks(), expectedBlocks()) << "Wrong block count for the large array";
  size_type wrong_values = wrongValues(array);
  EXPECT_EQ(wrong_values, 0u) << "The array contains " << wrong_values << " wrong values";

  // Fill the last block with values.
  fillLastBlock(array);
  EXPECT_EQ(array.size(), expectedBlocks() * BlockArray::BLOCK_SIZE) << "Wrong size for the filled array";
  EXPECT_EQ(array.blocks(), expectedBlocks()) << "Wrong block count for the filled array";
  array.push_back(13);
  EXPECT_EQ(array.size(), expectedBlocks() * BlockArray::BLOCK_SIZE + 1) << "Wrong size after the final insertion";
  EXPECT_EQ(array.blocks(), expectedBlocks() + 1) << "Wrong block count after the final insertion";
}

TEST_F(BlockArrayTest, ClearArrayStart)
{
  BlockArray array;
  fillArray(array);

  // Clear individual blocks at block start.
  for(size_type block = 0; block < array.blocks(); block++)
  {
    EXPECT_NE(array.data[block], nullptr) << "Block " << block << " has already been cleared" << std::endl;
    array.clearUntil(block * BlockArray::BLOCK_SIZE);
    EXPECT_NE(array.data[block], nullptr) << "Block " << block << " was deleted" << std::endl;
    if(block > 0)
    {
      EXPECT_EQ(array.data[block - 1], nullptr) << "Block " << block << " was not deleted" << std::endl;
    }
  }
}

TEST_F(BlockArrayTest, ClearArrayEnd)
{
  BlockArray array;
  fillArray(array);

  // Clear individual blocks at block end.
  for(size_type block = 0; block < array.blocks(); block++)
  {
    EXPECT_NE(array.data[block], nullptr) << "Block " << block << " has already been cleared" << std::endl;
    array.clearUntil((block + 1) * BlockArray::BLOCK_SIZE - 1);
    EXPECT_NE(array.data[block], nullptr) << "Block " << block << " was deleted" << std::endl;
    if(block > 0)
    {
      EXPECT_EQ(array.data[block - 1], nullptr) << "Block " << block << " was not deleted" << std::endl;
    }
  }
}

TEST_F(BlockArrayTest, Serialization)
{
  BlockArray array;
  std::string filename = TempFile::getName("BlockArray");

  // Empty array.
  testSerialization(array, filename, "Empty");

  // Filled array.
  fillArray(array);
  testSerialization(array, filename, "Filled");

  TempFile::remove(filename);
}

//------------------------------------------------------------------------------

class GapArrayTest : public ::testing::Test
{
public:
  std::vector<edge_type> small_array;
  std::vector<edge_type> large_array;

  GapArrayTest()
  {
  }

  void SetUp() override
  {
    this->small_array = 
    {
      { 1, 2 }, { 1, 4 }, { 4, 2 }, { 4, 1 }, { 0, 0 }, { 3, 5 }, { 4, 2 }, { 6, 1 }, { 1, 2 }
    };
  }

  static void buildAndCheckArray(const std::vector<edge_type>& values, const std::string& test_name)
  {
    std::vector<edge_type> correct_values = values;
    GapArray<BlockArray> array(correct_values);
    checkArray(array, correct_values, test_name);
  }

  template<class ByteArray>
  static void checkArray(GapArray<ByteArray>& array, const std::vector<edge_type>& values, const std::string& test_name)
  {
    ASSERT_EQ(array.size(), values.size()) << test_name << ": wrong size";
    EXPECT_EQ(array.empty(), (array.size() == 0)) << test_name << ": wrong empty() result";

    typename GapArray<ByteArray>::iterator iter(array);
    checkIterator(iter, values, test_name);
  }

  template<class ByteArray>
  static void checkBuffered(GapArray<ByteArray>& array, const std::vector<edge_type>& values, const std::string& test_name)
  {
    ASSERT_EQ(array.size(), values.size()) << test_name << ": wrong size";
    EXPECT_EQ(array.empty(), (array.size() == 0)) << test_name << ": wrong empty() result";

    typename GapArray<ByteArray>::iterator iter(array);
    ProducerBuffer<typename GapArray<ByteArray>::iterator> buffer(iter, MEGABYTE);
    checkIterator(buffer, values, test_name);
  }

  template<class Iterator>
  static void checkIterator(Iterator& iter, const std::vector<edge_type>& values, const std::string& test_name)
  {
    std::vector<edge_type> correct_values = values;
    sequentialSort(correct_values.begin(), correct_values.end());

    bool early_end = false;
    size_type wrong_values = 0, end_at = 0;
    for(size_type i = 0; i < correct_values.size(); i++)
    {
      if(iter.end()) { early_end = true; end_at = i; break; }
      if(*iter != correct_values[i]) { wrong_values++; }
      ++iter;
    }
    EXPECT_FALSE(early_end) << test_name << ": Array ended early (" << end_at << " / " << correct_values.size() << ")";
    EXPECT_TRUE(iter.end()) << test_name << ": Array is too large";
    EXPECT_EQ(wrong_values, 0u) << test_name << ": " << wrong_values << " wrong values";
  }

  static void testSerialization(const std::vector<edge_type>& values, const std::string& filename, const std::string& test_name)
  {
    std::vector<edge_type> buffer = values;
    GapArray<BlockArray> output(buffer);
    output.write(filename);

    GapArray<sdsl::int_vector_buffer<8>> input;
    gbwt::open(input, filename, values.size());
    checkArray(input, values, test_name);
  }
};

TEST_F(GapArrayTest, BasicTests)
{
  // Empty array.
  std::vector<edge_type> empty_values;
  GapArray<BlockArray> empty_array;
  checkArray(empty_array, empty_values, "Empty");

  // Sorted array.
  std::vector<edge_type> sorted_values = small_array;
  sequentialSort(sorted_values.begin(), sorted_values.end());
  buildAndCheckArray(sorted_values, "Sorted");

  // Unsorted array.
  buildAndCheckArray(small_array, "Unsorted");
}

TEST_F(GapArrayTest, BufferedReading)
{
  std::vector<edge_type> buffer = small_array;
  GapArray<BlockArray> array(buffer);
  checkBuffered(array, small_array, "Buffered");
}

TEST_F(GapArrayTest, Serialization)
{
  std::string filename = TempFile::getName("GapArray");

  // Empty array.
  std::vector<edge_type> empty_values;
  testSerialization(empty_values, filename, "Empty");

  // Non-empty array.
  testSerialization(small_array, filename, "Non-empty");

  TempFile::remove(filename);
}

TEST_F(GapArrayTest, Merge)
{
  // Initialize vectors.
  std::vector<edge_type> first_array = small_array;
  std::vector<edge_type> second_array =
  {
    { 1, 3 }, { 1, 4 }, { 0, 1 }, { 3, 5 }, { 3, 0 }, { 3, 7 }, { 6, 0 }, { 5, 1 }, { 5, 3 }
  };
  std::vector<edge_type> correct_values;
  correct_values.insert(correct_values.end(), first_array.begin(), first_array.end());
  correct_values.insert(correct_values.end(), second_array.begin(), second_array.end());

  // Create and merge GapArrays.
  GapArray<BlockArray> first(first_array);
  GapArray<BlockArray> second(second_array);
  GapArray<BlockArray> merged(first, second);
  checkArray(merged, correct_values, "Merged");
}

TEST_F(GapArrayTest, MultipleFiles)
{
  // Two files.
  std::vector<std::string> filenames;
  filenames.push_back(TempFile::getName("GapArray"));
  filenames.push_back(TempFile::getName("GapArray"));
  std::vector<range_type> node_ranges =  { { 0, 3 }, { 4, 6 } };

  // Build and split the array.
  std::vector<edge_type> correct_values = small_array;
  GapArray<BlockArray> array(correct_values);
  std::vector<size_type> value_counts;
  array.write(filenames, node_ranges, value_counts);

  // Basic tests.
  ASSERT_EQ(value_counts.size(), filenames.size()) << "Wrong number of value counts";
  size_type total_values = 0;
  for(size_type count : value_counts) { total_values += count; }
  ASSERT_EQ(total_values, correct_values.size()) << "Total value count is wrong";

  // Test each file separately.
  std::vector<edge_type>::iterator array_iter = correct_values.begin();
  for(size_type file = 0; file < filenames.size(); file++)
  {
    GapArray<sdsl::int_vector_buffer<8>> part;
    open(part, filenames[file], value_counts[file]);
    GapArray<sdsl::int_vector_buffer<8>>::iterator part_iter(part);
    size_type values = 0, wrong_values = 0;
    bool early_end = false;
    while(array_iter != correct_values.end() && array_iter->first <= node_ranges[file].second)
    {
      if(part_iter.end()) { early_end = true; break; }
      if(*part_iter != *array_iter) { wrong_values++; }
      values++;
      ++array_iter; ++part_iter;
    }
    EXPECT_FALSE(early_end) << "File " << file << " ended early";
    EXPECT_TRUE(part_iter.end()) << "File " << file << " is too large";
    EXPECT_EQ(values, value_counts[file]) << "File " << file << ": wrong number of values";
    EXPECT_EQ(wrong_values, 0u) << "File " << file << ": " << wrong_values << " wrong values";

    // Skip the values that were missing from the file.
    while(array_iter != correct_values.end() && array_iter->first <= node_ranges[file].second) { ++array_iter; }
  }

  // Cleanup.
  for(std::string& filename : filenames) { TempFile::remove(filename); }
}

TEST_F(GapArrayTest, Large)
{
  // Create and check GapArray.
  initLargeArray(large_array);
  GapArray<BlockArray> array(large_array);
  checkArray(array, large_array, "Large");

  // Ensure that the array was cleared, except for the last block
  for(size_type block = 0; block + 1< array.data.blocks(); block++)
  {
    EXPECT_EQ(array.data.data[block], nullptr) << "Block " << block << " was not deleted";
  }
}

//------------------------------------------------------------------------------

class RankArrayTest : public ::testing::Test
{
public:
  std::vector<edge_type> first, second;

  RankArrayTest()
  {
    Verbosity::set(Verbosity::SILENT);
  }

  void SetUp() override
  {
    this->first = 
    {
      { 1, 2 }, { 1, 4 }, { 4, 2 }, { 4, 1 }, { 0, 0 }, { 3, 5 }, { 4, 2 }, { 6, 1 }, { 1, 2 }
    };
    this->second = 
    {
      { 1, 3 }, { 1, 4 }, { 0, 1 }, { 3, 5 }, { 3, 0 }, { 3, 7 }, { 6, 0 }, { 5, 1 }, { 5, 3 }
    };
  }

  static void addFile(RankArray& array, const std::vector<edge_type>& values)
  {
    std::vector<edge_type> buffer = values;
    GapArray<BlockArray> data(buffer);
    size_type file_num = 0;
    std::string filename;
    std::tie(filename, file_num) = array.addFile();
    array.addValueCount(file_num, data.size());
    data.write(filename); data.clear();
  }

  static void checkArray(RankArray& array, const std::vector<edge_type>& values, const std::string& test_name)
  {
    std::vector<edge_type> correct_values = values;
    sequentialSort(correct_values.begin(), correct_values.end());

    bool early_end = false;
    size_type wrong_values = 0, end_at = 0;
    array.open();
    for(size_type i = 0; i < correct_values.size(); i++)
    {
      if(array.end()) { early_end = true; end_at = i; break; }
      if(*array != correct_values[i]) { wrong_values++; }
      ++array;
    }
    EXPECT_FALSE(early_end) << test_name << ": Array ended early (" << end_at << " / " << correct_values.size() << ")";
    EXPECT_TRUE(array.end()) << test_name << ": Array is too large";
    EXPECT_EQ(wrong_values, 0u) << test_name << ": " << wrong_values << " wrong values";
    array.close();
  }

  // Note that this sorts correct_values.
  static void checkBuffer(RankArray& array, std::vector<edge_type>& correct_values, const std::string& test_name)
  {
    parallelQuickSort(correct_values.begin(), correct_values.end());

    bool early_end = false;
    size_type wrong_values = 0, end_at = 0;
    ProducerBuffer<RankArray> buffer(array);
    for(size_type i = 0; i < correct_values.size(); i++)
    {
      if(buffer.end()) { early_end = true; end_at = i; break; }
      if(*buffer != correct_values[i]) { wrong_values++; }
      ++buffer;
    }
    EXPECT_FALSE(early_end) << test_name << ": Array ended early (" << end_at << " / " << correct_values.size() << ")";
    EXPECT_TRUE(buffer.end()) << test_name << ": Array is too large";
    EXPECT_EQ(wrong_values, 0u) << test_name << ": " << wrong_values << " wrong values";
    array.close();
  }
};

TEST_F(RankArrayTest, BasicTests)
{
  // Empty array.
  std::vector<edge_type> empty_values;
  RankArray array;
  checkArray(array, empty_values, "Empty");

  // Single file.
  addFile(array, first);
  checkArray(array, first, "Single");

  // Multiple files.
  std::vector<edge_type> correct_values;
  correct_values.insert(correct_values.end(), first.begin(), first.end());
  correct_values.insert(correct_values.end(), second.begin(), second.end());
  addFile(array, second);
  checkArray(array, correct_values, "Multiple");
}

TEST_F(RankArrayTest, BufferedReading)
{
  // Empty array.
  std::vector<edge_type> correct_values;
  RankArray array;
  checkBuffer(array, correct_values, "Empty");

  // Single file.
  std::vector<edge_type> data;
  initLargeArray(data);
  addFile(array, data);
  correct_values.insert(correct_values.end(), data.begin(), data.end());
  data.clear();
  checkBuffer(array, correct_values, "Single");

  // Multiple files.
  initLargeArray(data, 0x42424242);
  addFile(array, data);
  correct_values.insert(correct_values.end(), data.begin(), data.end());
  data.clear();
  checkBuffer(array, correct_values, "Multiple");
}

TEST_F(RankArrayTest, MergeBuffers)
{
  // Merge parameters.
  MergeParameters parameters;
  parameters.setPosBufferSize(1);
  parameters.setThreadBufferSize(4);
  parameters.setMergeBuffers(2);
  parameters.setMergeJobs(2);

  // Input data and buffers.
  std::vector<std::vector<edge_type>> data(2);
  initLargeArray(data[0]);
  initLargeArray(data[1], 0x42424242);

  // Determine the node ranges.
  node_type max_node = 0;
  for(size_type i = 0; i < data.size(); i++)
  {
    for(size_type j = 0; j < data[i].size(); j++)
    {
      if(data[i][j].first > max_node) { max_node = data[i][j].first; }
    }
  }
  std::vector<range_type> node_ranges = Range::partition(range_type(0, max_node), parameters.merge_jobs);

  // Create the merge buffers.
  MergeBuffers buffers(data[0].size() + data[1].size(), data.size(), parameters, node_ranges);

  // Insert the data.
  //#pragma omp parallel for schedule(static)
  for(size_type thread = 0; thread < data.size(); thread++)
  {
    for(edge_type position : data[thread])
    {
      buffers.insert(position, thread);
    }
  }
  buffers.flush();

  // Generate the correct values.
  std::vector<edge_type> correct_values;
  for(size_type i = 0; i < data.size(); i++)
  {
    correct_values.insert(correct_values.end(), data[i].begin(), data[i].end());
    data[i].clear();
  }
  parallelQuickSort(correct_values.begin(), correct_values.end());

  // Test each RankArray separately.
  std::vector<edge_type>::iterator array_iter = correct_values.begin();
  for(size_type i = 0; i < buffers.ra.size(); i++)
  {
    bool early_end = false;
    size_type wrong_values = 0;
    ProducerBuffer<RankArray> buffer(*(buffers.ra[i]));
    while(array_iter != correct_values.end() && array_iter->first <= node_ranges[i].second)
    {
      if(buffer.end()) { early_end = true; break; }
      if(*buffer != *array_iter) { wrong_values++; }
      ++array_iter; ++buffer;
    }
    EXPECT_FALSE(early_end) << "Rank array " << i << " ended early";
    EXPECT_TRUE(buffer.end()) << "Rank array " << i << " is too large";
    EXPECT_EQ(wrong_values, 0u) << "Rank array " << i << ": " << wrong_values << " wrong values";

    // Skip the values that were missing from the file.
    while(array_iter != correct_values.end() && array_iter->first <= node_ranges[i].second) { ++array_iter; }
  }
}

//------------------------------------------------------------------------------

} // namespace
