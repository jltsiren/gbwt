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

#include <gbwt/bwtmerge.h>

using namespace gbwt;

namespace
{

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

// FIXME
// Test GapArray
// Test GapArray iterators
// Test RankArray

//------------------------------------------------------------------------------

} // namespace
