/*
  Copyright (c) 2019 Jouni Siren

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

#include <gbwt/support.h>

using namespace gbwt;

namespace
{

//------------------------------------------------------------------------------

void
initArray(std::vector<size_type>& array, size_type seed = 0xDEADBEEF)
{
  constexpr size_type TOTAL_VALUES = KILOBYTE;
  constexpr size_type UNIVERSE_SIZE = MEGABYTE;

  array.clear();
  array.reserve(TOTAL_VALUES);
  std::mt19937_64 rng(seed);
  for(size_type i = 0; i < TOTAL_VALUES; i++)
  {
    array.push_back(rng() % UNIVERSE_SIZE);
  }
  sequentialSort(array.begin(), array.end());
}

//------------------------------------------------------------------------------

TEST(SDIteratorTest, Select)
{
  std::vector<size_type> array;
  initArray(array);
  sdsl::sd_vector<> v(array.begin(), array.end());

  // select() returns the original values.
  bool ok = true;
  size_type failed_at = 0, got = 0;
  for(size_type i = 0; i < array.size(); i++)
  {
    SDIterator iter(v, i + 1);
    if(*iter != array[i]) { ok = false; failed_at = i; got = *iter; break; }
  }

  EXPECT_TRUE(ok) << "At " << failed_at << ": expected " << array[failed_at] << ", got " << got;
}

TEST(SDIteratorTest, Iterator)
{
  std::vector<size_type> array;
  initArray(array);
  sdsl::sd_vector<> v(array.begin(), array.end());
  SDIterator iter(v, 1);

  // Right size.
  ASSERT_EQ(iter.size(), array.size()) << "The number of values is wrong";

  // Iterate over the values.
  bool ok = true;
  size_type failed_at = 0, got = 0, rank = 0;
  size_type found_values = 0;
  while(!(iter.end()) && found_values < array.size())
  {
    if(*iter != array[found_values] || iter.rank() != found_values)
    {
      ok = false;
      failed_at = found_values; got = *iter; rank = iter.rank();
      break;
    }
    found_values++;
    ++iter;
  }

  ASSERT_TRUE(ok) << "Expected " << array[failed_at] << " at " << failed_at << ", got " << got << " at " << rank;
  EXPECT_TRUE(iter.end()) << "The iterator finds too many values";
  EXPECT_EQ(found_values, array.size()) << "Expected " << array.size() << found_values << ", got " << found_values;
}

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

//------------------------------------------------------------------------------

} // namespace
