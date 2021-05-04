/*
  Copyright (c) 2019, 2021 Jouni Siren

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

class StringArrayTest : public ::testing::Test
{
public:
  void check_array(const StringArray& array, const std::vector<std::string>& truth) const
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
      view_type view = array.view(i);
      std::string from_view(view.first, view.second);
      EXPECT_EQ(from_view, correct) << "Incorrect view of string " << i;
    }
  }

  void try_remove(const std::vector<std::string>& original, size_type i) const
  {
    std::vector<std::string> copy = original;
    if(i < original.size()) { copy.erase(copy.begin() + i); }
    StringArray truth(copy);

    StringArray removed(original);
    removed.remove(i);
    EXPECT_EQ(removed, truth) << "Remove failed for " << i << " / " << original.size();
  }

  void check_file_size(const StringArray& original, std::ifstream& in) const
  {
    size_type expected_size = original.simple_sds_size() * sizeof(sdsl::simple_sds::element_type);
    size_type file_size = fileSize(in);
    ASSERT_EQ(expected_size, file_size) << "Incorrect file size";
  }
};

TEST_F(StringArrayTest, DefaultEmptyArray)
{
  std::vector<std::string> truth;
  StringArray array;
  this->check_array(array, truth);
}

TEST_F(StringArrayTest, FromEmptyArray)
{
  std::vector<std::string> truth;
  StringArray array(truth);
  this->check_array(array, truth);
}

TEST_F(StringArrayTest, NonEmptyArray)
{
  std::vector<std::string> truth
  {
    "first", "second", "third", "fourth"
  };
  StringArray array(truth);
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
  StringArray array(original.size(), [](size_type i) -> bool
  {
    return ((i & 1) != 0);
  }, [&original](size_type i) -> size_type
  {
    return original[i].length();
  }, [&original](size_type i) -> view_type
  {
    return str_to_view(original[i]);
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

TEST_F(StringArrayTest, SerializeEmpty)
{
  std::vector<std::string> truth;
  StringArray original(truth);

  std::string filename = TempFile::getName("string-array");
  std::ofstream out(filename, std::ios_base::binary);
  original.serialize(out);
  out.close();

  StringArray copy;
  std::ifstream in(filename, std::ios_base::binary);
  copy.load(in);
  in.close();
  ASSERT_EQ(copy, original) << "Serialization changed the empty array";

  TempFile::remove(filename);
}

TEST_F(StringArrayTest, CompressEmpty)
{
  std::vector<std::string> truth;
  StringArray original(truth);

  std::string filename = TempFile::getName("string-array");
  std::ofstream out(filename, std::ios_base::binary);
  original.simple_sds_serialize(out);
  out.close();

  StringArray copy;
  std::ifstream in(filename, std::ios_base::binary);
  this->check_file_size(original, in);
  copy.simple_sds_load(in);
  in.close();
  ASSERT_EQ(copy, original) << "Compression changed the empty array";

  TempFile::remove(filename);
}

TEST_F(StringArrayTest, SerializeNonEmpty)
{
  std::vector<std::string> truth
  {
    "first",
    "second",
    "third",
    "fourth"
  };
  StringArray original(truth);

  std::string filename = TempFile::getName("string-array");
  std::ofstream out(filename, std::ios_base::binary);
  original.serialize(out);
  out.close();

  StringArray copy;
  std::ifstream in(filename, std::ios_base::binary);
  copy.load(in);
  in.close();
  ASSERT_EQ(copy, original) << "Serialization changed the non-empty array";

  TempFile::remove(filename);
}

TEST_F(StringArrayTest, CompressNonEmpty)
{
  std::vector<std::string> truth
  {
    "first",
    "second",
    "third",
    "fourth"
  };
  StringArray original(truth);

  std::string filename = TempFile::getName("string-array");
  std::ofstream out(filename, std::ios_base::binary);
  original.simple_sds_serialize(out);
  out.close();

  StringArray copy;
  std::ifstream in(filename, std::ios_base::binary);
  this->check_file_size(original, in);
  copy.simple_sds_load(in);
  in.close();
  ASSERT_EQ(copy, original) << "Compression changed the non-empty array";

  TempFile::remove(filename);
}

TEST_F(StringArrayTest, CompressWithEmptyString)
{
  // Here we test that the compression still works when there is an empty
  // string in the middle and the sd_vector used for the offsets contains
  // duplicate values.
  std::vector<std::string> truth
  {
    "first",
    "second",
    "",
    "fourth"
  };
  StringArray original(truth);

  std::string filename = TempFile::getName("string-array");
  std::ofstream out(filename, std::ios_base::binary);
  original.simple_sds_serialize(out);
  out.close();

  StringArray copy;
  std::ifstream in(filename, std::ios_base::binary);
  this->check_file_size(original, in);
  copy.simple_sds_load(in);
  in.close();
  ASSERT_EQ(copy, original) << "Compression changed the array with an empty string in the middle";

  TempFile::remove(filename);
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

} // namespace
