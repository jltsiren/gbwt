#include <gbwt/bwtmerge.h>
#include <gbwt/dynamic_gbwt.h>

#include <thrust_sort.cuh>

#include <gtest/gtest.h>

#include <iostream>
#include <random>
#include <vector>

using namespace gbwt;

namespace {
//------------------------------------------------------------------------------

int
getMax(int array[], int n)
{
  int max = array[0];
  for (int i = 1; i < n; i++)
    if (array[i] > max)
      max = array[i];
  return max;
}

void
countingSort(int array[], int size, int place)
{
  const int max = 10;
  int output[size];
  int count[max];

  for (int i = 0; i < max; ++i)
    count[i] = 0;
  for (int i = 0; i < size; i++)
    count[(array[i] / place) % 10]++;
  for (int i = 1; i < max; i++)
    count[i] += count[i - 1];
  for (int i = size - 1; i >= 0; i--) {
    output[count[(array[i] / place) % 10] - 1] = array[i];
    count[(array[i] / place) % 10]--;
  }

  for (int i = 0; i < size; i++)
    array[i] = output[i];
}

void
radixSort(int array[], int size)
{
  int max = getMax(array, size);

  for (int place = 1; max / place > 0; place *= 10) {
    countingSort(array, size, place);
  }
}
void
printSequences(std::vector<std::vector<int>>& seqs)
{
  for (size_t i = 0; i < seqs.size(); i++) {
    for (size_t j = 0; j < seqs[i].size(); j++) {
      std::cout << seqs[i][j] << " ";
    }
    std::cout << "\n";
  }
}
void
printArray(std::vector<int>& array)
{
  for (size_t i = 0; i < array.size(); i++) {
    std::cout << array[i] << " ";
  }
  std::cout << "\n";
}
void
sortSequences(std::vector<std::vector<int>>& seqs,
              std::vector<std::vector<int>>& sortedSeqs)
{
  int pos = 1;
  bool isAllEmpty = false;
  std::vector<int> seqLastRank(seqs.size());
  std::iota(std::begin(seqLastRank), std::end(seqLastRank), 0);
  std::vector<int> seqCurrRank(seqLastRank);
  while (!isAllEmpty) {
    isAllEmpty = true;
    std::vector<int> sortedPosition;
    for (int& i : seqLastRank) {
      if (pos == seqs[i][0]) {
        // update seq id to each position
        sortedPosition.push_back(i);
        // update seqCurrRank
        std::vector<int>::iterator pos_seq_id =
          std::find(seqCurrRank.begin(), seqCurrRank.end(), i);
        if (pos_seq_id != seqCurrRank.end()) {
          seqCurrRank.erase(pos_seq_id);
          seqCurrRank.push_back(i);
        }
        // pop out visited
        seqs[i].erase(seqs[i].begin());
      }
      if (!seqs[i].empty()) {
        isAllEmpty = isAllEmpty && seqs[i].empty();
      }
    }
    seqLastRank.swap(seqCurrRank);
    sortedSeqs.push_back(sortedPosition);
    pos++;
  }
}
text_type
getTextBuffer(const std::vector<vector_type>& sequences)
{
  size_type node_width = MILLION;
  text_type input_buffer(64, 0, node_width);
  size_type input_tail = 0;
  for (auto sequence : sequences) {
    for (auto node : sequence) {
      input_buffer[input_tail] = node;
      input_tail++;
    }
    input_buffer[input_tail] = ENDMARKER;
    input_tail++;
  }
  return input_buffer;
}

std::vector<Sequence>
getVectorOfSequences(text_type& text)
{
  std::vector<Sequence> seqs;

  size_type text_length = 17;
  bool seq_start = true;
  int index_sequences = 0;
  for (size_type i = 0; i < text_length; i++) {
    if (seq_start) {
      seqs.push_back(Sequence(text, i, index_sequences));
      seq_start = false;
      index_sequences++;
    }
    if (text[i] == ENDMARKER) {
      seq_start = true;
    }
  }
  return seqs;
}

//------------------------------------------------------------------------------
TEST(RadixSortTest, Scalar)
{
  int array[] = { 121, 432, 564, 23, 1, 45, 788 };
  int expected_array[] = { 1, 23, 45, 121, 432, 564, 788 };
  size_t n = sizeof(array) / sizeof(array[0]);

  radixSort(array, n);

  for (size_t i = 0; i < n; i++) {
    EXPECT_EQ(array[i], expected_array[i]) << "Wrong value at offset " << i;
  }
}
TEST(RadixSortTest, Sequences)
{
  std::vector<std::vector<int>> seqs{ { 1, 2, 4, 6, 7 },
                                      { 1, 2, 5, 7 },
                                      { 1, 3, 4, 5, 7 } };

  std::vector<std::vector<int>> sortedSeqs;
  sortSequences(seqs, sortedSeqs);

  std::vector<std::vector<int>> ans{ { 0, 1, 2 }, { 0, 1 }, { 2 },
                                     { 0, 2 },    { 1, 2 }, { 0 },
                                     { 1, 2, 0 } };

  EXPECT_EQ(sortedSeqs.size(), ans.size())
    << "Sequences size are not consistent " << sortedSeqs.size() << "and"
    << ans.size();
  for (size_t i = 0; i < sortedSeqs.size(); i++) {
    EXPECT_EQ(sortedSeqs[i].size(), ans[i].size())
      << "Sizes are not consistent " << sortedSeqs[i].size() << "and"
      << ans[i].size() << "at offset " << i;
    for (size_t j = 0; j < sortedSeqs[i].size(); j++) {
      EXPECT_EQ(sortedSeqs[i][j], ans[i][j])
        << "Wrong value at offset " << i << "," << j;
    }
  }
}
TEST(RadixSortTest, SerialRadixSort_1)
{
  // Test Case shown on paper
  std::vector<vector_type> test_seqs{ { 1, 2, 4, 6, 7 },
                                      { 1, 2, 5, 7 },
                                      { 1, 3, 4, 5, 7 } };
  text_type text_buffer = getTextBuffer(test_seqs);
  std::vector<Sequence> vec_seqs = getVectorOfSequences(text_buffer);

  std::vector<std::vector<std::pair<size_type, node_type>>> sorted_seqs;

  // serial version
  sortAllSequencesAllPosition(vec_seqs, sorted_seqs, text_buffer);

  std::vector<std::vector<std::pair<size_type, node_type>>> ans{
    { { 0, 2 }, { 1, 2 }, { 2, 3 } }, // Node 1
    { { 0, 4 }, { 1, 5 } },
    { { 2, 4 } },
    { { 0, 6 }, { 2, 5 } },
    { { 1, 7 }, { 2, 7 } },
    { { 0, 7 } },
    { { 1, 0 }, { 2, 0 }, { 0, 0 } }, // Node 7
  };

  EXPECT_EQ(sorted_seqs.size(), ans.size())
    << "Sequences size are not consistent " << sorted_seqs.size() << "and"
    << ans.size();
  for (size_t i = 0; i < sorted_seqs.size(); i++) {
    EXPECT_EQ(sorted_seqs[i].size(), ans[i].size())
      << "Sizes are not consistent " << sorted_seqs[i].size() << "and"
      << ans[i].size() << "at offset " << i;
    for (size_t j = 0; j < sorted_seqs[i].size(); j++) {
      EXPECT_EQ(sorted_seqs[i][j], ans[i][j])
        << "Wrong value at offset " << i << "," << j;
    }
  }
}

TEST(RadixSortTest, SerialRadixSort_2)
{
  // Test Case shown on paper
  std::vector<vector_type> test_seqs{ { 1, 2, 3, 5 },
                                      { 1, 3, 4, 5 },
                                      { 1, 2, 4, 5 } };
  text_type text_buffer = getTextBuffer(test_seqs);
  std::vector<Sequence> vec_seqs = getVectorOfSequences(text_buffer);

  std::vector<std::vector<std::pair<size_type, node_type>>> sorted_seqs;

  // serial version
  sortAllSequencesAllPosition(vec_seqs, sorted_seqs, text_buffer);

  std::vector<std::vector<std::pair<size_type, node_type>>> ans{
    { { 0, 2 }, { 1, 3 }, { 2, 2 } }, // Node 1
    { { 0, 3 }, { 2, 4 } },
    { { 1, 4 }, { 0, 5} },
    { { 2, 5 }, { 1, 5 } },
    { { 0, 0 }, { 2, 0 }, {1, 0} },
  };

  EXPECT_EQ(sorted_seqs.size(), ans.size())
    << "Sequences size are not consistent " << sorted_seqs.size() << "and"
    << ans.size();
  for (size_t i = 0; i < sorted_seqs.size(); i++) {
    EXPECT_EQ(sorted_seqs[i].size(), ans[i].size())
      << "Sizes are not consistent " << sorted_seqs[i].size() << "and"
      << ans[i].size() << "at offset " << i;
    for (size_t j = 0; j < sorted_seqs[i].size(); j++) {
      EXPECT_EQ(sorted_seqs[i][j], ans[i][j])
        << "Wrong value at offset " << i << "," << j;
    }
  }
}
TEST(RadixSortTest, ThrustRadixSort)
{
  // Test Case shown on paper
  std::vector<vector_type> test_seqs{ { 1, 2, 4, 6, 7 },
                                      { 1, 2, 5, 7 },
                                      { 1, 3, 4, 5, 7 } };
  text_type text_buffer = getTextBuffer(test_seqs);
  std::vector<Sequence> vec_seqs = getVectorOfSequences(text_buffer);

  std::vector<size_type> start_pos;
  for (auto& sequence : vec_seqs) {
    start_pos.emplace_back(sequence.pos);
  }
  // cuda version
  auto sorted_seqs = radix_sort(text_buffer, start_pos, 7);

  std::vector<std::vector<std::pair<size_type, node_type>>> ans{
    { { 0, 2 }, { 1, 2 }, { 2, 3 } }, // Node 1
    { { 0, 4 }, { 1, 5 } },
    { { 2, 4 } },
    { { 0, 6 }, { 2, 5 } },
    { { 1, 7 }, { 2, 7 } },
    { { 0, 7 } },
    { { 1, 0 }, { 2, 0 }, { 0, 0 } }, // Node 7
  };

  EXPECT_EQ(sorted_seqs.size(), ans.size())
    << "Sequences size are not consistent " << sorted_seqs.size() << "and"
    << ans.size();
  for (size_t i = 0; i < sorted_seqs.size(); i++) {
    EXPECT_EQ(sorted_seqs[i].size(), ans[i].size())
      << "Sizes are not consistent " << sorted_seqs[i].size() << "and"
      << ans[i].size() << "at offset " << i;
    for (size_t j = 0; j < sorted_seqs[i].size(); j++) {
      EXPECT_EQ(sorted_seqs[i][j], ans[i][j])
        << "Wrong value at offset " << i << "," << j;
    }
  }
}
//-----------------------------------------------------------------

} // namespace
