#include <gbwt/dynamic_gbwt.h>
#include <gtest/gtest.h>

#include <iostream>
#include <random>
#include <vector>

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
