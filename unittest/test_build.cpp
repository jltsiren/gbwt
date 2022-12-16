#include <iostream>
#include <random>

#include <gbwt/dynamic_gbwt.h>

#include <gtest/gtest.h>

int getMax(int array[], int n) {
	int max = array[0];
	for (int i = 1; i < n; i++)
		if (array[i] > max) max = array[i];
	return max;
}

void countingSort(int array[], int size, int place) {
	const int max = 10;
	int output[size];
	int count[max];

	for (int i = 0; i < max; ++i) count[i] = 0;
	for (int i = 0; i < size; i++) count[(array[i] / place) % 10]++;
	for (int i = 1; i < max; i++) count[i] += count[i - 1];
	for (int i = size - 1; i >= 0; i--) {
		output[count[(array[i] / place) % 10] - 1] = array[i];
		count[(array[i] / place) % 10]--;
	}

	for (int i = 0; i < size; i++) array[i] = output[i];
}

void radixSort(int array[], int size) {
	int max = getMax(array, size);

	for (int place = 1; max / place > 0; place *= 10) {
		countingSort(array, size, place);
	}
}

TEST(RadixSortTest, Scalar) {
	int array[] = {121, 432, 564, 23, 1, 45, 788};
	int expected_array[] = {1, 23, 45, 121, 432, 564, 788};
	size_t n = sizeof(array) / sizeof(array[0]);
	
	radixSort(array, n);
  	
  	for(size_t i = 0; i < n; i++)
  	{
    	EXPECT_EQ(array[i], expected_array[i]) << "Wrong value at offset " << i;
  	}
}

