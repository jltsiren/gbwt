//g++ gtest-add.cpp -o testAdd -Igtest/include -Lgtest/lib -lgtest -lpthread
#include <gtest/gtest.h>

int myadd(int a, int b) {
    return a + b;
}

TEST(testCase, test1) {
    EXPECT_EQ(myadd(2, 3), 5);
}

