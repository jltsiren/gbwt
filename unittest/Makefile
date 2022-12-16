CXX := g++-11
CXX_FLAGS := -Igtest/include -Lgtest/lib -lgtest -lpthread

TEST_BUILD: test_build.cpp main.cpp 
			$(CXX) test_build.cpp main.cpp -o TEST_BUILD $(CXX_FLAGS)

.PHONY: clean
clean:
		@echo "cleaning"
		rm -rf TEST_BUILD
