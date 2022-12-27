#SDSL_DIR=../../sdsl-lite
include $(SDSL_DIR)/Make.helper

GBWT_DIR=..
LIBRARY=$(GBWT_DIR)/lib/libgbwt.a

MY_CXX_FLAGS=-std=c++17
CXX_FLAGS=$(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(MY_CXX_OPT_FLAGS) -I$(GBWT_DIR)/include -I$(INC_DIR)

HEADERS=$(wildcard $(GBWT_DIR)/include/gbwt/*.h)
PROGRAMS=test_build

# Multithreading with OpenMP.
PARALLEL_FLAGS=-fopenmp -pthread
LIBS=-L$(LIB_DIR) -lsdsl -ldivsufsort -ldivsufsort64 -lgtest -lgtest_main

.PHONY: all clean test
all:$(PROGRAMS)

%.o:%.cpp $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -c $<

$(PROGRAMS):%:%.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

test:$(PROGRAMS)
	$(foreach PROG,$(PROGRAMS),./$(PROG) || exit 1;)

clean:
	rm -f $(PROGRAMS) *.o
