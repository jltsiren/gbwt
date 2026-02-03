SDSL_DIR=../../sdsl-lite
include $(SDSL_DIR)/Make.helper

MAIN_DIR=..
LIBRARY=$(MAIN_DIR)/lib/libgbwt.a

# Multithreading with OpenMP.
PARALLEL_FLAGS=-fopenmp -pthread

# Directories for dependencies.
INCLUDES=-I$(MAIN_DIR)/include -I$(INC_DIR)
LIBS=-L$(LIB_DIR) -lgbwt -lhandlegraph -lsdsl -ldivsufsort -ldivsufsort64 -lgtest -lgtest_main

# Use pkg-config to find system dependencies.
ifeq ($(shell pkg-config --exists libzstd && echo 1), 1)
    $(info Found libzstd.)
    INCLUDES += $(shell pkg-config --cflags libzstd)
    LIBS += $(shell pkg-config --libs libzstd)
else
    $(error Could not find libzstd. Please update PKG_CONFIG_PATH to include zstd development files)
endif

# Apple Clang does not support OpenMP directly, so we need special handling.
ifeq ($(shell uname -s), Darwin)
    # The compiler complains about -fopenmp instead of missing input.
    ifeq ($(strip $(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l)), 1)
        $(info The compiler is Apple Clang that needs libomp for OpenMP support.)

        # The compiler only needs to do the preprocessing.
        PARALLEL_FLAGS=-Xpreprocessor -fopenmp -pthread

        # Find libomp installed by Homebrew or MacPorts.
        ifeq ($(shell if [ -e $(HOMEBREW_PREFIX)/include/omp.h ]; then echo 1; else echo 0; fi), 1)
            $(info Found libomp installed by Homebrew and linked to $(HOMEBREW_PREFIX).)
            INCLUDES += -I$(HOMEBREW_PREFIX)/include
            LIBS += -L$(HOMEBREW_PREFIX)/lib
        else ifeq ($(shell if [ -d $(HOMEBREW_PREFIX)/opt/libomp/include ]; then echo 1; else echo 0; fi), 1)
            $(info Found a keg-only libomp installed by Homebrew at $(HOMEBREW_PREFIX)/opt/libomp.)
            INCLUDES += -I$(HOMEBREW_PREFIX)/opt/libomp/include
            LIBS += -L$(HOMEBREW_PREFIX)/opt/libomp/lib
        else ifeq ($(shell if [ -d /opt/local/lib/libomp ]; then echo 1; else echo 0; fi), 1)
            $(info Found libomp installed by MacPorts at /opt/local.)
            INCLUDES += -I/opt/local/include/libomp
            LIBS += -L/opt/local/lib/libomp
        else
            $(error Could not find libomp. Please install it using Homebrew or MacPorts)
        endif

        # We also need to link it.
        LIBS += -lomp
    endif

    # gtest isn't in pkg-config but also we can't find its libs by default if
    # they are installed via Homebrew.
    ifeq ($(shell if [ -e $(HOMEBREW_PREFIX)/lib/libgtest.a ]; then echo 1; else echo 0; fi), 1)
        $(info Found gtest installed by Homebrew and linked to $(HOMEBREW_PREFIX).)
        INCLUDES += -I$(HOMEBREW_PREFIX)/include
        LIBS += -L$(HOMEBREW_PREFIX)/lib
    endif
endif

CXX_FLAGS=$(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(MY_CXX_OPT_FLAGS) $(INCLUDES)

HEADERS=$(wildcard $(MAIN_DIR)/include/gbwt/*.h)
PROGRAMS=test_utils test_support test_bwtmerge test_metadata test_queries test_variants

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
