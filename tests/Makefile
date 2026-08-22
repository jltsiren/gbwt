MAIN_DIR=..
SDSL_DIR=../../sdsl-lite
include $(SDSL_DIR)/Make.helper

# The build settings the library was built with. Optional, since `make clean`
# gets here without them.
-include $(MAIN_DIR)/Make.helper

LIBRARY=$(MAIN_DIR)/lib/libgbwt.a

# Directories for dependencies.
INCLUDES=-I$(MAIN_DIR)/include -I$(INC_DIR) $(GBWT_INCLUDES)
LIBS=-L$(LIB_DIR) -lsdsl -lgtest -lgtest_main $(GBWT_LIBS)

# gtest isn't in pkg-config, so it can't be found by default when installed
# via Homebrew.
ifeq ($(shell uname -s), Darwin)
    ifeq ($(shell if [ -e $(HOMEBREW_PREFIX)/lib/libgtest.a ]; then echo 1; else echo 0; fi), 1)
        $(info Found gtest installed by Homebrew and linked to $(HOMEBREW_PREFIX).)
        INCLUDES += -I$(HOMEBREW_PREFIX)/include
        LIBS += -L$(HOMEBREW_PREFIX)/lib
    endif
endif

CXX_FLAGS=$(MY_CXX_FLAGS) $(GBWT_PARALLEL_FLAGS) $(MY_CXX_OPT_FLAGS) $(INCLUDES)

HEADERS=$(wildcard $(MAIN_DIR)/include/gbwt/*.h)
PROGRAMS=test_utils test_support test_metadata test_gbwt test_bwtmerge test_queries test_variants

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
