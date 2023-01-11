SDSL_DIR=../sdsl-lite
include $(SDSL_DIR)/Make.helper

BUILD_BIN=bin
BUILD_LIB=lib
BUILD_OBJ=obj
BUILD_CUDA_OBJ=obj
SOURCE_DIR=src

# Multithreading with OpenMP.
PARALLEL_FLAGS=-fopenmp -pthread
LIBS=-L$(LIB_DIR) -lsdsl -ldivsufsort -ldivsufsort64

# Apple Clang does not support OpenMP directly, so we need special handling.
ifeq ($(shell uname -s), Darwin)
    # The compiler complains about -fopenmp instead of missing input.
    ifeq ($(strip $(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l)), 1)
        # The compiler only needs to do the preprocessing.
        PARALLEL_FLAGS = -Xpreprocessor -fopenmp -pthread

        # If HOMEBREW_PREFIX is specified, libomp probably cannot be found automatically.
        ifdef HOMEBREW_PREFIX
            PARALLEL_FLAGS += -I$(HOMEBREW_PREFIX)/include
            LIBS += -L$(HOMEBREW_PREFIX)/lib
        # Macports installs libomp to /opt/local/lib/libomp
        else ifeq ($(shell if [ -d /opt/local/lib/libomp ]; then echo 1; else echo 0; fi), 1)
            PARALLEL_FLAGS += -I/opt/local/include/libomp
            LIBS += -L/opt/local/lib/libomp
        endif

        # We also need to link it.
        LIBS += -lomp
    endif
endif

MY_CXX_FLAGS=-std=c++17 -funroll-loops -O3
CXX_FLAGS=$(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(MY_CXX_OPT_FLAGS) -Iinclude  -I$(INC_DIR)

NVCC=nvcc
NVCC_FLAGS=-arch=sm_75 -use_fast_math -use_fast_math -Iinclude -I$(INC_DIR) -O3
NVCC_LIBS=
CUDA_LINK_LIBS= -lcudart -L/usr/local/cuda/lib64

HEADERS=$(wildcard include/gbwt/*.h include/*.cuh)
LIBOBJS=$(addprefix $(BUILD_OBJ)/,algorithms.o bwtmerge.o cached_gbwt.o dynamic_gbwt.o fast_locate.o files.o gbwt.o internal.o metadata.o support.o test.o thrust_sort_optimize.o utils.o variants.o)
LIBRARY=$(BUILD_LIB)/libgbwt.a

PROGRAMS=$(addprefix $(BUILD_BIN)/,build_gbwt build_ri merge_gbwt benchmark metadata_tool remove_seq)
OBSOLETE=build_gbwt build_ri merge_gbwt benchmark metadata_tool remove_seq

.PHONY: all clean directories test
all: directories $(LIBRARY) $(PROGRAMS)

directories: $(BUILD_BIN) $(BUILD_LIB) $(BUILD_OBJ)

$(BUILD_BIN):
	mkdir -p $@

$(BUILD_LIB):
	mkdir -p $@

$(BUILD_OBJ):
	mkdir -p $@

$(BUILD_CUDA_OBJ)/%.o:$(SOURCE_DIR)/%.cu $(HEADERS)
	$(NVCC) $(NVCC_FLAGS) -c $< -o $@ $(NVCC_LIBS)

$(BUILD_OBJ)/%.o:$(SOURCE_DIR)/%.cpp $(HEADERS)
	$(MY_CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) $(CUDA_LINK_LIBS) -c -o $@ $<

$(LIBRARY):$(LIBOBJS)
	ar rcs $@ $(LIBOBJS)

$(BUILD_BIN)/%:$(BUILD_OBJ)/%.o $(LIBRARY)
	$(MY_CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS) $(CUDA_LINK_LIBS)

test:$(LIBRARY)
	cd tests && $(MAKE) test

clean:
	rm -rf $(BUILD_BIN) $(BUILD_LIB) $(BUILD_OBJ)
	rm -f *.o *.a $(OBSOLETE)
	cd tests && $(MAKE) clean
