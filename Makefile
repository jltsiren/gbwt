SDSL_DIR=../sdsl-lite
include $(SDSL_DIR)/Make.helper

BUILD_BIN=bin
BUILD_LIB=lib
BUILD_OBJ=obj
SOURCE_DIR=src

# --- Optional build-time behaviors, overridable on the command line (e.g.
# `make GBWT_USE_EXCEPTIONS=0`). Defaults match GBWT's traditional behavior.
# See include/gbwt/error_handling.h and include/gbwt/support.h (StringArray)
# for what each one changes.

# Report fatal errors by throwing C++ exceptions (1) or not (0). Building
# with 0 is needed to embed GBWT in exceptions-free code, such as Google's
# DeepVariant, which builds with Bazel and -fno-exceptions.
GBWT_USE_EXCEPTIONS=1
# When GBWT_USE_EXCEPTIONS=0, report fatal errors via Abseil's
# ABSL_LOG(FATAL) (1) instead of stderr + std::abort() (0). No effect when
# GBWT_USE_EXCEPTIONS=1.
GBWT_USE_ABSEIL_LOGGING=0
# Build with OpenMP parallelism (1) or without it (0).
GBWT_USE_OPENMP=1
# Enable StringArray's shared-memory-backed storage, which needs
# Boost.Interprocess and sdsl-lite's own SDSL_ENABLE_SHARED_MEMORY int_vector
# constructor (see SDSL_ENABLE_SHARED_MEMORY in $(SDSL_DIR)/Make.helper,
# which sdsl-lite's CMake build generates). Defaults to on if both are
# available, off otherwise; set explicitly to override the autodetection
# either way.
GBWT_ENABLE_SHARED_MEMORY=$(shell if echo '\#include <boost/interprocess/managed_shared_memory.hpp>' | $(MY_CXX) -std=c++17 -E -x c++ - >/dev/null 2>&1 && [ "$(SDSL_ENABLE_SHARED_MEMORY)" = "ON" ]; then echo 1; else echo 0; fi)

ifeq ($(GBWT_ENABLE_SHARED_MEMORY), 1)
    ifneq ($(SDSL_ENABLE_SHARED_MEMORY), ON)
        $(error GBWT_ENABLE_SHARED_MEMORY=1 needs sdsl-lite built with -DSDSL_ENABLE_SHARED_MEMORY=ON (see SDSL_ENABLE_SHARED_MEMORY in $(SDSL_DIR)/Make.helper); rebuild sdsl-lite with that option, or unset GBWT_ENABLE_SHARED_MEMORY)
    endif
endif

DEFINES=
ifeq ($(GBWT_USE_EXCEPTIONS), 0)
    DEFINES += -DGBWT_NO_EXCEPTIONS
    ifeq ($(GBWT_USE_ABSEIL_LOGGING), 1)
        DEFINES += -DGBWT_USE_ABSEIL_LOGGING
    endif
endif

# Multithreading with OpenMP.
ifeq ($(GBWT_USE_OPENMP), 1)
    PARALLEL_FLAGS=-fopenmp -pthread -DGBWT_USE_OPENMP
else
    # Without OpenMP, the "#pragma omp ..." directives sprinkled through the
    # source are inert (silently ignored by the compiler), but -Wall turns on
    # -Wunknown-pragmas, which would otherwise warn about every one of them.
    PARALLEL_FLAGS=-pthread -Wno-unknown-pragmas
endif

# Directories for dependencies.
INCLUDES=-Iinclude -I$(INC_DIR)
LIBS=-L$(LIB_DIR) -lsdsl

# Use pkg-config to find system dependencies.
ifeq ($(shell pkg-config --exists libzstd && echo 1), 1)
    $(info Found libzstd.)
    INCLUDES += $(shell pkg-config --cflags libzstd)
    LIBS += $(shell pkg-config --libs libzstd)
else
    $(error Could not find libzstd. Please update PKG_CONFIG_PATH to include zstd development files)
endif

ifeq ($(GBWT_USE_ABSEIL_LOGGING), 1)
    ifeq ($(shell pkg-config --exists absl_log absl_log_internal_message && echo 1), 1)
        $(info Found Abseil logging.)
        INCLUDES += $(shell pkg-config --cflags absl_log absl_log_internal_message)
        LIBS += $(shell pkg-config --libs absl_log absl_log_internal_message)
    else
        $(error GBWT_USE_ABSEIL_LOGGING=1 but could not find Abseil logging via pkg-config. Please update PKG_CONFIG_PATH)
    endif
endif

ifeq ($(GBWT_ENABLE_SHARED_MEMORY), 1)
    # GBWT_ENABLE_SHARED_MEMORY and SDSL_ENABLE_SHARED_MEMORY are not passed
    # as -D flags: see the CONFIG_HEADER rule below for why.
    # Boost.Interprocess is header-only, but its named shared memory objects
    # need librt on Linux (not on macOS, where the equivalent is in libc).
    ifneq ($(shell uname -s), Darwin)
        LIBS += -lrt
    endif
endif

# Apple Clang does not support OpenMP directly, so we need special handling.
ifeq ($(shell uname -s), Darwin)
    ifeq ($(GBWT_USE_OPENMP), 1)
    # The compiler complains about -fopenmp instead of missing input.
    ifeq ($(strip $(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l)), 1)
        $(info The compiler is Apple Clang that needs libomp for OpenMP support.)

        # The compiler only needs to do the preprocessing.
        PARALLEL_FLAGS=-Xpreprocessor -fopenmp -pthread -DGBWT_USE_OPENMP

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
    endif
endif

CXX_FLAGS=$(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(DEFINES) $(MY_CXX_OPT_FLAGS) $(INCLUDES)

# A -D flag naming GBWT_ENABLE_SHARED_MEMORY only reaches GBWT's own .cpp
# compiles, not a downstream consumer's (e.g. gbwtgraph's) compiles of these
# same headers, so a consumer would have no reliable way to know whether GBWT
# was actually built with it -- and StringArray's shared-memory-only members
# and constructors (see support.h) differ in shape depending on it. Baking
# the resolved value into an actual header instead lets any translation unit
# that includes GBWT's headers see it, whether or not its own build system
# knows to pass a matching flag.
CONFIG_HEADER=include/gbwt/config.h
HEADERS=$(sort $(wildcard include/gbwt/*.h) $(CONFIG_HEADER))
LIBOBJS=$(addprefix $(BUILD_OBJ)/,algorithms.o bwtmerge.o cached_gbwt.o dynamic_gbwt.o fast_locate.o files.o gbwt.o internal.o metadata.o support.o test.o utils.o variants.o)
LIBRARY=$(BUILD_LIB)/libgbwt.a

PROGRAMS=$(addprefix $(BUILD_BIN)/,build_gbwt build_ri merge_gbwt benchmark metadata_tool remove_seq)
OBSOLETE=build_gbwt build_ri merge_gbwt benchmark metadata_tool remove_seq

.PHONY: all clean directories test .FORCE
.FORCE:

all: directories $(LIBRARY) $(PROGRAMS)

directories: $(BUILD_BIN) $(BUILD_LIB) $(BUILD_OBJ)

$(BUILD_BIN):
	mkdir -p $@

$(BUILD_LIB):
	mkdir -p $@

$(BUILD_OBJ):
	mkdir -p $@

$(BUILD_OBJ)/%.o:$(SOURCE_DIR)/%.cpp $(HEADERS) | $(BUILD_OBJ)
	$(MY_CXX) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -c -o $@ $<

# Unlike Make.helper below, this is only overwritten when its content
# actually changes (via the .tmp/cmp/mv dance), since every compile that
# includes a GBWT header depends on it: touching it on every invocation
# would force a full rebuild of GBWT and everything downstream of it every
# single time, instead of only when a relevant flag actually changed.
$(CONFIG_HEADER): .FORCE
	@mkdir -p $(dir $(CONFIG_HEADER))
	@echo "// Generated by GBWT's Makefile; do not edit, and do not commit to git." > $@.tmp
	@echo "#ifndef GBWT_CONFIG_H" >> $@.tmp
	@echo "#define GBWT_CONFIG_H" >> $@.tmp
	@if [ "$(GBWT_ENABLE_SHARED_MEMORY)" = "1" ]; then \
	  echo "#define GBWT_ENABLE_SHARED_MEMORY" >> $@.tmp; \
	  echo "#define SDSL_ENABLE_SHARED_MEMORY" >> $@.tmp; \
	fi
	@echo "#endif" >> $@.tmp
	@cmp -s $@.tmp $@ || cp $@.tmp $@
	@rm -f $@.tmp

# Generated for tests/Makefile to include, so it builds against the exact
# same defines/flags/libraries this Makefile resolved for $(LIBRARY) --
# including autodetected and command-line-overridden ones -- instead of
# separately re-deriving them and potentially disagreeing with what the
# library was actually built with. Regenerated on every invocation (via the
# .FORCE prerequisite) to stay in sync with whatever this `make` command
# line resolved; not meant to be edited or committed to git.
Make.helper: .FORCE
	@echo "# Generated by GBWT's Makefile; do not edit, and do not commit to git." > $@
	@echo "MY_CXX=$(MY_CXX)" >> $@
	@echo "MY_CC=$(MY_CC)" >> $@
	@echo "MY_CXX_FLAGS=$(MY_CXX_FLAGS)" >> $@
	@echo "MY_CXX_OPT_FLAGS=$(MY_CXX_OPT_FLAGS)" >> $@
	@echo "GBWT_DEFINES=$(DEFINES)" >> $@
	@echo "GBWT_PARALLEL_FLAGS=$(PARALLEL_FLAGS)" >> $@
	@echo "GBWT_INCLUDES=$(INCLUDES)" >> $@
	@echo "GBWT_LIBS=$(LIBS)" >> $@

# Make.helper is an order-only prerequisite: it must exist and be current
# before we link, but its own rule always running shouldn't itself force
# $(LIBRARY) to be considered out of date (and relinked) on every `make`.
$(LIBRARY): $(LIBOBJS) | Make.helper $(BUILD_LIB)
	ar rcs $@ $(LIBOBJS)

$(BUILD_BIN)/%:$(BUILD_OBJ)/%.o $(LIBRARY) | $(BUILD_BIN)
	$(MY_CXX) $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)

test:$(LIBRARY)
	cd tests && $(MAKE) test

clean:
	rm -rf $(BUILD_BIN) $(BUILD_LIB) $(BUILD_OBJ)
	rm -f *.o *.a $(OBSOLETE) Make.helper $(CONFIG_HEADER)
	cd tests && $(MAKE) clean
