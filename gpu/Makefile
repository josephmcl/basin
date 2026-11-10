# CUDA port of basin.
#
#   module load gcc/13.1.0 cuda/12.4.1
#   make              # build everything
#   make test         # build and run the tests (needs a GPU node)
#   make clean
#
# Library code lives in source/ (*.cu and *.cpp, all compiled with nvcc);
# each test/test_*.cu becomes bin/test_*.

NVCC      ?= nvcc
ARCH      ?= sm_80
HOST_CXX  ?= g++

NVCCFLAGS ?= -O2 -g -lineinfo
NVCCFLAGS += -std=c++20 -arch=$(ARCH) -ccbin $(HOST_CXX) -Iinclude \
             -Xcompiler -Wall,-Wno-unknown-pragmas --expt-relaxed-constexpr

SRC_DIR  := source
TEST_DIR := test
OBJ_DIR  := build
BIN_DIR  := bin

LIB_SRCS := $(wildcard $(SRC_DIR)/*.cu) $(wildcard $(SRC_DIR)/*.cpp)
LIB_OBJS := $(patsubst $(SRC_DIR)/%,$(OBJ_DIR)/%.o,$(LIB_SRCS))

TEST_SRCS := $(wildcard $(TEST_DIR)/test_*.cu)
TEST_BINS := $(patsubst $(TEST_DIR)/%.cu,$(BIN_DIR)/%,$(TEST_SRCS))

HEADERS := $(wildcard include/*.h) $(wildcard include/*.cuh)

.PHONY: all test clean
all: $(TEST_BINS)

# Objects are only reached through pattern rules, so make would otherwise
# treat them as intermediate files and delete them after every link.
.SECONDARY: $(LIB_OBJS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/% $(HEADERS) | $(OBJ_DIR)
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

$(BIN_DIR)/test_%: $(TEST_DIR)/test_%.cu $(LIB_OBJS) $(HEADERS) | $(BIN_DIR)
	$(NVCC) $(NVCCFLAGS) $< $(LIB_OBJS) -o $@

$(OBJ_DIR) $(BIN_DIR):
	mkdir -p $@

# Tests read ../operator/, so run them from this directory.
test: $(TEST_BINS)
	@for t in $(TEST_BINS); do echo "== $$t"; ./$$t || exit 1; done

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)
