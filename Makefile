# Basin: Hybrid SBP-SAT Poisson Solver
#
# Targets:
#   make mkl      - Build MKL solver (main)
#   make generate - Build component generator (precompute M, LambdaA)
#   make gpu      - Build GPU version (future)
#   make clean    - Remove build artifacts

# Compiler configuration
CXX      ?= mpicxx
ICPX     ?= icpx
CXXFLAGS := -std=c++20 -Wall -Wpedantic -O3

# Intel oneAPI paths (override via environment or command line)
ONEAPI_ROOT ?= /home/jmclaug2/intel/oneapi/2025.1
MKL_ROOT    ?= $(ONEAPI_ROOT)/mkl/latest

# Directory structure
SRC_DIR    := source
INC_DIR    := include
OBJ_DIR    := build
BIN_DIR    := bin

# =============================================================================
# MKL Backend
# =============================================================================

MKL_SRC_DIR := $(SRC_DIR)/mkl
MKL_INC_DIR := $(INC_DIR)/mkl
MKL_OBJ_DIR := $(OBJ_DIR)/mkl

# Main solver
MKL_TARGET  := $(BIN_DIR)/basin-mkl
MKL_MAIN    := hybrid_sbp_sat_2d_poisson

# Component generator
GEN_TARGET  := $(BIN_DIR)/basin-generate
GEN_MAIN    := generate_components

MKL_CXX      := $(ONEAPI_ROOT)/bin/mpicxx
MKL_CXXFLAGS := $(CXXFLAGS) -cxx=$(ONEAPI_ROOT)/bin/icpx -qmkl -qopenmp -vec -mavx -fast

MKL_INCLUDES := -I$(MKL_INC_DIR) -I$(MKL_ROOT)/include
MKL_LIBS     := -L$(ONEAPI_ROOT)/lib \
                -lmkl_scalapack_lp64 \
                -lmkl_blacs_intelmpi_lp64 \
                -lmkl_intel_lp64 \
                -lmkl_core \
                -lmkl_intel_thread \
                -liomp5 \
                -lpthread -lm -ldl

# All sources except entry points
MKL_SOURCES := $(filter-out $(MKL_SRC_DIR)/$(MKL_MAIN).cpp $(MKL_SRC_DIR)/$(GEN_MAIN).cpp, \
               $(wildcard $(MKL_SRC_DIR)/*.cpp))
MKL_OBJECTS := $(MKL_SOURCES:$(MKL_SRC_DIR)/%.cpp=$(MKL_OBJ_DIR)/%.o)
MKL_HEADERS := $(wildcard $(MKL_INC_DIR)/*.h)

# Entry point objects
MKL_MAIN_OBJ := $(MKL_OBJ_DIR)/$(MKL_MAIN).o
GEN_MAIN_OBJ := $(MKL_OBJ_DIR)/$(GEN_MAIN).o

# =============================================================================
# GPU Backend (placeholder for future)
# =============================================================================

GPU_SRC_DIR := $(SRC_DIR)/gpu
GPU_INC_DIR := $(INC_DIR)/gpu
GPU_OBJ_DIR := $(OBJ_DIR)/gpu
GPU_TARGET  := $(BIN_DIR)/basin-gpu

# GPU_SOURCES := $(wildcard $(GPU_SRC_DIR)/*.cpp)
# GPU_OBJECTS := $(GPU_SOURCES:$(GPU_SRC_DIR)/%.cpp=$(GPU_OBJ_DIR)/%.o)

# =============================================================================
# Build Rules
# =============================================================================

.PHONY: mkl generate gpu clean help all

help:
	@echo "Basin: Hybrid SBP-SAT Poisson Solver"
	@echo ""
	@echo "Usage:"
	@echo "  make mkl      Build MKL solver"
	@echo "  make generate Build component generator (precompute M, LambdaA)"
	@echo "  make gpu      Build GPU version (not yet implemented)"
	@echo "  make clean    Remove build artifacts"
	@echo ""
	@echo "Configuration (override via environment):"
	@echo "  ONEAPI_ROOT  Intel oneAPI root (current: $(ONEAPI_ROOT))"
	@echo "  MKL_ROOT     MKL root directory (current: $(MKL_ROOT))"

all: mkl generate

# MKL solver target
mkl: $(MKL_TARGET)

$(MKL_TARGET): $(MKL_OBJECTS) $(MKL_MAIN_OBJ) | $(BIN_DIR)
	@echo "[LINK] $@"
	@$(MKL_CXX) $(MKL_CXXFLAGS) $(MKL_OBJECTS) $(MKL_MAIN_OBJ) -o $@ $(MKL_LIBS)

# Component generator target
generate: $(GEN_TARGET)

$(GEN_TARGET): $(MKL_OBJECTS) $(GEN_MAIN_OBJ) | $(BIN_DIR)
	@echo "[LINK] $@"
	@$(MKL_CXX) $(MKL_CXXFLAGS) $(MKL_OBJECTS) $(GEN_MAIN_OBJ) -o $@ $(MKL_LIBS)

# Object file compilation
$(MKL_OBJ_DIR)/%.o: $(MKL_SRC_DIR)/%.cpp $(MKL_HEADERS) | $(MKL_OBJ_DIR)
	@echo "[CXX]  $<"
	@$(MKL_CXX) $(MKL_CXXFLAGS) $(MKL_INCLUDES) -c $< -o $@

# GPU target (placeholder)
gpu:
	@echo "GPU backend not yet implemented"
	@echo "Create source/gpu/ and include/gpu/ to add GPU support"

# Directory creation
$(BIN_DIR):
	@mkdir -p $@

$(MKL_OBJ_DIR):
	@mkdir -p $@

$(GPU_OBJ_DIR):
	@mkdir -p $@

# Clean
clean:
	@echo "[CLEAN] Removing build artifacts"
	@rm -rf $(OBJ_DIR)
	@rm -rf $(BIN_DIR)

# Convenience: rebuild from scratch
rebuild: clean mkl
