# Makefile for building they hybridized SBP-SAT basin experiment. 
# 
#

include .env

# Avoid naming collisions by modifying this variable for different 
# hardware platforms. 
hardware = saturn

target     = basin-petsc
mkl_target = basin-mkl-$(hardware) 

compiler_binary = mpiicpx -march=native -qmkl -fast -fiopenmp -mavx
# compiler_binary = /storage/users/josephmcl/intel/oneapi/compiler/latest/bin/icpx -qmkl -fast -fiopenmp -mavx

source_ext = cpp
header_directory = include
source_directory = source
object_directory = object
binary_directory = .
test_directory   = test

sources := $(wildcard $(source_directory)/*.cpp)
headers := $(wildcard $(header_directory)/*.h)
objects := $(sources:$(source_directory)/%.cpp=$(object_directory)/%.o)
test_sources := $(wildcard $(test_directory)/*.cpp)
test_headers := $(wildcard $(test_directory)/*.h)
test_objects := $(test_sources:$(test_directory)/%.cpp=$(object_directory)/%.o)	
test_target  := test


# MKL implementation
mkl_header_directory := $(header_directory)/mkl
mkl_source_directory := $(source_directory)/mkl
mkl_object_directory := $(object_directory)/mkl

mkl_object_directory_absent = $(mkl_object_directory)-

mkl_sources := $(wildcard $(mkl_source_directory)/*.cpp)
mkl_headers := $(wildcard $(mkl_header_directory)/*.h)
mkl_objects := $(mkl_sources:$(mkl_source_directory)/%.cpp=$(mkl_object_directory)/%.o)

test_entrypoint := test
test_target := lab

nil := 
space := $(nil) $(nil)

# The testing framework needs access to both, the project objects and  
# test framework objects. However, it should not include any other 
# entrypoint objects.
objects_and_test_objects = $(nil)
objects_and_test_objects += $(filter-out \
	$(object_directory)/$(entrypoint).o, $(objects))
objects_and_test_objects += $(test_objects)

rm = rm -f


test_includes := $(includes) -I./$(test_directory)

gccflags = -std=c++2a -Wall -Wpedantic  -O3 -finput-charset=UTF-8
compiler_flags := $(gccflags) -g

mkl_includes  := -I./$(header_directory)/mkl -I/opt/intel/oneapi/vtune/2023.1.0/include 
mkl_libraries := -Wl,--no-as-needed -lmkl_intel_lp64 \
	-lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl -L/opt/intel/oneapi/vtune/2023.1.0/lib64

define speaker
	@echo [make:$$PPID] $(1)
	@$(1)
endef

define compiler
	@echo [make:$$PPID] $(compiler_binary) $(1)
	@$(compiler_binary) $(1)
endef

mkdir = mkdir -p

optz = -xHost -vec  -vecabi=cmdtarget -axAVX -fiopenmp-simd \
	-DKNOWN_TRIP_COUNT
#-qopt-report
# -xCORE-AVX512
#-qopt-report-stdout

$(binary_directory)/$(target): $(objects)
	$(call compile, $(objects) -o $@ $(libraries))

$(binary_directory)/$(test_target): $(objects) $(test_objects)
	@echo $(objects)
	@echo $(test_objects)
	@echo $(objects_and_test_objects)
	$(call compiler, $(objects_and_test_objects) -o $@ $(libraries))

$(objects): $(object_directory)/%.o: $(source_directory)/%.$(source_ext) 
	$(call compiler, $(compiler_flags) -c $< -o $@ $(includes)) 

$(mkl_object_directory_absent):
	$(call speaker, $(mkdir) $(mkl_object_directory))

$(test_objects): $(object_directory)/%.o : $(test_directory)/%.cpp
	$(call compiler, $(compiler_flags) -c $< -o $@ $(test_includes)) 

$(binary_directory)/$(mkl_target): $(mkl_objects) 
	$(call compiler, $(optz) $(mkl_objects) -o $@ $(mkl_libraries))
 
$(mkl_objects): $(mkl_object_directory)/%.o: \
$(mkl_source_directory)/%.$(source_ext) $(mkl_object_directory_absent)
	$(call compiler, $(compiler_flags)  -c $< -o $@ $(mkl_includes)) 

.PHONY: petsc 
petsc: $(binary_directory)/$(target)

.PHONY: mkl 
mkl: $(binary_directory)/$(mkl_target)

.PHONY: clean
clean:
	@$(rm) -rf $(objects) $(mkl_objects)
	@$(rm) -rf $(test_objects)

.PHONY: remove
remove: clean
	@$(rm) $(BINDIR)/$(target)

