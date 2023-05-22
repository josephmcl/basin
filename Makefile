# Makefile for building they hybridized SBP-SAT basin experiment. 
# 
#

include .env

target = main

# cc = g++-12
cc = g++-13

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


entrypoint      := main
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
#, \
#	$(join $(join $(objects), $(space)), $(test_objects)))

rm = rm -f

gccflags = -std=c++2a -Wall -Wpedantic -g -O3 \
		   -fconcepts-diagnostics-depth=3 -finput-charset=UTF-8 -fopenmp

compiler_flags := $(gccflags) 

openblas_include := -I$(OPENBLAS_INCLUDE)
petsc_include    := -I${PETSC_INCLUDE}
openmpi_include  := -I${OPENMPI_INCLUDE}
cernroot_include := -I${CERN_ROOT_INCLUDE}
openmp_include   := -I/opt/homebrew/Cellar/libomp/16.0.4/include

includes := -I./$(header_directory) $(petsc_include) \
	$(openmpi_include) $(openmp_include)
 # $(cernroot_include)

test_includes := $(includes) -I./$(test_directory)

cernroot_library := -L${CERN_ROOT_LIBRARY} #--with-debugging=yes

home_library := -L/opt/homebrew/lib 
pthread_library := -lpthread
petsc_library := /opt/homebrew/Cellar/petsc/3.19.1/lib/libpetsc.dylib
mpi_library := /opt/homebrew/Cellar/open-mpi/4.1.5/lib/libmpi.dylib
openmp_library := /opt/homebrew/Cellar/libomp/16.0.4/lib/libomp.dylib


 # libraries := -L/usr/local/opt/openblas/lib -lopenblas -lpthread
libraries := $(home_library) $(petsc_library) $(pthread_library) \
	$(mpi_library) $(openmp_library)
 # -lomp #$(cernroot_library) 

define speaker
	@echo [make:$$PPID] $(1)
	@$(1)
endef

$(binary_directory)/$(target): $(objects)
	$(call speaker,\
	$(cc) $(objects) -o $@ $(libraries))

$(binary_directory)/$(test_target): $(objects) $(test_objects)
	@echo $(objects)
	@echo $(test_objects)
	@echo $(objects_and_test_objects)
	$(call speaker,\
	$(cc) $(objects_and_test_objects) -o $@ $(libraries))

$(objects): $(object_directory)/%.o: $(source_directory)/%.$(source_ext)
	$(call speaker,\
	$(cc) $(compiler_flags) -c $< -o $@ $(includes)) 

$(test_objects): $(object_directory)/%.o : $(test_directory)/%.cpp
	$(call speaker,\
	$(cc) $(compiler_flags) -c $< -o $@ $(test_includes)) 


.PHONY: clean
clean:
	@$(rm) $(objects)
	@$(rm) $(test_objects)

.PHONY: remove
remove: clean
	@$(rm) $(BINDIR)/$(target)
