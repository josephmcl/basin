# Makefile for building they hybridized SBP-SAT basin experiment. 
# 
#

include .env

target = main
mkl_target = lab-mkl 

# cc = g++-12
# cc = g++-13

cc = g++ -mkl -DMKL

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


# New stuff of MKL implementation
mkl_impl_header_directory := $(header_directory)/mkl
mkl_impl_source_directory := $(source_directory)/mkl
mkl_impl_object_directory := $(object_directory)/mkl

mkl_impl_object_directory_absent = $(mkl_impl_object_directory)-

mkl_impl_sources := $(wildcard $(mkl_impl_source_directory)/*.cpp)
mkl_impl_headers := $(wildcard $(mkl_impl_header_directory)/*.h)
mkl_impl_objects := $(mkl_impl_sources:$(mkl_impl_source_directory)/%.cpp=$(mkl_impl_object_directory)/%.o)

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

gccflags = -std=c++2a -Wall -Wpedantic -g -O3 -finput-charset=UTF-8 

compiler_flags := $(gccflags) 

openblas_include := -I$(OPENBLAS_INCLUDE)
petsc_include    := -I${PETSC_INCLUDE}
openmpi_include  := -I${OPENMPI_INCLUDE}
cernroot_include := -I${CERN_ROOT_INCLUDE}
openmp_include   := -I/opt/homebrew/Cellar/libomp/16.0.4/include 
mkl_include      := -m64 -I${MKLROOT}/include 

includes := -I./$(header_directory)/common -I./$(header_directory)/mkl $(mkl_include) 
#  \
#	$(openmpi_include) 
# $(openmp_include)
# $(cernroot_include)

test_includes := $(includes) -I./$(test_directory)

cernroot_library := -L${CERN_ROOT_LIBRARY} #--with-debugging=yes

home_library := -L/opt/homebrew/lib 
pthread_library := -lpthread
petsc_library := ${PETSC_LIBRARY}
mpi_library := /opt/homebrew/Cellar/open-mpi/4.1.5/lib/libmpi.dylib
openmp_library := /opt/homebrew/Cellar/libomp/16.0.4/lib/libomp.dylib
mkl_library := /gpfs/packages/spack/spack-rhel8/opt/spack/linux-rhel8-broadwell/gcc-13.1.0/intel-oneapi-mkl-2023.1.0-zq7rmpuigpkajdv7mkddhublalbxjbhy/mkl/2023.1.0/lib
#-lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lm 

 # libraries := -L/usr/local/opt/openblas/lib -lopenblas -lpthread
libraries := $(petsc_library) $(pthread_library) $(mkl_library) 
#	$(mpi_library) -fopenmp 
# $(openmp_library)
 # -lomp #$(cernroot_library) 

oneapi_root := /gpfs/packages/spack/spack-rhel8/opt/spack/linux-rhel8-broadwell/gcc-13.1.0/intel-oneapi-compilers-2023.1.0-3d5dbsmapp7perx5ikhy4b2dwpkoiz7w/compiler/2023.1.0/linux/include
oneapi_sycl := $(oneapi_root)/sycl/
oneapi_lib := -L/gpfs/packages/spack/spack-rhel8/opt/spack/linux-rhel8-broadwell/gcc-13.1.0/intel-oneapi-compilers-2023.1.0-3d5dbsmapp7perx5ikhy4b2dwpkoiz7w/compiler/2023.1.0/linux/lib/

oneapi_include := -I$(oneapi_root) -I$(oneapi_sycl)

# More MKL Implementation stuff
mkl_impl_includes  := -I./$(header_directory)/common -I./$(header_directory)/mkl   $(mkl_include) # $(oneapi_include)
mkl_impl_libraries := $(pthread_library) -L$(mkl_library) -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lpthread -lm -ldl -lgomp -fopenmp # $(oneapi_lib) 
define speaker
	@echo [make:$$PPID] $(1)
	@$(1)
endef

mkdir = mkdir -p

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

$(mkl_impl_object_directory_absent):
	$(call speaker, $(mkdir) $(mkl_impl_object_directory))

$(test_objects): $(object_directory)/%.o : $(test_directory)/%.cpp
	$(call speaker,\
	$(cc) $(compiler_flags) -c $< -o $@ $(test_includes)) 

$(binary_directory)/$(mkl_target): $(mkl_impl_objects) 
	$(call speaker,\
	$(cc) $(mkl_impl_objects) -o $@ $(mkl_impl_libraries))
 
$(mkl_impl_objects): $(mkl_impl_object_directory)/%.o: $(mkl_impl_source_directory)/%.$(source_ext) $(mkl_impl_object_directory_absent)
	$(call speaker,\
	$(cc) $(compiler_flags) -c $< -o $@ $(mkl_impl_includes)) 

.PHONY: petsc 
petsc: $(binary_directory)/$(target)

.PHONY: mkl 
mkl: $(binary_directory)/$(mkl_target)

.PHONY: clean
clean:
	@$(rm) -rf $(objects) $(mkl_impl_objects)
	@$(rm) -rf $(test_objects)

.PHONY: remove
remove: clean
	@$(rm) $(BINDIR)/$(target)

