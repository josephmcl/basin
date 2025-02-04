include .env

application = poisson2d
platform    = ${PLATFORM}
hardware    = ${HARDWARE}

rocm_install_directory = ${ROCM_INSTALL_DIRECTORY}

target = $(application).$(platform).$(hardware) 
symlink = ./latest

rocm_compiler   := hipcc
native_compiler := gcc

ifeq ($(platform), rocm)
	compiler_binary = $(rocm_compiler) 
else
	compiler_binary = $(native_compiler)
endif

source_ext = cpp
header_directory = include
source_directory = source
object_directory = object
binary_directory = target
object_directory_absent = $(object_directory)-

sources := $(wildcard $(source_directory)/*.cpp)
headers := $(wildcard $(header_directory)/*.h)
objects := $(sources:$(source_directory)/%.cpp=$(object_directory)/%.o)

nil := 
space := $(nil) $(nil)

# Commands
rm = rm -f
mkdir = mkdir -p
ln = ln -sf


gccf_lags := -std=c++2a -Wall -Wpedantic  -O3 -finput-charset=UTF-8
compiler_flags := $(gccf_lags) -g

includes_local = -I./$(header_directory)
includes_rocblas = -I$(rocm_install_directory)/include/rocblas
includes_rocsolver = -I$(rocm_install_directory)/include/rocsolver
includes_rocalution = -I$(rocm_install_directory)/include/rocalution

includes  := $(includes_local) $(includes_rocalution)
libraries := -L$(rocm_install_directory) -lrocalution -lrocsolver

define speaker
	@echo [make:$$PPID] $(1)
	@$(1)
endef

define compiler
	@echo [make:$$PPID] $(compiler_binary) $(1)
	@$(compiler_binary) $(1)
endef

$(binary_directory)/$(target): $(objects)
	$(call compiler, $(objects) -o $@ $(libraries))

$(objects): $(object_directory)/%.o: \
$(source_directory)/%.$(source_ext) $(object_directory_absent)
	$(call compiler, $(compiler_flags)  -c $< -o $@ $(includes)) 


$(object_directory_absent):
	$(call speaker, $(mkdir) $(object_directory))

.PHONY: basin 
basin: $(binary_directory)/$(target)
	$(call speaker, $(ln) $(binary_directory)/$(target) $(symlink))

.PHONY: clean
clean:
	@$(rm) -rf $(objects) $(objects)

.PHONY: reset
reset: clean
	@$(rm) $(binary_directory)/$(target)