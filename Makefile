# Makefile for building they hybridized SBP-SAT basin experiment. 
# 
#

target = main

cc = g++-11

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

objects_&_test_objects := $(filter-out $(object_directory)/$(entrypoint).o, \
	$(join $(join$(objects), $(space)), $(test_objects)))

rm = rm -f

gccflags = -std=c++2a -Wall -Wpedantic -g \
		   -fconcepts-diagnostics-depth=1 -finput-charset=UTF-8

compiler_flags := $(gccflags)
includes := -I./$(header_directory) 
test_includes := $(includes) -I./$(test_directory)

define speaker
	@echo "[make]$(1)" 
	@$(1)
endef

$(binary_directory)/$(target): $(objects)
	$(call speaker,\
	$(cc) $(objects) -o $@)

$(binary_directory)/$(test_target): $(objects) $(test_objects)
	$(call speaker,\
	$(cc) $(objects_&_test_objects) -o $@)

$(objects): $(object_directory)/%.o : $(source_directory)/%.cpp
	$(call speaker,\
	$(cc) $(compiler_flags) -c $< -o $@ $(includes))

$(test_objects): $(object_directory)/%.o : $(test_directory)/%.cpp
	$(call speaker,\
	$(cc) $(compiler_flags) -c $< -o $@ $(test_includes))


.PHONY: clean
clean:
	@$(rm) $(objects)

.PHONY: remove
remove: clean
	@$(rm) $(BINDIR)/$(target)
