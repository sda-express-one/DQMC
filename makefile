# default to release mode
BUILD ?= release

CXX = g++
CXXFLAGS = -fopenmp -Wall -Wextra -Werror -Iinclude
LDFLAGS = -lyaml-cpp

# directories
SRC_DIR = src
INCLUDE_DIR = include
BUILD_DIR = build/$(BUILD)
BIN_DIR = bin

# target name and build directory
PROGRAM_NAME = DQMC.o
BUILD_DIR = build/

# source files  and corresponding object files
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
SOURCES += $(wildcard $(SRC_DIR)/utils/*.cpp)

# map src/file.cpp -> build/release/foo.o (preserving directory structure)
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SOURCES))

# build-specific flags
# -fsanitize=thread -fsanitize=address
DEBUG_FLAGS = -g -O0 -fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined -DDEBUG
RELEASE_FLAGS = -O3 -DNDEBUG

# select flags based on BUILD variable
ifeq ($(BUILD),debug)
    CXXFLAGS += $(DEBUG_FLAGS)
else
    CXXFLAGS += $(RELEASE_FLAGS)
endif

# final target with build directory
PROGRAM_NAME = DQMC
TARGET = $(BIN_DIR)/$(PROGRAM_NAME)

# default target
all: $(TARGET)

# link the program
$(TARGET): $(OBJECTS) | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

# compile source files (preserving subdirectory structure in build/)
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# create directories if they don't exist
$(BUILD_DIR) $(BIN_DIR): 
		mkdir -p $@

# convenience targets
debug:
	$(MAKE) BUILD=debug

release:
	$(MAKE) BUILD=release

clean:
	rm -f $(OBJECTS)

cleanall:
	rm -rf build/ $(BIN_DIR)

.PHONY: all debug release clean cleanall
