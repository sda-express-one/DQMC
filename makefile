# default to release mode
BUILD ?= release

CXX = g++
CXXFLAGS = -fopenmp -Wall -Wextra -Werror

# target name and build directory
PROGRAM_NAME = DQMC.o
BUILD_DIR = build/

# source files
SRC_DIR = src
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
SOURCES += $(wildcard $(SRC_DIR)/utils/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

# build-specific flags
DEBUG_FLAGS = -g -O0 -fsanitize=address -fsanitize=undefined -DDEBUG
RELEASE_FLAGS = -O3 -DNDEBUG

# select flags based on BUILD variable
ifeq ($(BUILD),debug)
    CXXFLAGS += $(DEBUG_FLAGS)
	BUILD_DIR = build/$(BUILD)
else
    CXXFLAGS += $(RELEASE_FLAGS)
endif

# final target with build directory
TARGET = $(BUILD_DIR)/$(PROGRAM_NAME)

# default target
all: $(TARGET)

# link the program
$(TARGET): $(OBJECTS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $^ -o $@

# compile source files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# create build directory if it doesn't exist
$(BUILD_DIR): 
		mkdir -p $@

# convenience targets

debug:
	$(MAKE) BUILD=debug

release:
	$(MAKE) BUILD=release

clean:
	rm -f $(OBJECTS)

cleanall:
	rm -f $(OBJECTS) $(TARGET)

.PHONY: all debug release clean cleanall