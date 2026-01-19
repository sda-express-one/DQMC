#!/bin/bash
# compile script for the project

SOURCE_FILE=("/src/main.cpp /src/Diagram.cpp /src/GreenFuncNph.cpp /src/GreenFuncNphBands.cpp /src/utils/MC_Benchmarking.cpp 
    /src/utils/progressbar.cpp")
OUTPUT_FILE="build/main.o"
OPTIMIZATION_LEVEL="-O3"
WARNINGS="-Wall -Wextra -Werror"

# show usage
usage() {
    echo " -s <source_file> Source file(s) to compile (default: $SOURCE_FILE)"
    echo " -o <output_file> Output executable file name"
    echo " -O0, -O1, -O2 Optimization level (default: $OPTIMIZATION_LEVEL)"
    echo " -d Enable debug mode (adds -g flag)"
    echo " -h  Show this help message"
    exit 0
}

while getopts "s:o:O:O0O1O2O3dh" opt; do
    case $opt in
        s) SOURCE_FILE="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        O) OPTIMIZATION_LEVEL="-O$OPTARG" ;;
        O0) OPTIMIZATION_LEVEL="-O0" ;;
        O1) OPTIMIZATION_LEVEL="-O1" ;;
        O2) OPTIMIZATION_LEVEL="-O2" ;;
        O3) OPTIMIZATION_LEVEL="-O3" ;;
        d) DEBUG_FLAGS="-g" ;;
        h) usage ;;
        \?) echo "Invalid option: -$OPTARG" >&2
            exit 1 ;;
    esac
done

# Check if source file exists
if [ ! -f "$SOURCE_FILE" ]; then
    echo "Error: Source file '$SOURCE_FILE' not found!"
    exit 1
fi

# Execute compilation
g++ $WARNINGS $OPTIMIZATION_LEVEL $DEBUG_FLAGS $SOURCE_FILE -o $OUTPUT_FILE

# Check if compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful! Output: $OUTPUT_FILE"
    echo "Run with: ./$OUTPUT_FILE"
else
    echo "Compilation failed!"
    exit 1
fi
