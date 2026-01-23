#!/bin/bash
cd src
for file in *.cpp; do
    g++ -fopenmp  -Wall -Wextra -Werror -O3 -c "$file" -o "../build/$(basename "${file%.cpp}.out")"
done
cd utils
for file in *.cpp; do
    g++ -fopenmp -Wall -Wextra -Werror -O3 -c "$file" -o "../../build/$(basename "${file%.cpp}.out")"
done

cd ../../build

# Link all object files
echo "Linking object files..."
g++ *.out -o DQMC.o

# Clean up object files (optional)
rm -f *.out