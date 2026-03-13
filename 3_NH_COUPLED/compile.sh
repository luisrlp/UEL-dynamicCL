#!/bin/bash

# Define the output file for concatenation
CONCAT_FILE="uel.f90"

# Concatenate all source files from the src directory
# cat src/element/*.f90 src/material/*.f90 src/global.f90 > $CONCAT_FILE
# cat src/element/*.f90 src/material/*.f90 > $CONCAT_FILE
cat src/element/*.f90 src/material/*.f90 > $CONCAT_FILE

# Compile global.f to generate global.mod
gfortran -c global.f90

# Compile the concatenated file and main.for
gfortran -o my_program main.f90 $CONCAT_FILE global.o

# Compile concat_file with global.f90 into folder test_in_abaqus
cat global.f90 $CONCAT_FILE > test_in_abaqus/uel.f90

echo "Compilation finished. Executable created: my_program"