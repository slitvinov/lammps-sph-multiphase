#! /bin/bash

# a first argument is an output variable name
name=$1
shift

# find minimum of arguments
min=$(echo $* | xargs -n1 | sort -g | head -n1)

# store it in a file
echo "variable ${name} equal ${min}" > in.${name}
