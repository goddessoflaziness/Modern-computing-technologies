#!/bin/bash

touch result.txt
echo "" > result.txt
for i in {3..10}
do
    input_value=$((2 ** i)) 
    echo "Size of matrix: $input_value" >> result.txt 
    ./build/solver $input_value >> result.txt
done
