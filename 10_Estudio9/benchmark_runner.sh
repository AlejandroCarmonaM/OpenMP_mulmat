#!/bin/bash

# Usage: ./benchmark_runner.sh <matrix_dimensions>
# Possible values for matrix_dimensions: 512, 1024, 2048, all
# Example: ./benchmark_runner.sh 512

# Retrieve input parameter
matrix_dimensions=$1

# Compile TMulMatCuaBlo-bench (assuming "compilar" is a script or command to compile)
./compilar TMulMatCuaBlo-bench

# Write CSV header
echo "dimension_matriz_n,tamano_bloque_b,numero_hilos_t,tiempo_medio,GFlops" > benchmark.csv

# Define arrays for dimensions, threads, and rows
matrix_sizes=(512 1024 2048)
thread_counts=(4 8 16)
block_dimms=(4 16 32)

# Function to run benchmarks
run_benchmarks() {
    local size=$1
    for block_dimm in "${block_dimms[@]}"; do
        for threads in "${thread_counts[@]}"; do
            ./TMulMatCuaBlo-bench $size $block_dimm $threads >> benchmark.csv
        done
    done
}

# Check input parameter and run benchmarks accordingly
if [ "$matrix_dimensions" = "all" ]; then
    for size in "${matrix_sizes[@]}"; do
        run_benchmarks $size
    done
else
    # Ensure the input is one of the expected sizes before running
    if [[ " ${matrix_sizes[*]} " =~ " ${matrix_dimensions} " ]]; then
        run_benchmarks $matrix_dimensions
    else
        echo "Invalid matrix dimension. Please use one of the following: 512, 1024, 2048, all"
        exit 1
    fi
fi
