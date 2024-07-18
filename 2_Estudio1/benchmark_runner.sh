#!/bin/bash

# Usage: ./benchmark_runner.sh <matrix_dimensions>
# Possible values for matrix_dimensions: 512, 1024, 2048, all
# Example: ./benchmark_runner.sh 512

# Retrieve input parameter
matrix_dimensions=$1

# Compile MulMatCua-bench
./compilar MulMatCua-bench

# Write CSV header
echo "dimension_matriz_n,numero_filas_consecutivas_F,numero_hilos_t,tiempo_medio,GFlops" > benchmark.csv

# Define arrays for dimensions, threads, and rows
matrix_sizes=(512 1024 2048)
thread_counts=(1 4 8 16)
row_counts=(1 4 16 32)

# Function to run benchmarks
run_benchmarks() {
    local size=$1
    for rows in "${row_counts[@]}"; do
        for threads in "${thread_counts[@]}"; do
            ./MulMatCua-bench $size $rows $threads >> benchmark.csv
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
