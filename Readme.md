# OpenMP Mulmat

## Description

This project demonstrates matrix multiplication using OpenMP for parallel processing. It includes implementations for both block and non-block approaches to optimize performance based on different matrix sizes and hardware configurations.

It also includes a Documentation file explaining all implementations and comparing their performances in the same machine. However it is written in spanish as are the comments in code.

## Features

- Parallel matrix multiplication using OpenMP
- Block-based multiplication for optimized performance
- Support for dynamic and static task assignment
- Configurable matrix sizes and block dimensions
- Performance benchmarking with various configurations

## Installation

To compile the project, you will need a C compiler with OpenMP support (such as `gcc`). Clone the repository and use the provided `compilar` file to build the project:

```sh
git clone https://github.com/AlejandroCarmonaM/OpenMP_mulmat.git
cd OpenMP-Matmul
./compilar <program-name>
```

## Usage

You can run the matrix multiplication programs with the following commands:

- For the sequential block matrix multiplication:

```sh
./SecMulMatBlo <matrix_size> <block_size> <num_threads>
```

- For the parallel block matrix multiplication:
```sh
./MulMatCuaBlo <matrix_size> <block_size> <num_threads>
```

- For the parallel task-based block matrix multiplication:
```sh
./TMulMatCuaBlo <matrix_size> <block_size> <num_threads>
```

- Example:
```sh
./TMulMatCuaBlo 1024 16 8
```