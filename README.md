# Circuitscape

Linux: [![Build Status](https://travis-ci.org/ranjanan/Circuitscape.jl.svg?branch=master)](https://travis-ci.org/ranjanan/Circuitscape.jl)
[![Coverage Status](https://coveralls.io/repos/ranjanan/Circuitscape.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/ranjanan/Circuitscape.jl?branch=master)
[![codecov.io](http://codecov.io/github/ranjanan/Circuitscape.jl/coverage.svg?branch=master)](http://codecov.io/github/ranjanan/Circuitscape.jl?branch=master)

Circuitscape is an open-source program that uses circuit theory to model connectivity 
in heterogeneous landscapes. Its most common applications include modeling movement and gene flow 
of plants and animals, as well as identifying areas important for connectivity conservation. 

Circuitscape has now been rewritten in [Julia](https://julialang.org) for better performance and scalability. Julia is modern open-source language for scientific computing. 

This work is based on the original [Circuitscape](https://github.com/Circuitscape/Circuitscape) project by Brad McRae, Viral B. Shah 
and Tanmay Mohapatra. 

## The New Circuitscape - Modern, Fast and Flexible 

The new Circuitscape is built entirely in the Julia language, a new
programming language for technical computing. Julia is built from the
ground up to be [fast, scalable and efficient](http://julialang.org/benchmarks).

We benchmarked `Circuitscape.jl` with the Python version to obtain the
following results. 

![benchmark](http://raw.githubusercontent.com/Circuitscape/Circuitscape/readme/benchmark/benchmark.png)

These benchmarks were run on a Linux (Ubuntu) server machine with the following specs: 
* Name: Intel(R) Xeon(R) Silver 4114 CPU 
* Clock Speed: 2.20GHz
* Number of cores: 20  
* RAM: 384 GB

The best performing of all modes is the new Julia-CHOLMOD solver mode. 

### New Solver Mode - CHOLMOD

Julia-CHOLMOD is a new solver mode used in the new Circuitscape. It performs a [cholesky
decomposition]() on the graph constructed, and performs a batched back substitution
to solve the system. It plugs into the [CHOLMOD]() library, which is part of the
SuiteSparse collection of high performance sparse direct solver routines. 

The cholesky decomposition is a direct solver method, unlike 

### Parallel, everywhere 

The old Circuitscape had limited support for parallelism, which worked on Mac and
Linux, but didn't work on Windows. 

Julia as a programming language is built from the ground up to be parallel,
and as a result the new Circuitscape supports parallelism on all three
platforms: Windows, Mac and Linux. 

## Installation 

You will need to [install](https://julialang.org/downloads/) Julia on your system first. 

## Installation

1. First install Julia. 

2. Once you start Julia, install Circuitscape by: 

```julia
julia> Pkg.add("Circuitscape")
```

If you want the latest development version, you can additionally do: 

```julia
julia> Pkg.checkout("Circuitscape")
```

Check if all the tests are passing by doing the following:

```julia
julia> Pkg.test("Circuitscape")
```

## Usage

The current interface to Circuitscape is through the Julia terminal. 

```julia
julia> using Circuitscape # loads the package into your environment
julia> compute("path/to/config/file.ini")
```

## Contributing

If you have encounter any issues or would like to ask a question, please file 
a report [here](https://github.com/ranjanan/Circuitscape.jl/issues).
Contributions in the form of 
[pull requests](https://github.com/ranjanan/Circuitscape.jl/pulls) are also welcome! 
