# Circuitscape


| **Documentation** | **Chat** | **Build Status**| **Changelog**|
|:-----------------------------------------------------:|:------------------------------------:|:-----------:|:-------:|
| [![docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://circuitscape.org/docs) [![docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://docs.circuitscape.org/Circuitscape.jl/latest/) | [![gitter](https://badges.gitter.im/Circuitscape/Circuitscape.jl.png)](https://gitter.im/Circuitscape/Circuitscape.jl) | [![Build Status](https://github.com/Circuitscape/Circuitscape.jl/workflows/CI/badge.svg)](https://github.com/Circuitscape/Circuitscape.jl/actions?query=workflow%3ACI) [![Build status](https://ci.appveyor.com/api/projects/status/4a8u8985hq2mt569?svg=true)](https://ci.appveyor.com/project/ranjanan/circuitscape-jl) [![codecov.io](http://codecov.io/github/Circuitscape/Circuitscape.jl/coverage.svg?branch=master)](http://codecov.io/github/Circuitscape/Circuitscape.jl?branch=master) | [![news](https://img.shields.io/static/v1?label=version&message=v5.7.1&color=orange)](https://github.com/Circuitscape/Circuitscape.jl/releases) |

Circuitscape is an open-source program that uses circuit theory to model connectivity 
in heterogeneous landscapes. Its most common applications include modeling movement and gene flow 
of plants and animals, as well as identifying areas important for connectivity conservation. 

Circuitscape has now been rewritten in [Julia](https://julialang.org) for better performance and scalability. Julia is a modern open-source language for scientific computing. 

This work is based on the original [Circuitscape](https://github.com/Circuitscape/Circuitscape) project by Brad McRae, Viral B. Shah 
and Tanmay Mohapatra. 

[Check out the documentation.](https://circuitscape.org/docs)


## The New Circuitscape - Modern, Fast and Flexible

The new Circuitscape is built entirely in the Julia language, a new
programming language for technical computing. Julia is built from the
ground up to be [fast](http://julialang.org/benchmarks). As such, this offers a
number of advantages over the previous version, and these are detailed below.

### Faster and More Scalable

We benchmarked `Circuitscape.jl` (v0.1.0) with the Python version (v4.0.5) to obtain the
following results. We started up Circuitscape with 16 parallel processes,
and used benchmark problems from the standard Circuitscape 
[benchmark suite.](https://github.com/Circuitscape/BigTests)

<img src="docs/src/benchmark/benchmark.png" width=650 height=450>

These benchmarks were run on a Linux (Ubuntu) server machine with the following specs: 
* Name: Intel(R) Xeon(R) Silver 4114 CPU 
* Clock Speed: 2.20GHz
* Number of cores: 20  
* RAM: 384 GB

From the benchmark, we see that the new version is upto *4x faster*
on 16 processes. However, the best performing bar in the chart is 
_Julia-CHOLMOD_, which is a new feature introduced.

### New Solver Mode - CHOLMOD

Julia-CHOLMOD is a new solver mode used in the new Circuitscape. It performs a [cholesky
decomposition](https://en.wikipedia.org/wiki/Cholesky_decomposition) on the graph 
constructed, and performs a batched back substitution
to compute the voltages. It plugs into the 
[CHOLMOD](http://faculty.cse.tamu.edu/davis/suitesparse.html) library, 
which is part of the SuiteSparse collection of high performance sparse 
matrix algorithms.

To use the this new mode, include a line in your Circuitscape 
INI file:
```
solver = cholmod
```

The cholesky decomposition is a direct solver method, unlike the algebraic
multigrid method used by default in both the old and the new version.
The advantage with this new direct method is that it can be much faster than
the iterative solution, within a particular problem size. 

*Word of caution*: The cholesky decomposition is not practical
to use beyond a certain problem size because of phenomenon called
[fill-in](https://algowiki-project.org/en/Cholesky_method#Reordering_to_reduce_the_number_of_fill-in_elements), which results in loss of sparsity and large memory consumption.

### Parallel, everywhere 

The old Circuitscape had limited support for parallelism, which worked on Mac and
Linux, but didn't work on Windows. 

Julia as a programming language is built from the ground up to be parallel,
and as a result the new Circuitscape natively supports parallelism on all three
platforms.

### Single Precision (Experimental)

The new Circuitscape introduces the ability to run problems in
single precision as opposed to the standard double precision.

Single precision usually takes much less memory, but trades off
against solution accuracy. 

Use this new feature by including a line in your config file:
```
precision = single
```

## Installation 

1. You will need to [install](https://julialang.org/downloads/) Julia on your system first.  

2. Once you start Julia, install Circuitscape by: 

```julia
julia> using Pkg
julia> Pkg.add("Circuitscape")
```

If you want the latest development version, you can additionally do: 

```julia
julia> Pkg.add(PackageSpec(name="Circuitscape", rev="master"))
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

## Notes on INI files 

Circuitscape takes as input INI files, which contain paths to the raster map, sources, grounds,
and other inputs, as well as flags for each run. If you're using the [GUI](https://circuitscape.org/downloads/)
the INI file will automatically be generated for you and then fed into Circuitscape. But if you're 
using the Julia prompt, then you must write one yourself. The easiest way to do this is to copy 
an INI file [from the tests](https://github.com/Circuitscape/Circuitscape.jl/tree/master/test/input) and then modify it depending on your problem. 

## Citation

A preprint is available here: https://proceedings.juliacon.org/papers/10.21105/jcon.00058. You can also use the following BibTeX entry to cite this package: 
```bibtex
@article{Anantharaman2020,
  doi = {10.21105/jcon.00058},
  url = {https://doi.org/10.21105/jcon.00058},
  year = {2020},
  publisher = {The Open Journal},
  volume = {1},
  number = {1},
  pages = {58},
  author = {Ranjan Anantharaman and Kimberly Hall and Viral B. Shah and Alan Edelman},
  title = {Circuitscape in Julia: High Performance Connectivity Modelling to Support Conservation Decisions},
  journal = {Proceedings of the JuliaCon Conferences}
}

```
