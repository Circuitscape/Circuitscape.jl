# Circuitscape

Linux: [![Build Status](https://travis-ci.org/ranjanan/Circuitscape.jl.svg?branch=master)](https://travis-ci.org/ranjanan/CircuitScape.jl)
[![Coverage Status](https://coveralls.io/repos/ranjanan/Circuitscape.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/ranjanan/CircuitScape.jl?branch=master)
[![codecov.io](http://codecov.io/github/ranjanan/Circuitscape.jl/coverage.svg?branch=master)](http://codecov.io/github/ranjanan/CircuitScape.jl?branch=master)

Circuitscape is an open-source program that uses circuit theory to model connectivity 
in heterogeneous landscapes. Its most common applications include modeling movement and gene flow 
of plants and animals, as well as identifying areas important for connectivity conservation. 

Circuitscape has now been rewritten in [Julia](https://julialang.org) for better performance and scalability. Julia is modern open-source language for scientific computing. 

## Requirements

You will need to [install](https://julialang.org/downloads/) Julia on your system first. 

## Installation

1. First install Julia. 

2. Once you start Julia, install Circuitscape by: 

```julia
Pkg.add("Circuitscape")
```

If you want the latest development version, you can additionally do: 

```julia
Pkg.checkout("Circuitscape")
```

## Usage

The current interface to Circuitscape is through the Julia terminal. 

```julia
using Circuitscape # loads the package into your environment
compute("path/to/config/file.ini")
```

## Contributing

If you have encounter any issues or would like to ask a question, please file 
a report [here](https://github.com/ranjanan/Circuitscape.jl/issues).
Contributions in the form of 
[pull requests](https://github.com/ranjanan/Circuitscape.jl/pulls) are also welcome! 
