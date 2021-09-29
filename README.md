# JToric

[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://oscar-system.github.io/JToric.jl/dev)
[![CI](https://github.com/oscar-system/JToric.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/HereAround/JToric/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/oscar-system/JToric.jl/branch/master/graph/badge.svg?token=z9gKk5v3O6)](https://codecov.io/gh/oscar-system/JToric.jl)


## What is JToric about?

This package aims to make functionality for computing with toric varieties available in `Julia` and the computer algebra system `Oscar` (https://oscar.computeralgebra.de/). In particular, the functionality provided by the `ToricVarieties_project` (<https://github.com/homalg-project/ToricVarieties_project>) as well as `polymake` (https://polymake.org/doku.php) will be made available.

## Documentation

The documentation is available at https://oscar-system.github.io/JToric.jl/dev/.

## Installation

First, install `Julia` on your computer. Next, this development version of `JToric` should be placed outside of the `.julia` folder of your home folder. If not done yet, please move it to such a location now and register and build it as follows:

    using Pkg
    Pkg.develop( path = "path/to/your/JToric" )
    Pkg.build( "JToric" )
   
