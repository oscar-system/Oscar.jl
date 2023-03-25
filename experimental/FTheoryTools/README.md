# FtheoryTools.jl

## Goal

The package FTheoryTools.jl aims to automate a number of recuring and (at least in part) tedious computations in F-theory model building. Specifically we focus on the following setups:
* 4d F-theory compactifications,
* defined by a global and singular Weierstrass model as codimension 1 locus of a toric ambient space $Y$,
* which can be crepantly resolved.
Some of the techniques/algorithms extend naturally to more general settings. For example, it is not at all necessary to restrict to 4-dimensional settings (or alternatively, base spaces of dimension 3). Indeed, the current implementation does allow for *arbitrary* base dimension. For more extensions that we might address in the future, please take a look at the section "possibly future extensions" below.

We aim for the following workflow:
* User input:
    * Weierstrass polynomial $P_W$,
    * Data defining the toric ambient space $Y$ (if applicable),
    * Choice of resolved phase (if applicable),
    * Generating sections (for $U(1)$ symmetries).
* Output:
    * Singular loci in codimension 1, 2 and 3.
    * Defining data of resolved geometry.
    * (Pictures of) fibre diagrams of resolved fibre over the originally singular loci, including intersections of $U(1)$-sections.
    * Gauge group.
    * Topological data (e.g. Euler number).

## Status

This project just began. We hope to have a first working version by the end of the year 2022.

## Possible future extensions

Future extensions include, but are not necessarily limited to, the following:
* Specify a $G_4$-flux and work-out the chiral spectra.
* Specify a gauge potential and work out (candidates for) the line bundles whose cohomologies encode the vector-like spectra.
* Other singularity types (non-minimal, terminal, etc.)
* Base blowups for singularity resolution.


## Dependencies

We base this project on [OSCAR](https://oscar.computeralgebra.de/) for general functionality on toric spaces and (possibly even more importantly) polynomial operations. The latter are based on [Singular](https://www.singular.uni-kl.de/) and [Singular.jl](https://github.com/oscar-system/Singular.jl), respectively.


## Code coverage and tests

* [![CI](https://github.com/Julia-meets-String-Theory/FTheoryTools.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/Julia-meets-String-Theory/FTheoryTools.jl/actions/workflows/CI.yml),
* [![pages-build-deployment](https://github.com/Julia-meets-String-Theory/FTheoryTools.jl/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/Julia-meets-String-Theory/FTheoryTools.jl/actions/workflows/pages/pages-build-deployment),
* [![codecov](https://codecov.io/gh/Julia-meets-String-Theory/FTheoryTools.jl/branch/master/graph/badge.svg?token=T5456HQGYZ)](https://codecov.io/gh/Julia-meets-String-Theory/FTheoryTools.jl).


## Installation instructions for Linux

1. Install `Julia` on your computer. The latest version can be found [here](https://julialang.org/downloads/).
2. Download this development version of `FTheoryTools.jl`. Those interested in contributing should instead clone this repository:
```
    git clone https://github.com/Julia-meets-String-Theory/FTheoryTools.jl.git
```
3. Place your clone/download in a location outside of the `.julia` folder of your home folder.
4. Finally, register and build `FTheoryTools.jl` as follows:
```julia
    using Pkg
    Pkg.add("Oscar")
    Pkg.develop(path="path/to/your/FTheoryTools.jl")
    Pkg.build("FTheoryTools")
```
5. To use your `FTheoryTools` installation, you finally want to invoke
```julia
    using Oscar, FTheoryTools
```

Note that due to the modularity of the `Julia` programing language, access to the `Oscar`
functionality is not available in the REPL unless you issue `using("Oscar")` explicitly.
This is the reason for why this command appears in the examples of our documentation.

You can instruct `Julia` to automatically conduct the above steps upon starting, so that you
immediately have access to `FTheoryTools` (and `Oscar`). To this end perform the following steps:
1. Navigate to the *hidden* `julia` folder in your home directory (e.g. via `cd .julia/` in the terminal).
2. In this folder, create (if not already existing) the folder `config`.
3. In the `config` folder create (unless existing) the file `startup.jl`.
4. Add the following lines to the file `startup.jl`:
```julia
using Pkg
using Oscar
Pkg.develop(path="path/to/your/FTheoryTools.jl")
using FTheoryTools
```


## Documentation

For detailed information about the implemented functionality, please take a look at the most recent [documentation](https://julia-meets-string-theory.github.io/FTheoryTools.jl/dev/).


## Bugs and feature requests

If you want to report a bug or request a feature, please do it by raising a [github issue](https://github.com/Julia-meets-String-Theory/FTheoryTools.jl/issues).


## Contributions

Contributions are highly appreciated. Please notice that:
* Contributions must be done by opening a [pull request](https://github.com/Julia-meets-String-Theory/FTheoryTools.jl/pulls).
* Pull requests must pass a number of checks that are automatically conducted by our test suite, before they can be merged. A further approval by a code owner is appreciated.
* Code is expected to be in agreement with the [Oscar style guide](https://oscar-system.github.io/Oscar.jl/stable/DeveloperDocumentation/styleguide/).


## Contact

This software is work in progress of
* [Martin Bies](https://martinbies.github.io/),
* [Andrew Turner](https://apturner.net/).

We appreciate insightful discussions with [Mirjam Cvetič](https://live-sas-physics.pantheon.sas.upenn.edu/people/standing-faculty/mirjam-cvetic).

If you are interested in contributing, please feel free to reach out to us for more details.


## Citation

If `FTheoryTools` has been used in your research, please use the bibtex below or download it [here](https://github.com/Julia-meets-String-Theory/FTheoryTools.jl/blob/master/bibtex.bib):
```bibtex
@Misc{FTheoryTools,
  author       = "Bies, Martin and Turner, Andrew P.",
  howpublished = {\url{https://github.com/Julia-meets-String-Theory/FTheoryTools.jl}},
  title        = {FTheoryTools -- Julia tools for F-Theory compactifications},
  year         = "2022"
}
```


## Funding

The work of Martin Bies is supported by SFB-TRR 195 ``Symbolic Tools in Mathematics and their Application of the German Research Foundation (DFG)``.

The work of Andrew Turner is supported by DOE (HEP) Award DE-SC001352.
