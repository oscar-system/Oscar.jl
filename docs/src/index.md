# Welcome to FTheoryTools

## Goal

The package FTheoryTools.jl aims to automate a number of recurring and (at least in part) tedious computations in F-theory model building
as described in large detail in [Wei18](@cite). Specifically, we focus on the following setups:
* 4d F-theory compactifications,
* defined by a global and singular Weierstrass model as codimension 1 locus of a toric ambient space $Y$,
* which can be crepantly resolved.
For possible future extensions, see below.

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


## Documentation

For detailed information about the implemented functionality, please take a look at the most recent [documentation](https://julia-meets-string-theory.github.io/FTheoryTools.jl/dev/).


## Tests

We aim for a high test coverage. This should be conducted with good and demanding examples, which at the same time do not exceed certain
computational times. To this end, we use a particular toric 3-fold, which we have identified to serve this purpose.
```@docs
TestBase()
```


## Bugs and feature requests

If you want to report a bug or request a feature, please do it by raising a [github issue](https://github.com/Julia-meets-String-Theory/FTheoryTools.jl/issues).
```@docs
FTheoryTools.version
```


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
  author       = "Bies, Martin and Andrew, P. Turner",
  howpublished = {\url{https://github.com/Julia-meets-String-Theory/FTheoryTools.jl}},
  title        = {FTheoryTools -- Julia tools for F-Theory compactifications},
  year         = "2022"
}
```


## Funding

The work of Martin Bies is supported by SFB-TRR 195 `Symbolic Tools in Mathematics and their Application of the German Research Foundation (DFG)`.

The work of Andrew Turner is supported by DOE (HEP) Award DE-SC001352.
