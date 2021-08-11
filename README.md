# JToric

[![Build Status](https://github.com/HereAround/JToric.jl/workflows/CI/badge.svg)](https://github.com/HereAround/JToric.jl/actions)
[![Coverage](https://codecov.io/gh/HereAround/JToric.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/HereAround/JToric.jl)

## What is JToric about?

This package aims to make the functionality of the `ToricVarieties_project` (see https://github.com/homalg-project/ToricVarieties_project) available in `Julia`.


## Crucial dependency

This crucially hings on a functional Gap4-package `NConvex`. Unfortuantely, the current version of this package cannot be loaded in `Julia`, as it depends on `CddInterface` and `NormalizInterface`, which each are currently not ladable in `Julia`.


## Installation

A quick way to satisfy the dependency of `JToric` on a functional `NConvex`, is to replace `NConvex` with a development version, which does not depend on `CddInterface` and `NormalizInterface`. Such a version is prepared in the `martindevel` branch of the repository https://github.com/HereAround/NConvex.git. To use this version, please proceed as follows:

1. Start `Julia`
2. Issue the tollowing two lines of code:

    `using GAP`
    `gap_location = GAP.GAPROOT`

2. Remember this path gap_location.
3. Exit julia.
4. `cd $(gap_location)/pkg`
5. `git clone https://github.com/HereAround/NConvex.git`
6. `git fetch origin martindevel:martindevel`

Restart `Julia` and, if necessary, install the packages `Singular.jl`, `Polymake.jl`, `Oscar.jl` and `CapAndHomalg.jl`:

7. Start julia
8. `using Pkg`
9. `Pkg.add( "Singular" )`
10. `Pkg.add( Polymake )`
11. `Pkg.add( Oscar )`
12. `Pkg.add( CapAndHomalg )`

Finally, remember that your development version of `JToric` should be placed outside of the `.julia` folder of your home folder. If not done yet, please move it to such a location now and register it in `julia` via:

13. `Pkg.develop( path = "path/to/your/JToric" )`

This completes the installation.


## Towards a functional CddInterface

Eventually, we hope to have a functional `CddInterface` and `NormalizInterface`, so that no changes to `NConvex` are necessary. Currently, a working installation of `CddInterface` can be achieved as follows:

1. Start Julia
2. Issue the following:

    `using GAP`
    `gap_location = GAP.GAPROOT`

3. Remember this path `gap_location`.
4. Exit julia.
5. `cd $(gap_location)/pkg`
6. `git clone https://github.com/homalg-project/CddInterface.git`
7. `git reset --hard b806a6f93788c38b0323f642bdb220287b9fc41e`
8. `cd CddInterface`
9. Open the file `PackageInfo.g` and edit the version number by hand to "2020.06.24".
10. Install CddInterface by issuing
      `./install.sh $(gap_location)`

This should provide a functional `CddInterface`.

## Towards a functional NormalizInterface

 Even after having a functional `CddInterface`, you cannot load `NConvex` unless you also establish a functional `NormalizInterface`. It an be anticipated that steps similar to the ones performed for `CddInterface` will achieve this. However, this has not yet been tested.
