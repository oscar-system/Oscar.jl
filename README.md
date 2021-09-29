# JToric

[![CI](https://github.com/oscar-system/JToric.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/HereAround/JToric/actions/workflows/CI.yml)
[![codecov](https://codecov.io/gh/oscar-system/JToric.jl/branch/master/graph/badge.svg?token=z9gKk5v3O6)](https://codecov.io/gh/oscar-system/JToric.jl)


## What is JToric about?

This package aims to make the functionality of the `ToricVarieties_project` (see <https://github.com/homalg-project/ToricVarieties_project>) available in `Julia`.

## Documentation

The documentation is available at https://oscar-system.github.io/JToric.jl/dev/.

## Crucial dependency

This crucially hinges on a functional Gap4-package `NConvex`. Unfortunately, the current version of this package cannot be loaded in `Julia`, as it depends on `CddInterface` and `NormalizInterface`, which each are currently not loadable in `Julia`.


## Installation

First, install `Julia` on your computer. Next, this development version of `JToric` should be placed outside of the `.julia` folder of your home folder. If not done yet, please move it to such a location now and register it in `julia` and build it, via:

    using Pkg
    Pkg.develop( path = "path/to/your/JToric" )
    Pkg.build( "JToric" )


## Details of the build step

As mentioned above, `JToric` crucially depends on `NConvex`.  A quick way to satisfy this dependency, is to replace `NConvex` with a development version, which is functional in `Julia`, i.e. does not depend on `CddInterface` and `NormalizInterface`. Such a version is prepared in the `martindevel` branch of the repository https://github.com/HereAround/NConvex.git. The build step tries to ensure that the `gap` used in `Julia` employs this very version of `NConvex`.

Explicitly, you can reproduce the build process by following the following steps:

1. Start `Julia`
2. Issue the following code:

        using GAP; gap_location = GAP.GAPROOT

2. Remember this path `gap_location`.
3. Exit julia.
4. `cd $(gap_location)/pkg`
5. `git clone -b martindevel https://github.com/HereAround/NConvex.git`

In `$(gap_location)/pkg` you should now have the development version of `NConvex` which does not depend on `CddInterface` and `NormalizInterface`. Since its version is higher than the distributed `NConvex`, this version should be used by `gap` in `julia`.


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

 Even after having a functional `CddInterface`, you cannot load the distributed version of `NConvex` (but yes, the above development version can) unless you also establish a functional `NormalizInterface`. It can be anticipated that steps similar to the ones performed for `CddInterface` will achieve this. However, this has not yet been tested.
