```@contents
Pages = ["architecture.md"]
```

```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

# Architecture

This page aims to give a short technical overview of the architecture of Oscar.
A more in-depth overview of the various components of Oscar is given
[on the Oscar homepage](https://oscar.computeralgebra.de/about/).

## Julia packages

Oscar is developed as a pure julia package Oscar.jl and builds on the features
and interfaces provided by the julia packages for the cornerstones:
  - ANTIC: [AbstractAlgebra.jl](https://github.com/Nemocas/AbstractAlgebra.jl),
           [Hecke.jl](https://github.com/thofma/Hecke.jl),
           [Nemo.jl](https://github.com/nemocas/Nemo.jl)
  - [GAP.jl](https://github.com/oscar-system/GAP.jl)
  - [Polymake.jl](https://github.com/oscar-system/Polymake.jl)
  - [Singular.jl](https://github.com/oscar-system/Singular.jl)

The packages are integrated into the julia package manager and will be
installed automatically as dependencies of Oscar. They can be accessed directly
by their names once Oscar is loaded.

The current versions of these packages can be inspected with the `oscar_versioninfo` command:

```@repl oscar
oscar_versioninfo()
```

## Binary packages for non-julia libraries

In addition to the pure julia packages, Oscar builds on many non-julia libraries
for the cornerstones ([FLINT](http://flintlib.org/), [GAP](https://gap-system.org/),
[polymake](https://polymake.org), [Singular](https://www.singular.uni-kl.de/))
and their dependencies.

All dependencies have been integrated into the [BinaryBuilder](https://binarybuilder.org/)
ecosystem which provides precompiled binaries for all supported architectures.
The build-recipes are maintained in julia's [Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil)
repository. These are used to automatically build binary artifacts together with the
corresponding `jll` packages which allow automatic installation of all required
binaries as dependencies for Oscar.

Both Polymake.jl and Singular.jl use [CxxWrap.jl](https://github.com/JuliaInterop/CxxWrap.jl)
together with the corresponding `libcxxwrap-julia` library as an intermediate layer between the
julia packages and the C / C++ libraries.


The `oscar_versioninfo` function can also include the versions of all binary packages that are
maintained by the Oscar developers:

```@repl oscar
oscar_versioninfo(jll=true)
```

For a full list of all dependencies of the current project please use
`using Pkg; Pkg.status(mode=PKGMODE_MANIFEST)` or the Pkg REPL command `]st -m`.


```@docs
oscar_versioninfo
```
