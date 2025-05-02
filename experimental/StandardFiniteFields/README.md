# StandardFiniteFieldsOSCAR.jl

This is code written for the OSCAR computer algebra system which aims to
construct finite fields in a standardized way.
The fields are constructed in a way that avoids the usage of Conway polynomials.

The ideas used are described by Lübeck in [[1]](#1)
and some of the code is based on Lübeck's StandardFF GAP package.


## Usage

This code runs in Julia using the OSCAR computer algebra package,
which requires Julia 1.6 or newer. It can be installed as follows.

```julia
julia> using Pkg
julia> Pkg.add("Oscar")
julia> using Oscar
```

For more detailed information, please consult the [installation
instructions](https://www.oscar-system.org/install/) on the OSCAR website.

## References
<a id="1">[1]</a>
Lübeck, F. Standard generators of finite fields and their cyclic subgroups, Journal of Symbolic Computation (2023)(https://arxiv.org/pdf/2107.02257.pdf).
