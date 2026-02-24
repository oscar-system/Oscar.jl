```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Introduction

Implementation of the algorithms from [Hof26](@cite) for computing unit groups and the first K-group of finite rings and finite associative algebras.

## Status

This part of OSCAR is in an experimental state; please see [Adding new projects to experimental](@ref) for what this means.

## Contact

You can ask questions in the [OSCAR Slack](https://www.oscar-system.org/community/#slack).

Alternatively, you can [raise an issue on github](https://www.oscar-system.org/community/#how-to-report-issues).

## Examples

```jldoctest
julia> R, = finite_ring(GF(2)[symmetric_group(4)]); # this the modular group ring F_2[S_4]

julia> U, mU = unit_group(R)
(Finitely presented group of order 3145728, Map: U -> finite ring)

julia> R, = finite_ring(GF(2)[small_group(4, 2)]);

julia> Oscar.k1(R)
((Z/2)^3, Map: (Z/2)^3 -> finite ring)
```
