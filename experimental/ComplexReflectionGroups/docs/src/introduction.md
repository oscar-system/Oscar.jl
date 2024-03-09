```@meta
CurrentModule = Oscar
```

# Introduction

This project implements complex reflection groups [LT09](@cite).

## Features

The list of features at the moment is:

* Explicit matrix group models of complex reflection groups, especially various models for
  the exceptional groups: the unitary models from Lehrer & Taylor [LT09](@cite), the
  invariant models from Marin & Michel [MM10](@cite) as implemented in CHEVIE
  [Mic15](@cite), the integral models as implemented by Taylor in Magma [BCP97](@cite).

* A structure for the type of a complex reflection group, i.e., its
  $\mathrm{GL}_n(\mathbb{C})$-conjugacy class as labeled by Shephard & Todd [ST54](@cite).
  This allows to implement functions and data for complex reflection groups that do not
  depend on an explicit representative like the order of the group, the number of
  reflections, etc. The explicit matrix group models know their type and this avoids heavy
  matrix group computations for data that is already in the literature.

* A structure that encapsulates data of a complex reflection and functions determining such
  data for a reflection given by a matrix: root, coroot, hyperplane, non-trivial eigenvalue,
  order, whether it is unitary.


## Showcase

```jldoctest
julia> W = complex_reflection_group_type([33, (8,4,6)])
Complex reflection group type G33 x G(8,4,6)

julia> order(W)
2446118092800

julia> num_reflections(W)
171

julia> complex_reflection_group(W)
Direct product of
 Matrix group of degree 5 over cyclotomic field of order 3
 Matrix group of degree 6 over cyclotomic field of order 8
```

## Future plans

In the future we also want to implement symplectic reflection groups [Coh80](@cite) and all
sorts of objects associated to reflections groups like the world of rational Cherednik
algebras, symplectic reflection algebras, Calogero--Moser spaces [EG02](@cite). One goal is
to obtain new results about symplectic singularities [Bea00](@cite).

## Contact

Please direct questions about this part of OSCAR to the following people:

* [Ulrich Thiel](https://ulthiel.com/math)

## Acknowledgements

This work is a contribution to the SFB-TRR 195 '[Symbolic Tools in Mathematics and their
Application](https://www.computeralgebra.de/sfb/)' of the German Research Foundation (DFG).
