```@meta
CurrentModule = Oscar
DocTestSetup = Oscar.doctestsetup()
```

# Introduction

This project implements finite complex reflection groups, following the viewpoint of
[LT09](@cite): a reflection group is a linear representation, not only an abstract finite
group.

Three levels of data are kept separate:

* For an irreducible group in $\mathrm{GL}_n(\mathbb C)$, a
  `ComplexReflectionGroupType` describes its conjugacy class in
  $\mathrm{GL}_n(\mathbb C)$. Shephard--Todd symbols are classification labels for these
  classes. For an essential reducible group, the type stores the normalized unordered
  product of the irreducible types. It stores no matrices, basis, or trivial fixed summand.
  Thus the present object describes essential conjugacy types; a future type for an
  arbitrary ambient realization must additionally record the dimension of its fixed space.
* A **reference marking** adds conventional labels to the components and
  reflecting-hyperplane orbits, and chooses an orientation of each cyclic pointwise
  stabilizer. Its nonidentity powers label the reflection classes. Representative words use
  the generator alphabet of the CHEVIE reference model; these labels are not invariants of
  an unmarked type.
* A **concrete realization** is a matrix group returned by
  `complex_reflection_group`. The selected model determines its matrices,
  coefficient field, generators, and any supplied invariant Hermitian form.

## Features

The list of features at the moment is:

* Explicit matrix group realizations, especially several models for the exceptional groups:
  the standard-unitary models from Lehrer--Taylor [LT09](@cite), the invariant models of
  Marin--Michel [MM10](@cite) as implemented in CHEVIE [Mic15](@cite), and models with
  algebraic-integer matrix entries as implemented by Taylor in Magma [BCP97](@cite).

* Type-level invariants, including order, rank, degrees, and the numbers of reflections and
  reflecting hyperplanes. The labels for the conjugacy types follow the Shephard--Todd
  classification [ST54](@cite). Concrete realizations retain their type, avoiding
  matrix-group computations for data already determined by the classification.

* Type-level and concrete descriptions of reflecting-hyperplane orbits and reflection
  conjugacy classes. These record the hierarchy from a hyperplane orbit, through its cyclic
  pointwise stabilizer, to the classes represented by the nonidentity powers of a marked
  generator. The same reference orientation is transported across the built-in matrix
  models.

* A structure encapsulating the data of a complex reflection and functions recovering that
  data from a matrix: root, coroot, hyperplane, nontrivial eigenvalue, order, and whether the
  matrix preserves the standard Hermitian form.


## Showcase

```jldoctest
julia> T = complex_reflection_group_type([33, (8, 4, 6)])
Complex reflection group type G(8,4,6) x G33

julia> order(T)
2446118092800

julia> rank(T)
11

julia> number_of_reflections(T)
171

julia> complex_reflection_group(T)
Matrix group of degree 11
  over cyclotomic field of order 24
```

## Future plans

The component also contains initial constructions for symplectic reflection groups
[Coh80](@cite). Further development will connect the reflection data to rational Cherednik
algebras, symplectic reflection algebras, and Calogero--Moser spaces [EG02](@cite), with
applications to symplectic singularities [Bea00](@cite).

## Contact

Please direct questions about this part of OSCAR to the following people:

* [Ulrich Thiel](https://ulthiel.com/math)

## Acknowledgements

This work is a contribution to the SFB-TRR 195 '[Symbolic Tools in Mathematics and their
Application](https://www.computeralgebra.de/sfb/)' of the German Research Foundation (DFG).
