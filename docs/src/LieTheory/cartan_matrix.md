# Cartan matrices

Cartan matrices can be constructed from a Cartan type, and are represented as a square `ZZMatrix`.

Many functions taking a Cartan matrix as input (like [`root_system`](@ref) and [`weyl_group`](@ref)) will also accept a Cartan type as input. Both cases are semantically equivalent, but the latter may be more efficient.

!!! note
    The convention for Cartan matrices in OSCAR is $(a_{ij}) = (\langle \alpha_i^\vee, \alpha_j \rangle)$ for simple roots $\alpha_i$.

!!! note
    See [Cartan types](@ref) for our conventions on Cartan types and ordering of simple roots.

## Table of contents

```@contents
Pages = ["cartan_matrix.md"]
Depth = 2:5
```

## Constructors

```@docs
cartan_matrix(::Symbol, ::Int)
cartan_matrix(::Vector{Tuple{Symbol,Int}})
```


## Properties

```@docs
is_cartan_matrix(::ZZMatrix)
cartan_symmetrizer(::ZZMatrix)
cartan_bilinear_form(::ZZMatrix)
```


## Cartan types

Cartan types and, in particular, the ordering of simple roots are used differently in the literature.
OSCAR follows the conventions of [Bou02; Plates I-IX](@cite) and [Hum72; Thm 11.4](@cite), i.e. we consider the following Dynkin diagrams
for the Cartan types (where the ordering of simple roots is given by the numbering of the vertices):

- ``A_n`` (for ``n \geq 1``):
  ```
  1 - 2 - 3 - ... - n-1 - n
  ```

- ``B_n`` (for ``n \geq 2``):
  ```
  1 - 2 - 3 - ... - n-1 >=> n
  ```

- ``C_n`` (for ``n \geq 2``):
  ```
  1 - 2 - 3 - ... - n-1 <=< n
  ```

- ``D_n`` (for ``n \geq 4``):
  ```
                        n-1
                      /
  1 - 2 - 3 - ... - n-2
                      \
                        n
  ```

- ``E_6``:
  ```
  1 - 3 - 4 - 5 - 6
          |
          2
  ```

- ``E_7``:
  ```
  1 - 3 - 4 - 5 - 6 - 7
          |
          2
  ```

- ``E_8``:
  ```
  1 - 3 - 4 - 5 - 6 - 7 - 8
          |
          2
  ```

- ``F_4``:
  ```
  1 - 2 >=> 3 - 4
  ```

- ``G_2``:
  ```
  1 <<< 2
  ```

The following function is used to verify the validity of a Cartan type, and is thus used in input sanitization.
```@docs
is_cartan_type(::Symbol, ::Int)
```

Given a Cartan matrix, the following functions can be used to determine its Cartan type.
```@docs
cartan_type(::ZZMatrix)
cartan_type_with_ordering(::ZZMatrix)
```
