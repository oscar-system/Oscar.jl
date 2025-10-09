```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Matrix groups

## Introduction

A *matrix group* is a group that consists of invertible square matrices
over a common ring, the *base ring* of the group
(see [`base_ring(G::MatrixGroup)`](@ref)).

We distinguish between *matrices* and *elements of matrix groups*,
that is, matrix group elements are not regarded as matrices.
This distinction is necessary because we want that each element
of a matrix group has this group as its `parent`.
An advantage of this setup is that we know more about a matrix group
element than about a matrix, for example that it is square and invertible.

Apply `matrix` to the matrix group element $g$ in order to create
the corresponding matrix $m$.
Call $G(m)$ in order to create the group element corresponding to $m$.

```jldoctest matgroupxpl
julia> G = general_linear_group(3, 2)
GL(3,2)

julia> x, y = gens(G)
2-element Vector{MatrixGroupElem{FqFieldElem, FqMatrix}}:
 [1 1 0; 0 1 0; 0 0 1]
 [0 0 1; 1 0 0; 0 1 0]

julia> x
[1   1   0]
[0   1   0]
[0   0   1]

julia> m = matrix(x)
[1   1   0]
[0   1   0]
[0   0   1]

julia> m == x
false

julia> G(m) == x
true

julia> parent(x)
GL(3,2)

julia> parent(m)
Matrix space of 3 rows and 3 columns
  over prime field of characteristic 2
```

For convenience, functions for matrices such as `det` and `tr`
are defined also for matrix group elements,
see [the section on functions for elements of matrix groups](@ref "Functions for elements of matrix groups").

```jldoctest matgroupxpl
julia> det(x)
1

julia> tr(x)
1
```

Matrix groups in OSCAR have the type
[`MatrixGroup{RE<:RingElem, T<:MatElem{RE}}`](@ref),
their elements have the type
[`MatrixGroupElem{RE<:RingElem, T<:MatElem{RE}}`](@ref).

## Creation of matrix groups

The typical ways to *create* a matrix group are

- to create matrices that generate the group,
  and then to call `matrix_group` with these generators, or

- to start with [a classical matrix group](@ref "Classical groups"),
  and to create subgroups of it.

```@docs
matrix_group(R::Ring, m::Int, V::AbstractVector{T}; check::Bool=true) where T<:Union{MatElem,MatrixGroupElem}
```

Here are examples for the creation of matrix groups as subgroups
of classical groups.

```jldoctest matgroupxpl
julia> g = GL(3, GF(2))
GL(3,2)

julia> F = GF(4)
Finite field of degree 2 and characteristic 2

julia> G = GL(3, F)
GL(3,4)

julia> is_subset(g, G)
false

julia> base_ring(g)
Prime field of characteristic 2

julia> base_ring(G)
Finite field of degree 2 and characteristic 2

julia> h = change_base_ring(F, g)
Matrix group of degree 3
  over finite field of degree 2 and characteristic 2

julia> is_subset(h, G)
true

julia> flag, mp = is_subgroup(h, G)
(true, Hom: h -> G)

julia> x = gen(h, 1);

julia> parent(x) == h
true

julia> parent(mp(x))
GL(3,4)
```

## Classical groups

```@docs
general_linear_group(n::Int, F::Ring)
special_linear_group(n::Int, F::Ring)
symplectic_group(n::Int, F::Ring)
orthogonal_group(e::Int, n::Int, F::Ring)
special_orthogonal_group(e::Int, n::Int, F::Ring)
omega_group(e::Int, n::Int, F::Ring)
unitary_group(n::Int, q::Int)
special_unitary_group(n::Int, q::Int)
```

## Functions for matrix groups

Defining data of a matrix group are its `base_ring` and its `degree`.

```@docs
base_ring(G::MatrixGroup{RE}) where RE <: RingElem
degree(G::MatrixGroup)
map_entries(f, G::MatrixGroup)
```

## Functions for elements of matrix groups

The following functions delegate to the underlying matrix
of the given matrix group element.

```@docs
matrix(x::MatrixGroupElem)
base_ring(x::MatrixGroupElem)
nrows(x::MatrixGroupElem)
det(x::MatrixGroupElem)
tr(x::MatrixGroupElem)
multiplicative_jordan_decomposition(x::MatrixGroupElem{T}) where T <: FinFieldElem
is_semisimple(x::MatrixGroupElem{T}) where T <: FinFieldElem
is_unipotent(x::MatrixGroupElem{T}) where T <: FinFieldElem
```

## Sesquilinear forms

Sesquilinear forms are alternating and symmetric bilinear forms,
hermitian forms, and quadratic forms.

Sesquilinear forms in OSCAR have the type
[`SesquilinearForm{T<:RingElem}`](@ref).

A sesquilinear form can be created
[from its Gram matrix](@ref "Create sesquilinear forms from Gram matrices")
or [as an invariant form of a matrix group](@ref "Invariant forms").

## Create sesquilinear forms from Gram matrices

```@docs
alternating_form(B::MatElem{T}) where T <: FieldElem
symmetric_form(B::MatElem{T}) where T <: FieldElem
hermitian_form(B::MatElem{T}) where T <: FieldElem
quadratic_form(B::MatElem{T}) where T <: FieldElem
quadratic_form(f::MPolyRingElem{T}) where T <: FieldElem
```

## Invariant forms

The following functions compute (Gram matrices of) sesquilinear forms
that are invariant under the given matrix group.

```@docs
invariant_bilinear_forms(G::MatrixGroup{S,T}) where {S,T}
invariant_sesquilinear_forms(G::MatrixGroup{S,T}) where {S,T}
invariant_quadratic_forms(G::MatrixGroup{S,T}) where {S,T}
invariant_symmetric_forms(G::MatrixGroup{S,T}) where {S,T}
invariant_alternating_forms(G::MatrixGroup{S,T}) where {S,T}
invariant_hermitian_forms(G::MatrixGroup{S,T}) where {S,T}
invariant_bilinear_form(G::MatrixGroup)
invariant_sesquilinear_form(G::MatrixGroup)
invariant_quadratic_form(G::MatrixGroup)
preserved_quadratic_forms(G::MatrixGroup{S,T}) where {S,T}
preserved_sesquilinear_forms(G::MatrixGroup{S,T}) where {S,T}
orthogonal_sign(G::MatrixGroup)
```

## Functions for sesquilinear forms

```@docs
is_alternating(f::SesquilinearForm)
is_hermitian(f::SesquilinearForm)
is_quadratic(f::SesquilinearForm)
is_symmetric(f::SesquilinearForm)
corresponding_bilinear_form(B::SesquilinearForm)
corresponding_quadratic_form(B::SesquilinearForm)
gram_matrix(f::SesquilinearForm)
defining_polynomial(f::SesquilinearForm)
radical(f::SesquilinearForm{T}) where T
witt_index(f::SesquilinearForm{T}) where T
is_degenerate(f::SesquilinearForm{T}) where T
is_singular(f::SesquilinearForm{T}) where T
is_congruent(f::SesquilinearForm{T}, g::SesquilinearForm{T}) where T <: RingElem
isometry_group(f::SesquilinearForm{T}) where T
```

## Utilities for matrices

The inputs/outputs of the functions described in this section
are matrices not matrix group elements.

```@docs
pol_elementary_divisors(A::MatElem{T}) where T
generalized_jordan_block(f::T, n::Int) where T<:PolyRingElem
generalized_jordan_form(A::MatElem{T}; with_pol=false) where T
matrix(A::Vector{AbstractAlgebra.Generic.FreeModuleElem{T}}) where T <: RingElem
conjugate_transpose(x::MatElem{T}) where T <: FinFieldElem
complement(V::AbstractAlgebra.Generic.FreeModule{T}, W::AbstractAlgebra.Generic.Submodule{T}) where T <: FieldElem
permutation_matrix(R::NCRing, Q::AbstractVector{<:IntegerUnion})
is_alternating(M::MatrixElem)
is_hermitian(B::MatElem{T}) where T <: FinFieldElem
```

## Technicalities

```@docs
MatrixGroup{RE<:RingElem, T<:MatElem{RE}}
MatrixGroupElem{RE<:RingElem, T<:MatElem{RE}}
ring_elem_type(::Type{MatrixGroup{S,T}}) where {S,T}
mat_elem_type(::Type{MatrixGroup{S,T}}) where {S,T}
SesquilinearForm{T<:RingElem}
```
