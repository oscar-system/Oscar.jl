```@meta
CurrentModule = Oscar
CollapsedDocStrings = true
DocTestSetup = Oscar.doctestsetup()
```

# Matrix groups

## Introduction

A *matrix group* is a group that consists of invertible square matrices
over a common ring, the *base ring* of the group
(see [`base_ring(G::MatrixGroup{RE}) where RE <: RingElem`](@ref)).

Matrix groups in OSCAR have the type
[`MatrixGroup{RE<:RingElem, T<:MatElem{RE}}`](@ref),
their elements have the type
[`MatrixGroupElem{RE<:RingElem, T<:MatElem{RE}}`](@ref).

## Basic Creation

In order to *create* a matrix group,
one first creates matrices that generate the group,
and then calls `matrix_group` with these generators.

```jldoctest matgroupxpl
julia> mats = [[0 -1; 1 -1], [0 1; 1 0]]
2-element Vector{Matrix{Int64}}:
 [0 -1; 1 -1]
 [0 1; 1 0]

julia> matelms = map(m -> matrix(ZZ, m), mats)
2-element Vector{ZZMatrix}:
 [0 -1; 1 -1]
 [0 1; 1 0]

julia> g = matrix_group(matelms)
Matrix group of degree 2
  over integer ring

julia> describe(g)
"S3"

julia> t = matrix_group(ZZ, 2, typeof(matelms[1])[])
Matrix group of degree 2
  over integer ring

julia> t == trivial_subgroup(g)[1]
true

julia> F = GF(3); matelms = map(m -> matrix(F, m), mats)
2-element Vector{FqMatrix}:
 [0 2; 1 2]
 [0 1; 1 0]

julia> g = matrix_group(matelms)
Matrix group of degree 2
  over prime field of characteristic 3

julia> describe(g)
"S3"
```

Alternatively,
one can start with [a classical matrix group](@ref "Classical groups").

```jldoctest matgroupxpl
julia> g = GL(3, GF(2))
GL(3,2)

julia> F = GF(4)
Finite field of degree 2 and characteristic 2

julia> G = GL(3, F)
GL(3,4)

julia> is_subset(g, G)
false

julia> h = change_base_ring(F, g)
Matrix group of degree 3
  over finite field of degree 2 and characteristic 2

julia> flag, mp = is_subgroup(h, G)
(true, Hom: h -> G)

julia> mp(gen(h, 1)) in G
true
```

## Functions for matrix groups

```@docs
matrix_group(R::Ring, m::Int, V::AbstractVector{T}; check::Bool=true) where T<:Union{MatElem,MatrixGroupElem}
base_ring(G::MatrixGroup{RE}) where RE <: RingElem
degree(G::MatrixGroup)
centralizer(G::MatrixGroup{T}, x::MatrixGroupElem{T}) where T <: FinFieldElem
map_entries(f, G::MatrixGroup)
```

## Elements of matrix groups

The following functions delegate to the underlying matrix
of the given matrix group element.

```@docs
matrix(x::MatrixGroupElem)
base_ring(x::MatrixGroupElem)
nrows(x::MatrixGroupElem)
det(x::MatrixGroupElem)
tr(x::MatrixGroupElem)
multiplicative_jordan_decomposition(x::MatrixGroupElem)
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

## Utilities for matrices

(The inputs/outputs of the functions described in this section
are matrices not matrix group elements.)

```@docs
pol_elementary_divisors(A::MatElem{T}) where T
generalized_jordan_block(f::T, n::Int) where T<:PolyRingElem
generalized_jordan_form(A::MatElem{T}; with_pol=false) where T
matrix(A::Vector{AbstractAlgebra.Generic.FreeModuleElem{T}}) where T <: RingElem
upper_triangular_matrix(L::AbstractVector{T}) where {T <: NCRingElement}
lower_triangular_matrix(L::AbstractVector{T}) where {T <: NCRingElement}
conjugate_transpose(x::MatElem{T}) where T <: FinFieldElem
complement(V::AbstractAlgebra.Generic.FreeModule{T}, W::AbstractAlgebra.Generic.Submodule{T}) where T <: FieldElem
permutation_matrix(F::Ring, Q::AbstractVector{<:IntegerUnion})
is_alternating(B::MatElem)
is_hermitian(B::MatElem{T}) where T <: FinFieldElem
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

## Technicalities

```@docs
MatrixGroup{RE<:RingElem, T<:MatElem{RE}}
MatrixGroupElem{RE<:RingElem, T<:MatElem{RE}}
ring_elem_type(::Type{MatrixGroup{S,T}}) where {S,T}
mat_elem_type(::Type{MatrixGroup{S,T}}) where {S,T}
SesquilinearForm{T<:RingElem}
```
