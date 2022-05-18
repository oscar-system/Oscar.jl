```@meta
CurrentModule = Oscar
DocTestSetup = quote
  using Oscar
end
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["matgroup.md"]
```

# Matrix groups

```@docs
MatrixGroup{RE<:RingElem, T<:MatElem{RE}}
MatrixGroupElem{RE<:RingElem, T<:MatElem{RE}}
base_ring(G::MatrixGroup)
degree(G::MatrixGroup)
centralizer(G::MatrixGroup{T}, x::MatrixGroupElem{T}) where T <: FinFieldElem
```

## Elements of matrix groups

```@docs
matrix(x::MatrixGroupElem)
base_ring(x::MatrixGroupElem)
nrows(x::MatrixGroupElem)
det(x::MatrixGroupElem)
trace(x::MatrixGroupElem)
multiplicative_jordan_decomposition(x::MatrixGroupElem)
is_semisimple(x::MatrixGroupElem{T}) where T <: FinFieldElem
is_unipotent(x::MatrixGroupElem{T}) where T <: FinFieldElem
```

## Sesquilinear forms

```@docs
SesquilinearForm{T<:RingElem}
is_alternating_form(f::SesquilinearForm)
is_hermitian_form(f::SesquilinearForm)
is_quadratic_form(f::SesquilinearForm)
is_symmetric_form(f::SesquilinearForm)
alternating_form(B::MatElem{T}) where T <: FieldElem
symmetric_form(B::MatElem{T}) where T <: FieldElem
hermitian_form(B::MatElem{T}) where T <: FieldElem
quadratic_form(B::MatElem{T}) where T <: FieldElem
quadratic_form(f::MPolyElem{T}) where T <: FieldElem
corresponding_bilinear_form(B::SesquilinearForm)
corresponding_quadratic_form(B::SesquilinearForm)
gram_matrix(f::SesquilinearForm)
defining_polynomial(f::SesquilinearForm)
radical(f::SesquilinearForm{T}) where T
witt_index(f::SesquilinearForm{T}) where T
is_degenerate(f::SesquilinearForm{T}) where T
is_singular(f::SesquilinearForm{T}) where T
is_congruent(f::SesquilinearForm{T}, g::SesquilinearForm{T}) where T <: RingElem
```

## Invariant forms

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
isometry_group(f::SesquilinearForm{T}) where T
orthogonal_sign(G::MatrixGroup)
```

## Utilities for matrices (replace by available functions, or document elsewhere?)

```@docs
pol_elementary_divisors(A::MatElem{T}) where T
generalized_jordan_block(f::T, n::Int) where T<:PolyElem
generalized_jordan_form(A::MatElem{T}; with_pol=false) where T
matrix(A::Vector{AbstractAlgebra.Generic.FreeModuleElem{T}}) where T<: FieldElem
upper_triangular_matrix(L)
lower_triangular_matrix(L)
conjugate_transpose(x::MatElem{T}) where T <: FinFieldElem
complement(V::AbstractAlgebra.Generic.FreeModule{T}, W::AbstractAlgebra.Generic.Submodule{T}) where T <: FieldElem
permutation_matrix(F::Ring, Q::AbstractVector{<:IntegerUnion})
is_skewsymmetric_matrix(B::MatElem{T}) where T <: RingElem
is_hermitian_matrix(B::MatElem{T}) where T <: FinFieldElem
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
matrix_group(V::AbstractVector{T}) where T<:Union{MatElem,MatrixGroupElem}
```
