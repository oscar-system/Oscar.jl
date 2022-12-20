```@meta
CurrentModule = Oscar
```

```@setup oscar
using Oscar
```

```@contents
Pages = ["matroids.md"]
```

# Matroids

## Introduction

Matroids are a fundamental combinatorial object with connections to various fields of mathematics.
It is an abstraction of linear independence in vector spaces and forests in graphs.
One way to define a *matroid* is via the following two sets of data:

- a finite *ground set* $E := \{1,\ldots,n\}$ and
- a nonempty finite set $\mathcal{B} \subseteq \mathcal{P}(E)$ of *bases* satisfying an exchange property.

There are however many equivalent ways to define a matroid.
One can also define a matroid via its *circuits*, *hyperplanes*, a *graph*, or a *matrix*.
For a detailed introduction of matroids we refer to the textbook [Oxl11](@cite).

## Construction

```@docs
matroid_from_bases(bases::AbstractVector{T}, nelements::IntegerUnion; check::Bool=true) where T<:GroundsetType
matroid_from_nonbases(nonbases::AbstractVector{T}, nelements::IntegerUnion; check::Bool=true) where T<:GroundsetType
matroid_from_circuits(circuits::AbstractVector{T}, nelements::IntegerUnion) where T<:GroundsetType
matroid_from_hyperplanes(hyperplanes::AbstractVector{T}, nelements::IntegerUnion) where T<:GroundsetType
matroid_from_matrix_columns(A::MatrixElem; check::Bool=true)
matroid_from_matrix_rows(A::MatrixElem, ; check::Bool=true)
cycle_matroid(g::Graph)
bond_matroid(g::Graph)
cocycle_matroid(g::Graph)
Matroid(pm_matroid::Polymake.BigObjectAllocated, E::GroundsetType=Vector{Integer}(1:pm_matroid.N_ELEMENTS))
matroid_from_revlex_basis_encoding(rvlx::String, r::IntegerUnion, n::IntegerUnion)
```

## Examples

```@docs
uniform_matroid(r::IntegerUnion,n::IntegerUnion)
fano_matroid()
non_fano_matroid()
non_pappus_matroid()
pappus_matroid()
vamos_matroid()
all_subsets_matroid(r::Int)
projective_plane(q::Int)
projective_geometry(r::Int, q::Int; check::Bool=false)
affine_geometry(r::Int, q::Int; check::Bool=false)
```

### Modifying matroids
```@docs
dual_matroid(M::Matroid)
direct_sum(M::Matroid, N::Matroid)
deletion(M::Matroid,set::GroundsetType)
restriction(M::Matroid, set::GroundsetType)
contraction(M::Matroid,set::GroundsetType)
minor(M::Matroid, set_del::GroundsetType, set_cont::GroundsetType)
principal_extension(M::Matroid, set::GroundsetType, elem::ElementType)
free_extension(M::Matroid, elem::ElementType)
series_extension(M::Matroid, old::ElementType, new::ElementType)
parallel_extension(M::Matroid, old::ElementType, new::ElementType)
```

## Properties
```@docs
matroid_groundset(M::Matroid)
length(M::Matroid)
rank(M::Matroid)
rank(M::Matroid, set::GroundsetType)
bases(M::Matroid)
nonbases(M::Matroid)
circuits(M::Matroid)
hyperplanes(M::Matroid)
flats(M::Matroid, r::Union{Int,Nothing}=nothing)
cyclic_flats(M::Matroid, r::Union{Int,Nothing}=nothing)
closure(M::Matroid, set::GroundsetType)
nullity(M::Matroid, set::GroundsetType)
fundamental_circuit(M::Matroid, basis::GroundsetType, elem::ElementType)
fundamental_cocircuit(M::Matroid, basis::GroundsetType, elem::ElementType)
independent_sets(M::Matroid)
spanning_sets(M::Matroid)
cobases(M::Matroid)
cocircuits(M::Matroid)
cohyperplanes(M::Matroid)
corank(M::Matroid, set::GroundsetType)
is_clutter(sets::AbstractVector{T}) where T <: GroundsetType
is_regular(M::Matroid)
is_binary(M::Matroid)
is_ternary(M::Matroid)
n_connected_components(M::Matroid)
connected_components(M::Matroid)
is_connected(M::Matroid)
loops(M::Matroid)
coloops(M::Matroid)
is_loopless(M::Matroid)
is_coloopless(M::Matroid)
is_simple(M::Matroid)
direct_sum_components(M::Matroid)
connectivity_function(M::Matroid, set::GroundsetType)
is_vertical_k_separation(M::Matroid,k::IntegerUnion, set::GroundsetType) 
is_k_separation(M::Matroid,k::IntegerUnion, set::GroundsetType)
vertical_connectivity(M::Matroid)
girth(M::Matroid, set::GroundsetType=M.groundset)
tutte_connectivity(M::Matroid)
tutte_polynomial(M::Matroid)
characteristic_polynomial(M::Matroid)
reduced_characteristic_polynomial(M::Matroid)
revlex_basis_encoding(M::Matroid)
is_isomorphic(M1::Matroid, M2::Matroid)
is_minor(M::Matroid, N::Matroid)
```

### Chow Rings
```@docs
chow_ring(M::Matroid; ring::Union{MPolyRing,Nothing}=nothing, extended::Bool=false)
augmented_chow_ring(M::Matroid)
```

## Matroid strata and realization spaces

For a matroid ``M``, of rank ``d`` on ``\{1,\ldots,n\}`` realizable over a ring ``K``, 
its *matroid realization space* ``\mathsf{Gr}(M)`` is the locally closed subscheme of the Grassmannian
``\mathsf{Gr}(d,K^n)`` of ``d``-dimensional subspaces of ``K^n`` realizing ``M``. Precisely, 

\begin{equation}
    \mathsf{Gr}(M) = \{F \in \mathsf{Gr}(d,K^n) \, : \, p_{I}(F) \neq 0 \text{ iff } I \in \mathcal(B)(M)\}
\end{equation}

where ``\mathcal{B}(M)`` denotes the set of bases of ``M`` and  ``p_{I}(F)`` denotes the ``I``th Pl\"ucker 
coordinate of the linear subspace ``F \subset K^n``. The coordinate ring of ``\mathsf{Gr}(M)`` can be computed
using matrix coordinates, see Construction 2.2 of [Cor21](@cite). The following function computes this ring. 

```@docs
matroid_stratum_matrix_coordinates(M::Matroid, B::Vector{Int64}, F::AbstractAlgebra.Ring = ZZ)
```

When the matroid ``M`` is connected, the diagonal torus of ``\mathsf{PGL}(n)`` acts freely on 
``\mathsf{Gr}(M)`` and its quotient is the *realization space* ``R(M)`` of ``M``. There are two main differences between the coordinate rings ``K[Gr(M)]`` and ``K[R(M)]``.
 
 * To comupte ``K[R(M)]``, we may assume that ``d+1`` columns ``A = \{a_1,\ldots,a_{d+1}\}``  of the reference matrix ``X`` are the unit vectors ``e_1,\ldots,e_n`` and ``e_1+ \cdots + e_n``. Note that every ``d`` element subset of ``A`` must be a basis of ``M``, and the existence of such a ``A`` follows from the connectedness hypothesis on ``M``. 
 * The columns of the reference matrix ``X`` may be treated as projective coordiantes.

The coordinate ring of ``R(M)`` is computed in the following function. 


```@docs
matroid_realization_space(M::Matroid, A::Vector{Int64}, F::AbstractAlgebra.Ring=ZZ)
```
