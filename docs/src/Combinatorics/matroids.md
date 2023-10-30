```@meta
CurrentModule = Oscar
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
bond_matroid(g::Graph)
cocycle_matroid(g::Graph)
cycle_matroid
Matroid(pm_matroid::Polymake.BigObjectAllocated, E::GroundsetType=Vector{Integer}(1:pm_matroid.N_ELEMENTS))
matroid_from_bases(bases::AbstractVector{T}, nelements::IntegerUnion; check::Bool=true) where T<:GroundsetType
matroid_from_circuits(circuits::AbstractVector{T}, nelements::IntegerUnion) where T<:GroundsetType
matroid_from_hyperplanes(hyperplanes::AbstractVector{T}, nelements::IntegerUnion) where T<:GroundsetType
matroid_from_matrix_columns(A::MatrixElem; check::Bool=true)
matroid_from_matrix_rows(A::MatrixElem, ; check::Bool=true)
matroid_from_nonbases(nonbases::AbstractVector{T}, nelements::IntegerUnion; check::Bool=true) where T<:GroundsetType
matroid_from_revlex_basis_encoding(rvlx::String, r::IntegerUnion, n::IntegerUnion)
```

## Examples

```@docs
affine_geometry(r::Int, q::Int; check::Bool=false)
all_subsets_matroid(r::Int)
fano_matroid()
moebius_kantor_matroid()
non_fano_matroid()
non_pappus_matroid()
pappus_matroid()
perles_matroid()
projective_geometry(r::Int, q::Int; check::Bool=false)
projective_plane(q::Int)
R10_matroid()
uniform_matroid(r::IntegerUnion,n::IntegerUnion)
vamos_matroid()
```

### Modifying matroids
```@docs
contraction(M::Matroid,set::GroundsetType)
deletion(M::Matroid,set::GroundsetType)
direct_sum(M::Matroid, N::Matroid)
dual_matroid(M::Matroid)
free_extension(M::Matroid, elem::ElementType)
minor(M::Matroid, set_del::GroundsetType, set_cont::GroundsetType)
parallel_extension(M::Matroid, old::ElementType, new::ElementType)
principal_extension(M::Matroid, set::GroundsetType, elem::ElementType)
restriction(M::Matroid, set::GroundsetType)
series_extension(M::Matroid, old::ElementType, new::ElementType)
```

## Properties
```@docs
automorphism_group(M::Matroid)
bases(M::Matroid)
characteristic_polynomial(M::Matroid)
circuits(M::Matroid)
closure(M::Matroid, set::GroundsetType)
cobases(M::Matroid)
cocircuits(M::Matroid)
cohyperplanes(M::Matroid)
coloops(M::Matroid)
connected_components(M::Matroid)
connectivity_function(M::Matroid, set::GroundsetType)
corank(M::Matroid, set::GroundsetType)
cyclic_flats(M::Matroid, r::Union{Int,Nothing}=nothing)
direct_sum_components(M::Matroid)
flats(M::Matroid, r::Union{Int,Nothing}=nothing)
fundamental_circuit(M::Matroid, basis::GroundsetType, elem::ElementType)
fundamental_cocircuit(M::Matroid, basis::GroundsetType, elem::ElementType)
girth(M::Matroid, set::GroundsetType=M.groundset)
hyperplanes(M::Matroid)
independent_sets(M::Matroid)
is_binary(M::Matroid)
is_clutter(sets::AbstractVector{T}) where T <: GroundsetType
is_coloopless(M::Matroid)
is_connected(M::Matroid)
is_isomorphic(M1::Matroid, M2::Matroid)
is_k_separation(M::Matroid,k::IntegerUnion, set::GroundsetType)
is_loopless(M::Matroid)
is_minor(M::Matroid, N::Matroid)
is_regular(M::Matroid)
is_simple(M::Matroid)
is_ternary(M::Matroid)
is_vertical_k_separation(M::Matroid,k::IntegerUnion, set::GroundsetType) 
length(M::Matroid)
loops(M::Matroid)
matroid_base_polytope(M::Matroid)
matroid_groundset(M::Matroid)
n_connected_components(M::Matroid)
nonbases(M::Matroid)
nullity(M::Matroid, set::GroundsetType)
rank(M::Matroid)
rank(M::Matroid, set::GroundsetType)
reduced_characteristic_polynomial(M::Matroid)
revlex_basis_encoding(M::Matroid)
spanning_sets(M::Matroid)
tutte_connectivity(M::Matroid)
tutte_polynomial(M::Matroid)
vertical_connectivity(M::Matroid)
```


### Chow Rings
```@docs
chow_ring(M::Matroid; ring::Union{MPolyRing,Nothing}=nothing, extended::Bool=false)
augmented_chow_ring(M::Matroid)
```

## Matroid realization spaces


Let ``M`` be a matroid of rank ``d`` on a ground set ``E`` of size ``n``. Its *realization space* 
``\mathcal{R}(M)`` is an affine scheme that parameterizes all hyperplane arrangements that realize
the matroid ``M`` (up to the action of ``PGL(r)``).  We provide functions that determine the 
affine coordinate ring of ``\mathcal{R}(M)``. 


```@docs
is_realizable(M::Matroid; char::Union{Int,Nothing}=nothing, q::Union{Int,Nothing}=nothing)
defining_ideal(RS::MatroidRealizationSpace)
inequations(RS::MatroidRealizationSpace)
ambient_ring(RS::MatroidRealizationSpace)
realization_space(M::Matroid; B::Union{GroundsetType,Nothing} = nothing, 
    F::AbstractAlgebra.Ring = ZZ, saturate::Bool=false, reduce::Bool=true,
    char::Union{Int,Nothing}=nothing, q::Union{Int,Nothing}=nothing)
realization(M::Matroid; B::Union{GroundsetType,Nothing} = nothing, 
    F::AbstractAlgebra.Ring = ZZ, saturate::Bool=false, reduce::Bool=true,
    char::Union{Int,Nothing}=nothing, q::Union{Int,Nothing}=nothing)
realization(RS::MatroidRealizationSpace)
```

If ``B`` is the polynomial ring `ambient_ring(RS)`, ``I`` the ideal `defining_ideal(RS)`, and 
``U`` the multiplicative semigroup generated by `inequations(RS)`, then the coordinate ring of
the realization space ``\mathcal{R}(M)`` is isomorphic to ``U^{-1}B/I``.  
