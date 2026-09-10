module InjectiveResolutions
using ..Oscar
using ..Oscar: IntegerUnion # add other things that are not exported from Oscar here
# functions with new methods
import ..Oscar:
  _build_sparse_row,
  _degree_fast,
  _extend_free_resolution,
  _graded_kernel,
  _reduce,
  _saturation,
  annihilator,
  coefficient_ring,
  coefficients,
  cone,
  coordinates,
  coordinates_atomic,
  coordinates_via_transform,
  degree,
  dim,
  elem_type,
  evaluate,
  faces,
  free_resolution,
  gens,
  grading_group,
  hyperplanes,
  images_of_generators,
  in_atomic,
  intersect,
  inv,
  is_normal,
  is_pointed,
  is_subset,
  is_zm_graded,
  kernel,
  kernel_atomic,
  lift_std,
  ModuleGens,
  normal_form,
  one,
  oscar_free_module,
  oscar_generators,
  primitive_generator,
  singular_freemodule,
  singular_generators,
  singular_module,
  singular_poly_ring,
  sparse_row,
  standard_basis,
  SubModuleOfFreeModule,
  syzygy_module,
  twist,
  zero,
  zonotope

import ..Oscar.Singular:
  FreeModule,
  has_global_ordering,
  svector,
  Module


for i in names(Oscar)
  !isdefined(Oscar, i) && continue
  @eval import Oscar: $i
  #@eval export $i
end

for i in names(Oscar.Orderings)
  !isdefined(Oscar.Orderings, i) && continue
  @eval import Oscar.Orderings: $i
  #@eval export $i
end

#=
for i in names(Oscar; all=true)
  !isdefined(Oscar, i) && continue
  @eval import Hecke: $i
end
=#

import Base:
  +,
  -,
  *,
  ==,
  deepcopy_internal
# add more things here

## Functions visible on the outside
export monoid_algebra
export monoid_algebra_ideal
export faces
export hyperplanes
export saturation

export irreducible_resolution
export irreducible_decomposition

export injective_resolution
export old_injective_resolution

export local_cohomology
export local_cohomology_all
export zeroth_local_cohomology

export MonoidAlgebra
export MonoidAlgebraIdeal
export MonoidAlgebraElem
export AffineSemigroup
export affine_semigroup
export ambient_dimension
export semigroup_generators
export polyhedral_cone
export cone
export saturation_ideal
export saturation_map
export holes_module
export is_Q_graded

export compute_shift
export compute_shift_bound
export InjMod
export IndecInj
export monomial_matrix
export irreducible_hull
export kQ_module
export injective_hull
export Q_graded_part
export mod_quotient
export monoid_algebra_ideal
export FaceQ
export prime_of_face
export generators_W_H
export ZF_basis
export generates_Zd
export degrees_of_bass_numbers
export degrees_of_bass_numbers_bound
export graded_bass_numbers
export is_minimal
export in_intersection
export in_semigroup
export is_in_aZF
export coefficients_wrt_generators
export relevant_generators
export relevant_relations
export _get_irreducible_ideal
export underlying_element
export mod_saturate
export underlying_ideal
export _get_irreducible_ideal_unsaturated
export is_in_semigroup

#########################
# some composite types
#########################

include("MonoidAlgebra.jl")

struct IndecInj #indecomposable injective
  face::FaceQ
  vector::Vector{Int}
end
mutable struct InjMod #ZZ^d-graded injective module over monoid algebra
  monoid_algebra::MonoidAlgebra
  indec_injectives::Vector{IndecInj}
  Q_graded_part::Union{SubquoModule,Nothing}

  function InjMod(A::MonoidAlgebra,I::Vector{IndecInj})
    return new(A,I,nothing)
  end
end

mutable struct IrrSum #direct sum of modules k[Q]/W, where W is an irreducible ideal
  monoid_algebra::MonoidAlgebra
  indec_injectives::Vector{IndecInj}
  kQ_module::Union{SubquoModule,Nothing} #k[Q]/W

  function IrrSum(A::MonoidAlgebra,I::Vector{IndecInj})
    return new(A,I,nothing)
  end
end

function kQ_module(I::IrrSum) #irreducible sum as finitely generated k[Q]-module
  if I.kQ_module === nothing
    I.kQ_module = compute_Q_graded_part(I.monoid_algebra, I.indec_injectives)
  end
  return I.kQ_module
end

function Q_graded_part(I::InjMod)
  if I.Q_graded_part === nothing
    I.Q_graded_part = compute_Q_graded_part(I.monoid_algebra, I.indec_injectives)
  end
  return I.Q_graded_part
end

function compute_Q_graded_part(kQ::MonoidAlgebra, I::Vector{IndecInj})
  if isempty(I)
    F = graded_free_module(kQ, 0)
    return quo(F, [zero(F)])[1]
  end
  if is_normal(kQ)
    irreducible_ideals = [_get_irreducible_ideal(kQ, J) for J in I]
  else
    irreducible_ideals = [_get_irreducible_ideal_unsaturated(kQ, J) for J in I]
  end
  irreducible_comp = [quotient_ring_as_module(Ji) for Ji in irreducible_ideals]

  isone(length(irreducible_comp)) && return irreducible_comp[1]
  return direct_sum(irreducible_comp..., task=:none)
end

struct IrrRes #Q-graded irreducible resolution
  mod::SubquoModule
  irr_sums::Vector{IrrSum}
  cochain_maps::Vector{SubQuoHom}
  cochain_complex::ComplexOfMorphisms # if sequence not exact return trivial cochain_complex (M0 -> M0)
end

#direct sum of two injective modules over the same monoid algebra
function +(I::InjMod, J::InjMod)
  @req I.monoid_algebra == J.monoid_algebra "monoid algebras not the same"
  return InjMod(I.monoid_algebra, vcat(I.indec_injectives, J.indec_injectives))
end

#direct sum of two injective modules over the same monoid algebra
function +(I::IrrSum, J::IrrSum)
  @req I.monoid_algebra == J.monoid_algebra "monoid algebras not the same"
  return IrrSum(I.monoid_algebra, vcat(I.indec_injectives, J.indec_injectives))
end

struct InjRes #ZZ^d-graded injective resolution
  mod::SubquoModule
  inj_mods::Vector{InjMod}
  cochain_maps::Vector{MatElem}
  upto::Int
  Q_graded_part::IrrRes
  shift::Vector{Int} #not needed
end

struct MonomialMatrix{T <: Union{InjMod, IrrSum}} # monomial matrix as in [HM05]
  matrix::MatElem
  source::T
  target::T
  index::Int  # cohomological index i of the differential d^i
end

function Base.show(io::IO, mm::MonomialMatrix{InjMod})
  print(io, "monomial matrix for I^", mm.index, " -> I^", mm.index + 1)
end

function Base.show(io::IO, mm::MonomialMatrix{IrrSum})
  print(io, "monomial matrix for W^", mm.index, " -> W^", mm.index + 1)
end

function _show_indec(io::IO, J::IndecInj, ::Type{InjMod})
  println(io, "  k{", J.vector, " + F - Q}, where p_F = ", J.face.prime)
end
function _show_indec(io::IO, J::IndecInj, ::Type{IrrSum})
  println(io, "  k{", J.vector, " + F - Q}_Q, where p_F = ", J.face.prime)
end

function Base.show(io::IO, ::MIME"text/plain", mm::MonomialMatrix{T}) where T
  src_label, tgt_label = T == InjMod ? ("I^$(mm.index)", "I^$(mm.index + 1)") : ("W^$(mm.index)", "W^$(mm.index + 1)")
  println(io, "monomial matrix for $src_label -> $tgt_label")
  println(io, "source summands ($src_label):")
  for J in mm.source.indec_injectives
    _show_indec(io, J, T)
  end
  println(io, "target summands ($tgt_label):")
  for J in mm.target.indec_injectives
    _show_indec(io, J, T)
  end
  print(io, mm.matrix)
end

function Base.show(io::IO,J::InjMod)
  print(
    io, "injective module given by direct sum of ", length(J.indec_injectives)," indecomposable injectives"
  )
end

function Base.show(io::IO, ::MIME"text/plain", J::InjMod)
  println(io, "injective module given by direct sum of indecomposable injectives")
  for Ji in J.indec_injectives
    _show_indec(io, Ji, InjMod)
  end
  print(io, "over ", J.monoid_algebra)
end

function Base.show(io::IO,W::IrrSum)
  print(io, "irreducible sum given by direct sum of ", length(W.indec_injectives)," components")
end

function Base.show(io::IO, ::MIME"text/plain", W::IrrSum)
  println(io, "irreducible sum given by direct sum of components")
  for Ji in W.indec_injectives
    _show_indec(io, Ji, IrrSum)
  end
  print(io, "over ", W.monoid_algebra)
end

function Base.show(io::IO, ::MIME"text/plain", res::InjRes)
  println(io, "injective resolution ")
  println(io, "  ", join(["I^$i" for i in 0:res.upto], " -> "))
  println(io, "where ")
  for i in eachindex(res.inj_mods)
    println(io, " I^$(i-1) = direct sum of")
    for Ji in res.inj_mods[i].indec_injectives
      print(io, "  "); _show_indec(io, Ji, InjMod)
    end
  end
  println(io, "of ", res.mod)
  print(io, "over ", base_ring(res.mod))
end

function Base.show(io::IO, ::MIME"text/plain", res::IrrRes)
  println(io, "irreducible resolution ")
  println(io, "  ", join(["W^$i" for i in 0:(length(res.irr_sums) - 1)], " -> "))
  println(io, "where ")
  for i in eachindex(res.irr_sums)
    println(io, " W^$(i-1) = direct sum of")
    for Ji in res.irr_sums[i].indec_injectives
      print(io, "  "); _show_indec(io, Ji, IrrSum)
    end
  end
  println(io, "of ", res.mod)
  print(io, "over ", base_ring(res.mod))
end

@doc raw"""
    monomial_matrix(i::Int, res::IrrRes)
    monomial_matrix(i::Int, res::InjRes)

Return the differential $d^i$ of the (co)chain complex of `res` as a
[`MonomialMatrix`](@ref): the matrix of the map together with its row and column
labels, i.e. the `(face, degree)` data of the source and target summands.

Here `i` is the cohomological degree of the differential
$d^i \colon W^i \to W^{i+1}$ (resp. $J^i \to J^{i+1}$).
"""
function monomial_matrix(i::Int, res::IrrRes)
  n = length(res.irr_sums)
  @req 0 <= i <= n - 2 "cohomological index i must be in 0:$(n - 2) for this resolution"
  A = matrix(res.cochain_maps[i + 2])
  src = res.irr_sums[i + 1]
  tgt = res.irr_sums[i + 2]
  @assert nrows(A) == length(src.indec_injectives) "row count $(nrows(A)) ≠ number of source summands $(length(src.indec_injectives))"
  @assert ncols(A) == length(tgt.indec_injectives) "column count $(ncols(A)) ≠ number of target summands $(length(tgt.indec_injectives))"
  return MonomialMatrix(A, src, tgt, i)
end

function monomial_matrix(i::Int, res::InjRes)
  maxi = min(length(res.cochain_maps), length(res.inj_mods)) - 2
  @req 0 <= i <= maxi "cohomological index i must be in 0:$maxi for this resolution"
  A = res.cochain_maps[i + 2]
  src = res.inj_mods[i + 1]
  tgt = res.inj_mods[i + 2]
  @assert nrows(A) == length(src.indec_injectives) "row count $(nrows(A)) ≠ number of source summands $(length(src.indec_injectives))"
  @assert ncols(A) == length(tgt.indec_injectives) "column count $(ncols(A)) ≠ number of target summands $(length(tgt.indec_injectives))"
  return MonomialMatrix(A, src, tgt, i)
end

function Base.show(io::IO, ::MIME"text/plain", Ji::IndecInj)
  println(io, "indecomposable injective")
  println(io, "  k{", Ji.vector, " + F - Q},")
  print(io, "where p_F = ", Ji.face.prime)
end

function Base.show(io::IO, Ji::IndecInj)
  print(
      io, "indecomposable injective k{", Ji.vector, " + F - Q}, where p_F = ", Ji.face.prime
    )
end

@doc raw"""
    generators_W_H(kQ::MonoidAlgebra, H::HyperplaneQ, a::Vector{Int})

Given a monoid algebra  kQ = $k[Q]$, a hyperplane $H$ that bounds the polyhedral cone $\RR_{\geq 0}Q$ and a vector
$a \in \mathbb{Z}^d$, return a finite set $B$ such that

$(x^b \mid b \in B) \cong k\{(a + H_+^\circ)\cap Q\}.$

This is Algorithm 3.11. in [HM05](@cite).

!!! note
    The monoid algebra $k[Q]$ must be normal. 
"""
function generators_W_H(kQ::MonoidAlgebra, H::HyperplaneQ, a::Vector{Int})
  @assert torsion_free_rank(grading_group(kQ)) == length(a)
  @assert H in hyperplanes(kQ)

  F = intersect(H.hyperplane, cone(kQ))

  #get faces of Q intersecting F only origin in Q
  D = Vector{Polyhedron}()
  for f in faces(kQ)
    if dim(intersect(f.poly, F)) == 0
      push!(D, f.poly)
    end
  end

  B = Vector{Vector{Int}}()
  PaF = polyhedron(H.A, H.A*a+H.b) #a + RR h
  for d in D
    I = intersect(PaF, d) #(a + RR H)\cap RR_+D
    if dim(I) >= 0
      B_d = [
        a for
        a in lattice_points(I + zonotope(kQ)[1]) if (dim(intersect(convex_hull(a), PaF)) < 0)
      ]
      append!(B, B_d)
    end
  end
  return B
end

@doc raw"""
    degrees_of_bass_numbers(M::SubquoModule,i::Int)

Return the $\mathbb{Z}^d$-degrees of non-zero Bass numbers of $M$ up to cohomological degree $i$.
"""
function degrees_of_bass_numbers(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)
  R_Q = base_ring(M)

  # residue field
  I_m = ideal(R_Q, gens(R_Q))
  k = quotient_ring_as_module(I_m)

  # compute free resolution of k, apply hom(-,M) and take homology to obtain Ext^j(k,M)
  free_res = free_resolution(k; length=i+2)
  lifted   = hom(free_res.C[first(Hecke.map_range(free_res.C)):-1:1], M)

  degrees = Vector{Vector{Int}}()
  for j in 0:i
    E = simplify_light(homology(lifted, j))[1]
    for g in filter(!is_zero, gens(E))
      push!(degrees, degree(Vector{Int}, g))
    end
  end
  return unique(degrees) # filter duplicates
end

@doc raw"""
    degrees_of_bass_numbers_bound(M::SubquoModule, i::Int)

Return a finite set `D` of $\mathbb{Z}^d$-degrees such that every degree of a
non-zero Bass number of `M` up to cohomological degree `i` lies in `D + Q`:

$$D = \{\deg(g) - \deg(e) : g \in \mathrm{gens}(M),\ e \in \mathrm{gens}(F_j),\ 0 \le j \le i+1\},$$

where $F_\bullet$ is a graded free resolution of the residue field $k = R_Q/m$.

Since $\mathrm{Ext}^j(k,M)$ is a subquotient of $\mathrm{Hom}(F_j, M) = \bigoplus_e M(\deg e)$
and $\mathrm{supp}(M) \subseteq \bigcup_g (\deg(g) + Q)$, the degrees of the Bass
numbers are contained in `D + Q`. This is a cheap over-approximation of
[`degrees_of_bass_numbers`](@ref): it reads degrees off the resolution and avoids
applying $\mathrm{Hom}(-, M)$ and computing homology, at the cost of a possibly
larger set. Because $Q$ is a semigroup, any shift `a` with `D + a ⊆ Q` also
satisfies `(bass degrees) + a ⊆ Q`, so `D` may replace the exact Bass degrees
in [`compute_shift`](@ref); see [`compute_shift_bound`](@ref).
"""
function degrees_of_bass_numbers_bound(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)
  R_Q = base_ring(M)

  # residue field k = R_Q/m and its graded free resolution
  k = quotient_ring_as_module(ideal(R_Q, gens(R_Q)))
  free_res = free_resolution(k; length=i+2)

  # generator degrees of M
  deg_M = [degree(Vector{Int}, g) for g in filter(!is_zero, gens(M))]

  # generator degrees of F_0, ..., F_{i+1}
  deg_F = Vector{Vector{Int}}()
  for j in 0:(i + 1)
    Fj = free_res[j]
    for d in degrees_of_generators(Fj)
      push!(deg_F, Int[d[l] for l in 1:ngens(parent(d))])
    end
  end

  # D = { deg(g) - deg(e) }
  D = Set{Vector{Int}}()
  for dg in deg_M, de in deg_F
    push!(D, dg .- de)
  end
  return collect(D)
end

# Given a finite set of ℤ^d-degrees, return a shift a = j*c (j ≥ 0, c the sum of
# primitive ray generators of Q) such that every degree + a lies in Q.
function _shift_into_Q(kQ::MonoidAlgebra, degrees::Vector{Vector{Int}})
  Q = affine_semigroup(kQ)
  if is_normal(kQ)
    c = zonotope(kQ)[2]
    C = cone(kQ)
    in_Q = b -> is_subset(convex_hull(b), C) # for normal Q: Q = cone ∩ lattice
  else
    c = sum(gens(Q)) #TODO: is the best way?
    in_Q = b -> is_in_semigroup(kQ, b)
  end

  j = 0
  for b in degrees
    j_b = 0
    v = b
    while !in_Q(v)
      v = v + c
      j_b += 1
    end
    j = max(j, j_b)
  end
  return j*c
end

@doc raw"""
    compute_shift(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)

Let $M$ be finitely generated $\mathbb{Z}^d$-graded module over a monoid algebra $k[Q]$. This function computes $a\in \mathbb{Z}^d$
such that all $\mathbb{Z}^d$-degrees of non-zero Bass numbers of $M(-a)$ lie in $Q$.
"""
function compute_shift(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)
  #get all degrees of non-zero Bass numbers up to cohomological degree i
  return _shift_into_Q(base_ring(M), degrees_of_bass_numbers(M, i))
end

@doc raw"""
    compute_shift_bound(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)

Like [`compute_shift`](@ref), but uses the cheap over-approximation
[`degrees_of_bass_numbers_bound`](@ref) instead of the exact Bass-number degrees.
The returned shift is valid (all Bass degrees of `M(-a)` lie in `Q`) but may be
larger than the one from `compute_shift`.
"""
function compute_shift_bound(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)
  return _shift_into_Q(base_ring(M), degrees_of_bass_numbers_bound(M, i))
end

# MILP that minimises sum(a) subject to a + B[j] ∈ Q for all j, where a ≥ 0.
# Variables: [a (d); λ_0 (n); λ_1 (n); ...; λ_k (n)] with λ_i integer for non-normal Q.
# For normal Q, λ_i can stay continuous (cone ∩ ℤ^d = Q), giving a pure LP.
function _shift_milp(kQ::MonoidAlgebra, degrees::Vector{Vector{Int}})
  isempty(degrees) && return zeros(Int, ambient_dimension(kQ))
  A  = Matrix{Int64}(semigroup_generators(kQ))
  d, n = size(A)
  B  = [Vector{Int64}(b) for b in degrees]
  k  = length(B)

  DA  = block_diagonal_matrix(ZZ, [A for _ in 1:(k+1)])
  I_d = identity_matrix(ZZ, d)
  N   = hcat(vcat([-I_d for _ in 1:(k+1)]...), DA)
  M   = vcat(N, -N, -identity_matrix(ZZ, (k+1)*n + d))
  b   = vcat(B...)
  M_b = vcat(zeros(Int64, d), b, zeros(Int64, d), -b, zeros(Int64, (k+1)*n + d))

  P    = polyhedron(M, M_b)
  c    = vcat(ones(Int64, d), zeros(Int64, (k+1)*n))
  # For non-normal Q, make λ_i integer so we enforce semigroup (not just cone) membership.
  int_vars = is_normal(kQ) ? Int[] : collect(d+1 : d+(k+1)*n)
  milp = mixed_integer_linear_program(P, c; integer_variables=int_vars, convention = :min)
  _, opt = solve_milp(milp)
  opt === nothing && error("MILP shift infeasible — no shift puts all Bass-number degrees into Q")
  return [ceil(Int, opt[i]) for i in 1:d]
end

@doc raw"""
    compute_shift_milp(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)

Like [`compute_shift`](@ref), but computes the shift by minimising `sum(a)`
subject to `a + b ∈ Q` for every Bass-number degree `b` (LP for normal `Q`,
MILP with integer semigroup-membership variables for non-normal `Q`).
This can give a smaller shift than `compute_shift` when the optimal direction
is not a multiple of the sum of primitive ray generators.
"""
function compute_shift_milp(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)
  return _shift_milp(base_ring(M), degrees_of_bass_numbers(M, i))
end

@doc raw"""
    compute_shift_milp_bound(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)

Like [`compute_shift_milp`](@ref), but uses the cheap over-approximation
[`degrees_of_bass_numbers_bound`](@ref) for the Bass-number degrees.
"""
function compute_shift_milp_bound(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)
  return _shift_milp(base_ring(M), degrees_of_bass_numbers_bound(M, i))
end

@doc raw"""
    mod_quotient(M::SubquoModule,I::Ideal)

Computes the submodule

$(0 :_M I) := \{m \in M \mid m\cdot I = 0\}.$

# Examples
```jldoctest
julia> R_Q,(x,y) = graded_polynomial_ring(QQ,["x","y"]; weights = [[1,0],[0,1]])
(Graded multivariate polynomial ring in 2 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y])

julia> I = ideal(R_Q,[x^4,x^2*y^2,y^4])
Ideal generated by
  x^4
  x^2*y^2
  y^4

julia> M = quotient_ring_as_module(I)
Graded subquotient of graded submodule of R_Q^1 with 1 generator
  1: e[1]
by graded submodule of R_Q^1 with 3 generators
  1: x^4*e[1]
  2: x^2*y^2*e[1]
  3: y^4*e[1]

julia> m = ideal(R_Q,[x,y])
Ideal generated by
  x
  y

julia> Oscar.InjectiveResolutions.mod_quotient(M,m)
(Graded subquotient of graded submodule of R_Q^1 with 2 generators
  1: x*y^3*e[1]
  2: x^3*y*e[1]
by graded submodule of R_Q^1 with 3 generators
  1: x^4*e[1]
  2: x^2*y^2*e[1]
  3: y^4*e[1], Hom: graded subquotient of graded submodule of R_Q^1 with 2 generators
  1: x*y^3*e[1]
  2: x^3*y*e[1]
by graded submodule of R_Q^1 with 3 generators
  1: x^4*e[1]
  2: x^2*y^2*e[1]
  3: y^4*e[1] -> M)
```
"""
function mod_quotient(M::SubquoModule, I::Ideal)
  T = elem_type(M)
  R_I = quotient_ring_as_module(I)
  m = R_I[1] #generator of R_I
  H = hom(R_I, M)[1]

  if is_zero(H)
    Q_gens = Vector{T}()
  else
    Q_gens = [element_to_homomorphism(g)(m) for g in gens(H)]
  end
  return sub(M, Q_gens)
end

@doc raw"""
    mod_saturate(M::SubquoModule,I::Ideal)

Compute the saturation

$(0 :_M I^\infty) := \{m \in M \mid m\cdot I^n = 0\text{ for some }n\in \NN_{>0}\}.$

# Examples
```jldoctest
julia> R_Q,(x,y) = graded_polynomial_ring(QQ,["x","y"]; weights = [[1,0],[0,1]])
(Graded multivariate polynomial ring in 2 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y])

julia> I = ideal(R_Q,[x^4,x^2*y^2,y^4])
Ideal generated by
  x^4
  x^2*y^2
  y^4

julia> M = quotient_ring_as_module(I)
Graded subquotient of graded submodule of R_Q^1 with 1 generator
  1: e[1]
by graded submodule of R_Q^1 with 3 generators
  1: x^4*e[1]
  2: x^2*y^2*e[1]
  3: y^4*e[1]

julia> m = ideal(R_Q,[x,y])
Ideal generated by
  x
  y

julia> Oscar.InjectiveResolutions.mod_saturate(M,m)
Graded subquotient of graded submodule of R_Q^1 with 12 generators
  1: x*y^3*e[1]
  2: x^3*y*e[1]
  3: y^3*e[1]
  4: x*y^2*e[1]
  5: x^2*y*e[1]
  6: x^3*e[1]
  7: y^2*e[1]
  8: x*y*e[1]
  9: x^2*e[1]
  10: y*e[1]
  11: x*e[1]
  12: e[1]
by graded submodule of R_Q^1 with 3 generators
  1: x^4*e[1]
  2: x^2*y^2*e[1]
  3: y^4*e[1]
```
"""
function mod_saturate(M::SubquoModule, I::Ideal)
  M_sat, _ = mod_quotient(M, I)
  M_prev = M_sat #previous module quotient
  i = 2
  while true
    M_q, _ = mod_quotient(M, I^i)
    if M_prev == M_q
      break
    end
    M_sat = M_sat + M_q
    M_prev = M_q #update
    i = i + 1
  end
  return M_sat
end

@doc raw"""
    ZF_basis(M::SubquoModule, p::FaceQ)

Let $p = k\{Q\setminus F\}$ for some face $F$. This functions computes a $k[\mathbb{Z}F]$-basis of the quotient

$(0 :_M p)[\mathhbb{Z}F] = \{m \in M \mid m\cdot p = 0\}[\mathbb{Z}F].$
"""
function ZF_basis(N::SubquoModule{<:MonoidAlgebraElem}, p::FaceQ)
  kQ = base_ring(N)
  @assert kQ.algebra == base_ring(p.prime)

  T = elem_type(N)

  #compute quotient (0 :_M p)[\ZZ F]
  Np = mod_quotient(N, monoid_algebra_ideal(kQ,p.prime))[1]

  #initialize
  L = Np
  h_L = identity_map(L)
  B = Vector{T}() # empty vector of k[ZF]-basis 

  for g in filter(!is_zero, gens(Np))
    N_g = sub(L, [h_L(g)])[1] #submodule of N =(0 :_M p_F)/(y0,...,yn) generated by g

    if annihilator(N_g) == MonoidAlgebraIdeal(kQ, p.prime)
      push!(B, g)
    end
    N_B, _ = sub(Np, B)
    L, h_L = quo(Np, N_B) # update N
    if is_zero(L)
      break
    end
  end
  # return filter(!is_zero, B)
  return B
end

function evaluate(
    f::MonoidAlgebraElem{<:RingElem, PT}, 
    vals::Vector; 
    check::Bool=true
  ) where {PT <: MonoidAlgebra{<:RingElem, <:MPolyQuoRing}}
  return evaluate(underlying_element(f), vals; check)
end

function evaluate(a::MPolyQuoRingElem, vals::Vector; check::Bool=true)
  @check all(is_zero(evaluate(f, vals)) for f in gens(modulus(parent(a))))
  return evaluate(a.f, vals)
end

function evaluate(
    f::MonoidAlgebraElem{<:RingElem, PT}, 
    vals::Vector; 
    check::Bool=true
  ) where {PT <: MonoidAlgebra{<:RingElem, <:MPolyRing}}
  return evaluate(underlying_element(f), vals)
end

# get coefficients of m w.r.t. generators of M
function coefficients_wrt_generators(m::SubquoModuleElem{T}, N::SubquoModule{T}) where {T <: MonoidAlgebraElem}
  # @assert parent(m) == N "m is not an element of N"
  kQ = base_ring(N)
  k = coefficient_ring(kQ)
  m_amb = ambient_representative(m) ## why is this needed??
  coord_sparse = coordinates(N(m_amb))
  _coord = [evaluate(coord_sparse[i], [1 for _ in 1:ngens(kQ)]) for i in 1:ngens(N)]
  d = degree(Vector{Int}, m)
  for i in 1:length(_coord)
    if !is_zero(_coord[i])
      d_diff = d - degree(Vector{Int}, N[i])
      if !is_in_semigroup(kQ, d_diff) || is_zero(monomial_basis(kQ, d_diff)[1]*N[i])
        _coord[i] = k()
      end
    end
  end
  return _coord
end

#get all relevant generators of N, i.e., check (deg(b) + F) \cap (deg(g) + Q) ≠ ∅
function relevant_generators(N::SubquoModule{T}, p::FaceQ, b::SubquoModuleElem{T}) where {T <: MonoidAlgebraElem}
  kQ = base_ring(N)
  # b_poly = convex_hull(degree(Vector{Int},b))
  rel_gens = Vector{SubquoModuleElem}()
  for g in filter(!is_zero,gens(N))
    if in_intersection(kQ,degree(Vector{Int},g),degree(Vector{Int},b),p)
      push!(rel_gens,g)
    end 
  end
  return rel_gens
end

#get all b-relevant relations w.r.t. F
# check deg(r) + Q = b has a solution
function relevant_relations(N::SubquoModule{T},p::FaceQ, b::SubquoModuleElem{T}, c_b) where {T <: MonoidAlgebraElem}
  kQ = base_ring(N)
  R = relations(N)
  rel_rels = Vector{FreeModElem{elem_type(kQ)}}()
  for r in R
    if in_intersection(kQ,degree(Vector{Int},r),degree(Vector{Int},b),p)
      push!(rel_rels,r)
    end
  end
  return rel_rels
end


@doc raw"""
    coefficients(N::SubquoModule, p::FaceQ)

Return a subset Bp $\subseteq M$ and a $k$-matrix $\Lambda$ that defines an injective map

$(0 :_N p) \xrightarrow{\Lambda} \sum_{b\in Bp}k\{\deg(b) + F - Q\}.$

This fixes Algorithm 3.6. in [HM05](@cite).
"""
function coefficients(N::SubquoModule{T}, p::FaceQ) where {T <: MonoidAlgebraElem}
  kQ = base_ring(N)
  @assert base_ring(p.prime) == kQ.algebra

  k = coefficient_ring(kQ)

  # compute a k[ZF]-basis of (0 :_N p)[ZF]
  Bp = ZF_basis(N, p)
  if is_empty(Bp)
    return Bp, matrix(kQ, zeros(kQ, 1, 1))
  end

  if is_normal(kQ)
    return _coefficients_normal(N, p, Bp)
  else
    return _coefficients_non_normal(N, p, Bp)
  end
end

# Given, for a face p and a list of socle basis elements (their degrees `degs`,
# coefficient vectors `c_bs` w.r.t. the n generators of N, and per-element
# relation-coefficient matrices `C_rows`), compute the coefficient vectors
# lambda defining the embedding into the irreducible hull.
#
# Injectivity of the hull is a graded condition: the indecomposable injective
# k{deg(b)+F-Q} depends on deg(b) only modulo ZF, and the socle (0:_N p_F)
# decomposes over Z^d/ZF. So we work one ZF-class at a time. Within a class with
# socle elements b_1,...,b_s, all share the kernel K = ker(C) (relevant relations
# depend only on the class), and the c_{b_j}|_K are linearly independent. We pick
# lambda_1,...,lambda_s in K *dual* to the c_{b_j}, i.e. <lambda_i, c_{b_j}> =
# delta_{ij}, by solving one linear system. This makes the socle matrix the
# identity, hence the map injective. (Choosing each lambda_i only so that
# <lambda_i, c_{b_i}> != 0 -- one element at a time -- is NOT enough: the socle
# matrix can be singular even with linearly independent rows.)
function _dual_basis_lambdas(k, n::Int, p::FaceQ, degs, c_bs, C_rows)
  S = elem_type(k)
  lambda = Vector{Vector{S}}(undef, length(degs))

  # group socle elements by ZF-class of their degree
  groups = Vector{Tuple{Vector{Int}, Vector{Int}}}() # (representative degree, indices)
  for (idx, d_b) in enumerate(degs)
    gidx = findfirst(t -> is_in_aZF(t[1], p, d_b), groups)
    if gidx === nothing
      push!(groups, (d_b, Int[idx]))
    else
      push!(groups[gidx][2], idx)
    end
  end

  for (_, idxs) in groups
    s = length(idxs)
    # common kernel K = ker(C) for the class (rows from all its elements)
    rows = Vector{Vector{S}}()
    for i in idxs
      append!(rows, C_rows[i])
    end
    Kmat = isempty(rows) ? identity_matrix(k, n) : kernel(matrix(k, hcat(rows...)))
    m = nrows(Kmat)
    # P[j,t] = <c_{b_j}, kappa_t>, where kappa_t is the t-th row of Kmat
    P = matrix(k, s, m, S[sum(c_bs[idxs[j]][g] * Kmat[t, g] for g in 1:n) for j in 1:s for t in 1:m])
    # X (m x s) with P*X = I_s; rows of K^T*X = lambda are dual to the c_{b_j}
    X = Oscar.solve(P, identity_matrix(k, s); side = :right)
    Lambda = transpose(X) * Kmat # s x n
    for (jj, i) in enumerate(idxs)
      lambda[i] = S[Lambda[jj, g] for g in 1:n]
    end
  end
  return lambda
end

# Normal case: uses polyhedral intersection for relevance checks
function _coefficients_normal(N::SubquoModule{T}, p::FaceQ, Bp) where {T <: MonoidAlgebraElem}
  kQ = base_ring(N)
  k = coefficient_ring(kQ)
  R = relations(N)
  rel_data = [(r, convex_hull(degree(Vector{Int}, r)) + cone(kQ)) for r in R]

  degs = Vector{Vector{Int}}()
  c_bs = Vector{Vector{elem_type(k)}}()
  C_rows = Vector{Vector{Vector{elem_type(k)}}}()
  for b in Bp
    b_p = convex_hull(degree(Vector{Int}, b))
    b_ZF = b_p + p.poly + (-1)*p.poly # b + ZF
    #get coefficient vector w.r.t. generators of N
    b_amb = ambient_representative(b)
    _c_b = coordinates(N(b_amb))
    @assert all(is_zero(_c_b[i]) || is_homogeneous(_c_b[i]) for i in 1:ngens(N)) "non-homogeneous coordinate — all-ones evaluation is wrong here"
    c_b = [evaluate(_c_b[i], [1 for _ in 1:ngens(kQ)]) for i in 1:ngens(N)]

    #get all relevant generators of N, i.e., check (deg(b) + ZF) \cap (deg(g) + Q) ≠ ∅
    G_b = Vector{SubquoModuleElem}()
    G_b_indices = Vector{Int}()  # original generator indices into gens(N)
    for (g_idx, g_N) in enumerate(gens(N))
      is_zero(g_N) && continue
      g_p = convex_hull(degree(Vector{Int}, g_N))
      if dim(intersect(b_p + p.poly + (-1)*p.poly, g_p + cone(kQ))) >= 0
        push!(G_b, g_N)
        push!(G_b_indices, g_idx)
      end
    end
    x_Gb = [monomial_basis(kQ, degree(g))[1] for g in G_b]
    _N = sub(ambient_free_module(N), [ambient_representative(g) for g in G_b])[1]

    #get all b-relevant relations w.r.t. F
    C_bF = Vector{Vector{elem_type(k)}}()
    for (r,r_poly) in rel_data
      if dim(intersect(b_ZF, r_poly)) >= 0
        x_r = monomial_basis(kQ, degree(r))[1]
        a = lcm(x_Gb..., x_r)
        _r = (a//x_r).num*r

        if _r in _N
          _c_r = coordinates(_N(_r))
          c_r = Vector{elem_type(k)}()
          for i in 1:ngens(N)
            j = findfirst(==(i), G_b_indices)
            if j !== nothing
              push!(c_r, evaluate(_c_r[j], [1 for _ in 1:ngens(kQ)]))
            else
              push!(c_r, k())
            end
          end
          push!(C_bF, c_r)
        end
      end
    end

    push!(degs, degree(Vector{Int}, b))
    push!(c_bs, c_b)
    push!(C_rows, C_bF)
  end

  lambda = _dual_basis_lambdas(k, ngens(N), p, degs, c_bs, C_rows)
  return Bp, matrix(kQ, map(kQ, hcat(lambda...)))
end

# Non-normal case: uses MILP-based intersection checks
function _coefficients_non_normal(N::SubquoModule{T}, p::FaceQ, Bp) where {T <: MonoidAlgebraElem}
  kQ = base_ring(N)
  k = coefficient_ring(kQ)

  degs = Vector{Vector{Int}}()
  c_bs = Vector{Vector{elem_type(k)}}()
  C_rows = Vector{Vector{Vector{elem_type(k)}}}()
  # use saturation for monomial_basis since degrees may not be in Q
  kQsat = saturation(kQ)
  for b in Bp
    #get coefficient vector w.r.t. generators of N
    c_b = coefficients_wrt_generators(b, N)

    #get all relevant generators of N, i.e., check (deg(b) + ZF) \cap (deg(g) + Q) ≠ ∅
    rel_gens = relevant_generators(N, p, b)
    @assert !is_empty(rel_gens) "there are no relevant generators for b = $b and face p = $p"

    _N = sub(ambient_free_module(N), [ambient_representative(g) for g in rel_gens])[1]

    #get all b-relevant relations w.r.t. F
    rel_rels = relevant_relations(N, p, b, c_b)

    C_bF = Vector{Vector{elem_type(k)}}()
    x_rel_gens = [monomial_basis(kQsat, degree(g))[1] for g in rel_gens]
    for r in rel_rels
      x_r = monomial_basis(kQsat, degree(r))[1]
      a = lcm(x_rel_gens..., x_r)

      # compute scaling factor (a//x_r) and pull back to kQ
      scaling = divexact(a, x_r)
      # scale r: need to work in kQ, so express scaling as element of kQ if possible
      d_scaling = degree(Vector{Int}, scaling)
      if !is_in_semigroup(kQ, d_scaling) && !is_zero(d_scaling)
        continue
      end
      if is_zero(d_scaling)
        _r = r
      else
        _r = monomial_basis(kQ, d_scaling)[1] * r
      end

      if _r in _N
        _c_r = Oscar.coordinates(_r)
        c_r = Vector{elem_type(k)}()
        for i in 1:ngens(N)
          j = findfirst(g -> ambient_representative(g) == ambient_representative(N[i]), gens(N))
          if j == i
            push!(c_r, evaluate(_c_r[i], [1 for _ in 1:ngens(kQ)]))
          else
            push!(c_r, k())
          end
        end
        if !is_zero(c_r)
          push!(C_bF, c_r)
        end
      end
    end

    push!(degs, degree(Vector{Int}, b))
    push!(c_bs, c_b)
    push!(C_rows, C_bF)
  end

  if is_empty(Bp)
    return Bp, matrix(kQ, zeros(kQ, 1, 1))
  end
  lambda = _dual_basis_lambdas(k, ngens(N), p, degs, c_bs, C_rows)
  return Bp, matrix(kQ, map(kQ, hcat(lambda...)))
end

@doc raw"""
    irreducible_hull(Mi::SubquoModule, kQ::MonoidAlgebra, j=0)

Return an irreducible hull of $M$.
"""
function irreducible_hull(Mi::SubquoModule{<:MonoidAlgebraElem}, j=0)
  kQ = base_ring(Mi)
  T = elem_type(kQ)
  zero_Mi = (ideal(kQ, []) * Mi)[1]

  #initialize
  N = Mi
  summands = Vector{IndecInj}()
  lambda = Vector{dense_matrix_type(T)}()

  P = faces(kQ)
  for p in P
    Bp, lambda_p = coefficients(N, p)
    
    for b in Bp
      push!(summands, IndecInj(p, degree(Vector{Int}, b)))
    end

    if length(Bp) > 0 # we don't want to add zero vectors to lambda...
      push!(lambda, lambda_p)
    end

    M_sat = saturation(zero_Mi, monoid_algebra_ideal(kQ,p.prime))
    if !is_zero(p.prime) && !is_zero(M_sat)
      N, _ = quo(Mi, M_sat)
    end
    if is_zero(N) #TODO: should this be zero a some point? 
      break
    end
  end
  # return summands, hcat(lambda...)
  return IrrSum(kQ,summands), hcat(lambda...)
end

@doc raw"""
    irreducible_decomposition(I::MonoidAlgebraIdeal)

Return an irreducible decomposition of $I$.

# Examples
```jldoctest
julia> kQ = monoid_algebra([[1,0],[0,1]],QQ)
monoid algebra over rational field with cone of dimension 2

julia> x,y = gens(kQ)
2-element Vector{MonoidAlgebraElem{QQFieldElem, MonoidAlgebra{QQFieldElem, MPolyDecRing{QQFieldElem, QQMPolyRing}}}}:
 x_1
 x_2

julia> I = ideal(kQ,[x^4,x^2*y^2,y^4])
ideal over monoid algebra over rational field with cone of dimension 2 generated by x_1^4, x_1^2*x_2^2, x_2^4

julia> W = irreducible_decomposition(I)
2-element Vector{MonoidAlgebraIdeal{MonoidAlgebraElem{QQFieldElem, MonoidAlgebra{QQFieldElem, MPolyDecRing{QQFieldElem, QQMPolyRing}}}}}:
 ideal over monoid algebra over rational field with cone of dimension 2 generated by x_1^2, x_1^2*x_2, x_2^4, x_1*x_2^4
 ideal over monoid algebra over rational field with cone of dimension 2 generated by x_1^4, x_1^4*x_2, x_2^2, x_1*x_2^2

julia> I == intersect(W)
true
```
"""
function irreducible_decomposition(I::MonoidAlgebraIdeal)
  kQ = base_ring(I)

  J, _ = irreducible_hull(quotient_ring_as_module(I))
  if is_normal(kQ)
    return [_get_irreducible_ideal(kQ, I) for I in J.indec_injectives]
  else
    return [_get_irreducible_ideal_unsaturated(kQ, I) for I in J.indec_injectives]
  end
end

@doc raw"""
    _get_irreducible_ideal(kQ::MonoidAlgebra, J::IndecInj)

Given a monoid algebra $k[Q]$ and an indecomposable injective $J = k\{a + F - Q\}$ return the irreducible ideal $W\subseteq k[Q]$ such that $J_Q = k[Q]/W$ ($Q$-graded part of $J$). 

!!! note
    The monoid algebra $k[Q]$ must be normal. 
"""
function _get_irreducible_ideal(kQ::MonoidAlgebra, J::IndecInj)
  B_i = Vector{Vector{Vector{Int}}}()

  for h in hyperplanes(kQ)
    if is_subset(J.face.poly, h.hyperplane)
      B_h = generators_W_H(kQ, h, J.vector)
      push!(B_i, B_h)
    end
  end

  G_W = Vector{MPolyDecRingElem}()
  for b in B_i
    for bb in b
      a_v = Vector{ZZRingElem}()
      for a in bb
        push!(a_v, a)
      end
      push!(G_W, monomial_basis(kQ.algebra, a_v)[1])
    end
  end
  return ideal(kQ, G_W)
end

#compute the irreducible ideal (kQ unsaturated) (Algorithm 3.15 in HM05) 
function _get_irreducible_ideal_unsaturated(kQ::MonoidAlgebra, J::IndecInj)
  @assert base_ring(J.face.prime) == kQ.algebra
  @assert is_pointed(kQ) "k[Q] must be pointed"

  #get polyhedron a + ZF
  F = J.face.poly
  a = J.vector

  #get saturation of semigroup
  kQsat = saturation(kQ)

  #compute irreducible ideal in kQsat
  V = _get_irreducible_ideal(kQsat, J)

  # kQ as a kQsat-module
  I_kQ = saturation_ideal(kQ)

  W = intersect(I_kQ,V) #this is an ideal in kQsat

  B = [degree(Vector{Int},w) for w in filter(!is_zero,gens(W))]
  
  #intersection of p_D for all facets D of F
  I_D = Vector{MPolyQuoIdeal}()
  for d in facets(F)
    for p in faces(kQ)
      _facet = polyhedron([d],affine_hull(F))
      if _facet != F && p.poly == _facet
        push!(I_D, p.prime)
      end
    end
  end
  if dim(F) == 1 #the 0-dim face is always a facet of a 1-dim face
    push!(I_D,faces(kQ)[1].prime)
  end
  if length(I_D) > 0
    I = intersect(I_D...)
  else
    I = ideal(kQ.algebra,[])
  end

  #get W as an ideal in kQ
  _B = []
  for b in B # check if generators of W lie in kQ!!!
    if is_in_semigroup(kQ,b)
      push!(_B,b)
    end
  end
  W = ideal(kQ.algebra, [monomial_basis(kQ, b)[1] for b in _B])
  W_bar = quotient_ring_as_module(W)

  while !is_zero(W_bar)
    #get generators mod ZF
    W_F = mod_quotient(W_bar, J.face.prime)[1]
    if is_zero(W_F)
      break
    end
    i = 0
    for g in filter(!is_zero,gens(W_F))
      #check if deg(g) in a + ZF
      d_vec = degree(Vector{Int},g)
      if is_in_aZF(a,J.face,d_vec)
        continue
      end
      push!(_B,d_vec)
      i = i + 1
    end
    if i == 0
      break
    end

    #update W_bar
    W_bar = quotient_ring_as_module(ideal(kQ.algebra,[monomial_basis(kQ,d)[1] for d in _B]))

    if !is_zero(I)
       sat_W_bar = mod_saturate(W_bar,I)
      append!(_B,filter(!is_zero,[degree(g) for g in gens(sat_W_bar)]))

      # _W = ideal(kQ.algebra,[monomial_basis(kQ,d)[1] for d in _B])
      _G = filter(!is_zero,gens(sat_W_bar))

      #update W_bar
      W_bar,_ = quo(W_bar,_G)
    end
  end

  return ideal(kQ,[monomial_basis(kQ,b)[1] for b in _B])
end

# ── Combinatorial prototype for _get_irreducible_ideal_unsaturated (optimization 3) ──
#
# The key identity: b ∈ (W : p_F) \ W  iff
#   (a) b ∈ Q, b ∉ Q-span(_B)
#   (b) for every generator e of p_F:  b + e ∈ Q-span(_B)
# Candidates are {w − e : w ∈ _B, e ∈ prime_gens} ∩ Q, then filtered by (a) and (b).
# This replaces the mod_quotient Gröbner-basis call with semigroup-membership queries
# (which are already cached).  mod_saturate is kept as a Singular call for now.
#
# Toggle on/off with:  use_comb_unsaturated!(true)   /   use_comb_unsaturated!(false)

# ∃ w ∈ B, b − w ∈ Q  (i.e. b is in the ideal generated by B in Q-order)
function _in_Q_ideal(B_set::Set{Vector{Int}}, b::Vector{Int}, kQ::MonoidAlgebra)
    return any(is_in_semigroup(kQ, b .- w) for w in B_set)
end

# Returns degree vectors of new elements of (W : prime) \ W that are not in a + ZF.
function _comb_new_generators(B_set::Set{Vector{Int}}, prime_gens::Vector{Vector{Int}},
                               kQ::MonoidAlgebra, a::Vector{Int}, face::FaceQ)
    seen   = Set{Vector{Int}}()
    result = Vector{Int}[]
    for w in collect(B_set), e in prime_gens
        c = w .- e
        c ∈ seen && continue
        push!(seen, c)
        is_in_semigroup(kQ, c) || continue
        _in_Q_ideal(B_set, c, kQ) && continue                                    # already in W
        all(_in_Q_ideal(B_set, c .+ e2, kQ) for e2 in prime_gens) || continue   # (W:p_F) cond
        is_in_aZF(a, face, c) && continue                                        # in a + ZF
        push!(result, c)
    end
    return result
end

function _get_irreducible_ideal_unsaturated_comb(kQ::MonoidAlgebra, J::IndecInj)
    @assert base_ring(J.face.prime) == kQ.algebra
    @assert is_pointed(kQ) "k[Q] must be pointed"

    F  = J.face.poly
    a  = J.vector

    kQsat = saturation(kQ)
    V     = _get_irreducible_ideal(kQsat, J)
    I_kQ  = saturation_ideal(kQ)
    W     = intersect(I_kQ, V)

    B = [degree(Vector{Int}, w) for w in filter(!is_zero, gens(W))]

    I_D = Vector{MPolyQuoIdeal}()
    for d in facets(F)
        for p in faces(kQ)
            _facet = polyhedron([d], affine_hull(F))
            if _facet != F && p.poly == _facet
                push!(I_D, p.prime)
            end
        end
    end
    if dim(F) == 1
        push!(I_D, faces(kQ)[1].prime)
    end
    has_sat = !isempty(I_D)
    I = has_sat ? intersect(I_D...) : ideal(kQ.algebra, [])

    _B     = filter(b -> is_in_semigroup(kQ, b), B)
    _B_set = Set{Vector{Int}}(_B)

    prime_gens = [degree(Vector{Int}, g) for g in filter(!is_zero, gens(J.face.prime))]
    zero_vec   = zeros(Int, length(a))

    while zero_vec ∉ _B_set
        new_gens = _comb_new_generators(_B_set, prime_gens, kQ, a, J.face)
        isempty(new_gens) && break
        for c in new_gens
            push!(_B, c); push!(_B_set, c)
        end

        if has_sat
            W_bar     = quotient_ring_as_module(ideal(kQ.algebra, [monomial_basis(kQ, b)[1] for b in _B]))
            sat_W_bar = mod_saturate(W_bar, I)
            for g in filter(!is_zero, gens(sat_W_bar))
                d_vec = degree(Vector{Int}, g)
                _in_Q_ideal(_B_set, d_vec, kQ) && continue
                push!(_B, d_vec); push!(_B_set, d_vec)
            end
        end
    end

    return ideal(kQ, [monomial_basis(kQ, b)[1] for b in _B])
end

# workaround for homomorphism between monoid algebras
function hom(kQ1::MonoidAlgebra, kQ2::MonoidAlgebra, V::Vector)
  return hom(kQ1.algebra,kQ2.algebra,V)
end

function monomial_basis(kQ::MonoidAlgebra, a::Vector{Int})
  return monomial_basis(kQ.algebra, a)
end

@doc raw"""
    irreducible_resolution(M::SubquoModule{<:MonoidAlgebraElem}, i::Union{Int,Nothing}=nothing)

Return an irreducible resolution of $M$. If $i$ is specified then the resolution
is only computed up to cohomological degree $i$.

# Examples
```jldoctest
julia> kQ = monoid_algebra([[1,0],[0,1]],QQ)
monoid algebra over rational field with cone of dimension 2

julia> x,y = gens(kQ)
2-element Vector{MonoidAlgebraElem{QQFieldElem, MonoidAlgebra{QQFieldElem, MPolyDecRing{QQFieldElem, QQMPolyRing}}}}:
 x_1
 x_2

julia> I = ideal(kQ,[x^4,x^2*y^2,y^4])
ideal over monoid algebra over rational field with cone of dimension 2 generated by x_1^4, x_1^2*x_2^2, x_2^4

julia> M = quotient_ring_as_module(I)
Graded subquotient of graded submodule of kQ^1 with 1 generator
  1: 1*e[1]
by graded submodule of kQ^1 with 3 generators
  1: x_1^4*e[1]
  2: x_1^2*x_2^2*e[1]
  3: x_2^4*e[1]

julia> irr_res = irreducible_resolution(M)
irreducible resolution 
  W^0 -> W^1
where 
 W^0 = direct sum of
    k{[1, 3] + F - Q}_Q, where p_F = Ideal (x_2, x_1, x_1*x_2)
    k{[3, 1] + F - Q}_Q, where p_F = Ideal (x_2, x_1, x_1*x_2)
 W^1 = direct sum of
    k{[1, 1] + F - Q}_Q, where p_F = Ideal (x_2, x_1, x_1*x_2)
of Graded subquotient of graded submodule of kQ^1 with 1 generator
  1: 1*e[1]
by graded submodule of kQ^1 with 3 generators
  1: x_1^4*e[1]
  2: x_1^2*x_2^2*e[1]
  3: x_2^4*e[1]
over monoid algebra over rational field with cone of dimension 2
```
"""
function irreducible_resolution(M::SubquoModule{<:MonoidAlgebraElem}, i::Union{Int,Nothing}=nothing)
  kQ = base_ring(M)
  @assert generates_Zd(kQ) "The semigroup should generate ZZ^d."
  # @assert is_Q_graded(M) "M should be Q-graded."

  # if !is_normal(kQ)
  #   R_Q = saturation(kQ).algebra
  # else
    R_Q = kQ.algebra
  # end
  Mi = M # current module in resolution

  #initialize
  gi = identity_map(Mi)
  irreducible_sums = Vector{IrrSum}()
  cochain_maps = Vector{SubQuoHom}()

  j = 1
  while !is_zero(Mi) #until cokernel Mi is zero
    #compute irreducible hull
    Ji, _lambda = irreducible_hull(Mi, j)
    
    #get Q-graded part
    Wi = kQ_module(Ji)

    #multiply rows of lambda by degrees of generators of Mi
    m, n = size(_lambda)
    lambda = zero(_lambda)
    for ii in 1:m
      x_ii = monomial_basis(R_Q, degree(Mi[ii]))[1]
      for jj in 1:n
        lambda[ii, jj] = x_ii * _lambda[ii, jj]
        # lambda[ii, jj] = one(R_Q) * _lambda[ii, jj]
      end
    end

    #define injective map Mi -> Wi
    fi = hom(Mi, Wi, matrix(lambda))
    @assert is_injective(fi) "fi not injective"
    @assert is_welldefined(fi) "fi not well-defined"

    #get boundary map W{i-1} -> Wi
    hi = gi*fi

    #compute cokernel, filtering zero image generators to avoid trivial relations
    nz_img_gens = filter(!is_zero, [hi(g) for g in gens(domain(hi))])
    if isempty(nz_img_gens)
      Mi, gi = Wi, identity_map(Wi)
    else
      Mi, gi = quo(Wi, sub(Wi, nz_img_gens)[1])
    end

    push!(irreducible_sums, Ji)
    push!(cochain_maps, hi)

    # end at cohomological degree i
    if !isnothing(i) && j == i + 1
      break
    end
    j = j + 1
  end

  #get cochain complex
  C = cochain_complex(cochain_maps)

  return IrrRes(M, irreducible_sums, cochain_maps, C)
end

@doc raw"""
    injective_resolution(M::SubquoModule{<:MonoidAlgebraElem}, i::Int; shift::Symbol=:bound)

Return an injective resolution of $M$ up to cohomological degree i.

The keyword `shift` selects how the initial shift is computed:
* `:bound` (default) uses [`compute_shift_bound`](@ref).
* `:helm_miller` uses [`compute_shift`](@ref) (the shift described in
  Helm-Miller)

# Examples
```jldoctest
julia> kQ = monoid_algebra([[1,0],[0,1]],QQ)
monoid algebra over rational field with cone of dimension 2

julia> x,y = gens(kQ)
2-element Vector{MonoidAlgebraElem{QQFieldElem, MonoidAlgebra{QQFieldElem, MPolyDecRing{QQFieldElem, QQMPolyRing}}}}:
 x_1
 x_2

julia> I = ideal(kQ,[x^4,x^2*y^2,y^4])
ideal over monoid algebra over rational field with cone of dimension 2 generated by x_1^4, x_1^2*x_2^2, x_2^4

julia> M = quotient_ring_as_module(I)
Graded subquotient of graded submodule of kQ^1 with 1 generator
  1: 1*e[1]
by graded submodule of kQ^1 with 3 generators
  1: x_1^4*e[1]
  2: x_1^2*x_2^2*e[1]
  3: x_2^4*e[1]

julia> injective_resolution(M,2)
injective resolution 
  J^0 -> J^1 -> J^2
where 
 J^0 = direct sum of
    k{[1, 3] + F - Q}, where p_F = Ideal (x_2, x_1, x_1*x_2)
    k{[3, 1] + F - Q}, where p_F = Ideal (x_2, x_1, x_1*x_2)
 J^1 = direct sum of
    k{[-1, 3] + F - Q}, where p_F = Ideal (x_2, x_1, x_1*x_2)
    k{[1, 1] + F - Q}, where p_F = Ideal (x_2, x_1, x_1*x_2)
    k{[3, -1] + F - Q}, where p_F = Ideal (x_2, x_1, x_1*x_2)
 J^2 = direct sum of
    k{[-1, -1] + F - Q}, where p_F = Ideal (x_2, x_1, x_1*x_2)
of Graded subquotient of graded submodule of kQ^1 with 1 generator
  1: 1*e[1]
by graded submodule of kQ^1 with 3 generators
  1: x_1^4*e[1]
  2: x_1^2*x_2^2*e[1]
  3: x_2^4*e[1]
over monoid algebra over rational field with cone of dimension 2
```

```jldoctest
julia> kQ = monoid_algebra([[0,1],[1,1],[2,1]],QQ)
monoid algebra over rational field with cone of dimension 2

julia> x,y,z = gens(kQ)
3-element Vector{MonoidAlgebraElem{QQFieldElem, MonoidAlgebra{QQFieldElem, MPolyQuoRing{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}}}:
 x_1
 x_2
 x_3

julia> F = graded_free_module(kQ,2)
Graded free module monoid algebra over rational field with cone of dimension 2^2([0 0]) of rank 2 over monoid algebra over rational field with cone of dimension 2

julia> a = kQ[y y;0 x^2]
[x_2     x_2]
[  0   x_1^2]

julia> b = kQ[x^2*z 0; x^4*y 0; 0 x^5*y; 0 z^3]
[x_1^2*x_3           0]
[x_1^4*x_2           0]
[        0   x_1^5*x_2]
[        0       x_3^3]

julia> M = SubquoModule(F,a,b)
Graded subquotient of graded submodule of F with 2 generators
  1: x_2*e[1] + x_2*e[2]
  2: x_1^2*e[2]
by graded submodule of F with 4 generators
  1: x_1^2*x_3*e[1]
  2: x_1^4*x_2*e[1]
  3: x_1^5*x_2*e[2]
  4: x_3^3*e[2]

julia> injective_resolution(M,2)
injective resolution
  J^0 -> J^1 -> J^2
where
 J^0 = direct sum of
    k{[1, 4] + F - Q}, where p_F = Ideal (x_1, x_2, x_3, x_1*x_3)
    k{[5, 7] + F - Q}, where p_F = Ideal (x_1, x_2, x_3, x_1*x_3)
    k{[4, 7] + F - Q}, where p_F = Ideal (x_1, x_2, x_3, x_1*x_3)
    k{[0, 2] + F - Q}, where p_F = Ideal (x_2, x_3, x_1*x_3)
    k{[1, 2] + F - Q}, where p_F = Ideal (x_1, x_2, x_1*x_3)
 J^1 = direct sum of
    k{[1, 2] + F - Q}, where p_F = Ideal (x_1, x_2, x_3, x_1*x_3)
    k{[0, 3] + F - Q}, where p_F = Ideal (x_1, x_2, x_3, x_1*x_3)
    k{[-1, 3] + F - Q}, where p_F = Ideal (x_1, x_2, x_3, x_1*x_3)
    k{[5, 4] + F - Q}, where p_F = Ideal (x_1, x_2, x_3, x_1*x_3)
    k{[0, 5] + F - Q}, where p_F = Ideal (x_1, x_2, x_3, x_1*x_3)
    k{[4, 6] + F - Q}, where p_F = Ideal (x_1, x_2, x_3, x_1*x_3)
    k{[3, 6] + F - Q}, where p_F = Ideal (x_1, x_2, x_3, x_1*x_3)
    k{[-1, -2] + F - Q}, where p_F = Ideal (x_2, x_3, x_1*x_3)
    k{[-4, -2] + F - Q}, where p_F = Ideal (x_1, x_2, x_1*x_3)
 J^2 = direct sum of
    k{[0, 0] + F - Q}, where p_F = Ideal (x_1, x_2, x_3, x_1*x_3)
    k{[-1, 1] + F - Q}, where p_F = Ideal (x_1, x_2, x_3, x_1*x_3)
    k{[-1, 2] + F - Q}, where p_F = Ideal (x_1, x_2, x_3, x_1*x_3)
    k{[-2, 2] + F - Q}, where p_F = Ideal (x_1, x_2, x_3, x_1*x_3)
    k{[3, 5] + F - Q}, where p_F = Ideal (x_1, x_2, x_3, x_1*x_3)
    k{[2, 5] + F - Q}, where p_F = Ideal (x_1, x_2, x_3, x_1*x_3)
of Graded subquotient of graded submodule of F with 2 generators
  1: x_2*e[1] + x_2*e[2]
  2: x_1^2*e[2]
by graded submodule of F with 4 generators
  1: x_1^2*x_3*e[1]
  2: x_1^4*x_2*e[1]
  3: x_1^5*x_2*e[2]
  4: x_3^3*e[2]
over monoid algebra over rational field with cone of dimension 2
```
"""
function injective_resolution(M::SubquoModule{<:MonoidAlgebraElem}, i::Int; shift::Symbol=:bound)
  kQ = base_ring(M)
  @assert generates_Zd(kQ) "The semigroup should generate ZZ^d."

  G = grading_group(kQ)

  #compute irreducible resolution of shifted module
  if shift === :bound
    a_shift = compute_shift_bound(M, i+1)
  elseif shift === :helm_miller
    a_shift = compute_shift(M, i+1)
  elseif shift === :milp
    a_shift = compute_shift_milp(M, i+1)
  elseif shift === :milp_bound
    a_shift = compute_shift_milp_bound(M, i+1)
  else
    error("unknown shift strategy :$shift; expected :bound, :helm_miller, :milp, or :milp_bound")
  end
  
  M_a = twist(M, -G(a_shift))
  irr_res = irreducible_resolution(M_a,i)

  #get injective modules up to cohomological degree i, i.e. J^0, J^1, ...,J^i
  inj_modules = Vector{InjMod}()
  for j in 1:min(i + 1, length(irr_res.irr_sums))
    shifted_comp = map(
      indec -> IndecInj(indec.face, indec.vector - a_shift), irr_res.irr_sums[j].indec_injectives
    )
    push!(inj_modules, InjMod(kQ, shifted_comp))
  end

  #get all needed maps (as k-matrix or k[Q]-matrix?)
  cochain_maps = [
    matrix(irr_res.cochain_maps[k]) for k in eachindex(irr_res.cochain_maps) if 1 <= k <= i+1
  ]
  return InjRes(M, inj_modules, cochain_maps, length(inj_modules)-1, irr_res, a_shift)
end

# Generic rank of a SubquoModule M, i.e. dim_{Frac R} (M ⊗ Frac R). For f.g. M
# over a domain (which `k[Q]` is), M is torsion iff `annihilator(M) ≠ 0`, in
# which case the rank is 0. For non-torsion M we fall back to
# `rank` of the presentation matrix; that may hit "not implemented" over
# `MPolyQuoRing`-backed monoid algebras, in which case it errors.
function _generic_rank(M::SubquoModule)
  is_zero(M) && return 0
  !is_zero(annihilator(M)) && return 0
  pres = presentation(M)
  F0 = pres[0]
  A  = matrix(map(pres, 1))
  return Int(ngens(F0) - rank(A))
end

# Reduce `b ∈ ZZ^d` modulo the sublattice `ZF` spanned by the columns of `A`
# (the face generators). The indecomposable injective `E(F, a)` depends on `a`
# only modulo `ZF`, so degrees produced by `ZF_basis` are well-defined only up
# to that equivalence; this picks a canonical representative via HNF reduction.
function _reduce_mod_ZF(b::Vector{Int}, A::Union{Matrix{Int},Nothing})
  (A === nothing || size(A, 2) == 0) && return copy(b)
  H = hnf(transpose(matrix(ZZ, A)))           # rows = ZF basis in HNF
  bb = ZZRingElem[ZZ(b[k]) for k in 1:length(b)]
  for i in 1:nrows(H)
    j = 0
    for k in 1:ncols(H)
      if !iszero(H[i, k]); j = k; break; end
    end
    j == 0 && continue
    q = fdiv(bb[j], H[i, j])                   # floor div → canonical rep
    for k in 1:ncols(H)
      bb[k] -= q * H[i, k]
    end
  end
  return [Int(x) for x in bb]
end

@doc raw"""
    graded_bass_numbers(M::SubquoModule{<:MonoidAlgebraElem}, p::FaceQ, i::Int)

Return the graded Bass numbers of `M` at the prime `p = p_F` up to cohomological
degree `i`. The result is a `Vector` of length `i+1`; entry `j+1` is a dictionary
mapping each $\mathbb{Z}^d$-degree `a` (a canonical representative modulo
$\mathbb{Z}F$) to the multiplicity $\mu^j(p_F, M)_a$, i.e. the number of copies
of the indecomposable injective $E(F, a) = k\{a + F - Q\}$ in the `j`-th term of
the minimal injective resolution of `M`. Degrees are reduced modulo $\mathbb{Z}F$
because $E(F, a) \cong E(F, a+v)$ as graded $k[Q]$-modules for $v \in \mathbb{Z}F$
(via the $k[\mathbb{Z}F]$-action on $E_F = k\{F - Q\}$).
"""
function graded_bass_numbers(M::SubquoModule{<:MonoidAlgebraElem}, p::FaceQ, i::Int)
  kQ = base_ring(M)

  # The whole-cone face (p_F = (0)) is special. Over R_{p_F} = Frac(k[Q]) — a
  # field, since k[Q] is a domain — every module is projective, so
  #   μ^0((0), M) = rank(M),   μ^j((0), M) = 0 for j > 0.
  # The shift is also degree-independent (ZF = ZZ^d), hence stored under 0.
  # This short-circuit dodges the degenerate `mod_quotient(-, (0))` path that
  # currently crashes OSCAR's `hom(FreeMod, FreeMod)`.
  if is_zero(p.prime)
    d = rank(grading_group(kQ))
    r = _generic_rank(M)
    out = [Dict{Vector{Int},Int}() for _ in 0:i]
    r > 0 && (out[1][zeros(Int, d)] = r)
    return out
  end

  # k[Q]/p_F as a module, then Ext^j(k[Q]/p_F, M) for j = 0..i in one pass
  RpF    = quotient_ring_as_module(monoid_algebra_ideal(kQ, p.prime))
  fr     = free_resolution(RpF; length=i+2)
  lifted = hom(fr.C[first(Hecke.map_range(fr.C)):-1:1], M)

  out = Vector{Dict{Vector{Int},Int}}()
  for j in 0:i
    # for the whole-cone face p_F = (0), k[Q]/p_F is free and the homology can
    # come back as a FreeMod; ZF_basis expects a SubquoModule, so coerce.
    Hj    = homology(lifted, j)
    Ej    = simplify_light(Hj isa SubquoModule ? Hj : sub(Hj, gens(Hj))[1])[1]
    tally = Dict{Vector{Int},Int}()
    for b in ZF_basis(Ej, p)                 # k[ZF]-basis of (0 :_{Ej} p)[ZF]
      a = _reduce_mod_ZF(degree(Vector{Int}, b), p.A)
      tally[a] = get(tally, a, 0) + 1
    end
    push!(out, tally)
  end
  return out
end

@doc raw"""
    is_minimal(res::InjRes; verbose::Bool=false)

Check that the injective resolution `res` of `M = res.mod` is minimal, by
comparing, for every face `F` and cohomological degree `j`, the multiplicity of
each indecomposable injective $E(F, a)$ in `res.inj_mods[j+1]` against the
independently computed graded Bass number $\mu^j(p_F, M)_a$ (see
[`graded_bass_numbers`](@ref)).

Returns `true` iff all multiplicities agree. With `verbose=true`, logs each
mismatch (face, cohomological degree, claimed vs. expected multiplicities).
"""
function is_minimal(res::InjRes; verbose::Bool=false)
  M  = res.mod
  kQ = base_ring(M)
  ok = true

  for p in faces(kQ)
    bass = graded_bass_numbers(M, p, res.upto)
    for j in 0:res.upto
      claimed = Dict{Vector{Int},Int}()
      for ind in res.inj_mods[j+1].indec_injectives
        ind.face.poly == p.poly || continue
        a = _reduce_mod_ZF(ind.vector, p.A)
        claimed[a] = get(claimed, a, 0) + 1
      end
      if claimed != bass[j+1]
        ok = false
        verbose && @info "non-minimal term" face=p cohomological_degree=j claimed expected=bass[j+1]
      end
    end
  end
  return ok
end

@doc raw"""
    injective_resolution(I::MonoidAlgebraIdeal, i::Int; shift::Symbol=:bound)

Return an injective resolution of $M = k[Q]/I$ up to cohomological degree i.
See the module version for the `shift` keyword.
"""
function injective_resolution(I::MonoidAlgebraIdeal, i::Int; shift::Symbol=:bound)
  return injective_resolution(quotient_ring_as_module(I), i; shift=shift)
end

function injective_hull(M::SubquoModule{<:MonoidAlgebraElem})
  kQ = base_ring(M)
  @assert generates_Zd(kQ) "The semigroup should generate ZZ^d."
  # @assert is_Q_graded(M) "M should be Q-graded."

  R_Q = kQ.algebra
  G = grading_group(kQ)

  #compute irreducible hull of shifted module
  a_shift = compute_shift(M, 0)
  M_a = twist(M, -G(a_shift))
  J, _lambda = irreducible_hull(M_a, 0)

  #multiply rows of lambda by degrees of generators of M
  m, n = size(_lambda)
  lambda = zero(_lambda)
  for ii in 1:m
    for jj in 1:n
      lambda[ii, jj] = monomial_basis(R_Q, degree(M[ii]))[1] * _lambda[ii, jj]
    end
  end

  # undo shift to obtain injective hull of M
  E = InjMod(kQ,map(
    indec -> IndecInj(indec.face, indec.vector - a_shift), J.indec_injectives
  ))

  return E,lambda
end

# import local cohomology functions
include("LocalCohomology.jl")
include("ModuleFunctionality.jl")

end # module InjectiveResolutions

using .InjectiveResolutions

## Functions visible on the outside
export monoid_algebra
export faces
export hyperplanes
export saturation

export irreducible_resolution
export irreducible_decomposition

export injective_resolution
export old_injective_resolution

export local_cohomology
export local_cohomology_all
export zeroth_local_cohomology

export MonoidAlgebra
export MonoidAlgebraIdeal
export MonoidAlgebraElem
export AffineSemigroup
export affine_semigroup
export ambient_dimension
export semigroup_generators
export polyhedral_cone
export cone
export saturation_ideal
export saturation_map
export holes_module
export is_Q_graded

export compute_shift
export compute_shift_bound
export InjMod
export IndecInj
export monomial_matrix
export irreducible_hull
export kQ_module
export injective_hull
export Q_graded_part
export mod_quotient
export monoid_algebra_ideal
export FaceQ
export prime_of_face
export generators_W_H
export ZF_basis
export generates_Zd
export degrees_of_bass_numbers
export degrees_of_bass_numbers_bound
export graded_bass_numbers
export is_minimal
export in_intersection
export in_semigroup
export is_in_aZF
export coefficients_wrt_generators
export relevant_generators
export relevant_relations
export _get_irreducible_ideal
export underlying_element
export mod_saturate
export underlying_ideal
export _get_irreducible_ideal_unsaturated
export is_in_semigroup