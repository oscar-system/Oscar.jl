module InjectiveResolutions
using ..Oscar
using ..Oscar: IntegerUnion, ModuleGens, SubModuleOfFreeModule, OFPModule, Orderings
using ..Oscar: _extend_free_resolution, _graded_kernel, _presentation_minimal, _reduce
using ..Oscar: images_of_generators, oscar_free_module, oscar_generators
using ..Oscar: singular_freemodule
using ..Oscar: pretty, terse, is_terse, Lowercase, Indent, Dedent

# functions from OSCAR for which this module adds methods
import ..Oscar:
  _build_sparse_row,
  _degree_fast,
  _saturation,
  annihilator,
  augmentation_map,
  base_ring,
  canonical_unit,
  characteristic,
  coefficient_ring,
  coefficients,
  cochain_complex,
  cone,
  coordinates,
  coordinates_atomic,
  coordinates_via_transform,
  default_ordering,
  degree,
  dim,
  divexact,
  divides,
  domain,
  codomain,
  elem_type,
  evaluate,
  faces,
  free_resolution,
  gens,
  grading_group,
  hom,
  hyperplanes,
  ideal,
  in_atomic,
  intersect,
  inv,
  irreducible_decomposition,
  is_domain_type,
  is_exact,
  is_graded,
  is_homogeneous,
  is_nilpotent,
  is_normal,
  is_pointed,
  is_subset,
  is_unit,
  is_zero,
  is_zm_graded,
  kernel,
  kernel_atomic,
  krull_dim,
  lift_std,
  matrix,
  minimal_generating_set,
  monomial_basis,
  normal_form,
  number_of_generators,
  number_of_variables,
  one,
  parent,
  parent_type,
  primitive_generator,
  primitive_generator_with_scaling_factor,
  prune_with_map,
  radical,
  rand,
  rank,
  saturation,
  singular_generators,
  singular_module,
  singular_poly_ring,
  sparse_row,
  standard_basis,
  syzygy_module,
  twist,
  underlying_module,
  zero,
  zonotope

import ..Oscar.Singular:
  FreeModule,
  has_global_ordering,
  svector,
  Module

import Base:
  +,
  -,
  *,
  ==,
  deepcopy_internal,
  getindex

## Functions and types visible on the outside
export AffineSemigroup
export FaceQ
export IndecInj
export InjMod
export InjRes
export IrrRes
export IrrSum
export MonoidAlgebra
export MonoidAlgebraElem
export MonoidAlgebraIdeal
export MonomialMatrix
export affine_semigroup
export augmentation_map
export cochain_maps
export cohomological_degree
export degree_shift
export degrees_of_bass_numbers
export holes_module
export indecomposable_injectives
export injective_hull
export injective_modules
export injective_resolution
export irreducible_resolution
export irreducible_sums
export is_minimal
export local_cohomology
export local_cohomology_all
export monoid_algebra
export monoid_algebra_ideal
export monomial_matrix
export Q_graded_part
export sectors
export semigroup_generators
export zeroth_local_cohomology

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

# irreducible sum as a finitely generated k[Q]-module
function underlying_module(I::IrrSum)
  if I.kQ_module === nothing
    I.kQ_module = _compute_q_graded_part(I.monoid_algebra, I.indec_injectives)
  end
  return I.kQ_module
end

function Q_graded_part(I::InjMod)
  if I.Q_graded_part === nothing
    I.Q_graded_part = _compute_q_graded_part(I.monoid_algebra, I.indec_injectives)
  end
  return I.Q_graded_part
end

function _compute_q_graded_part(kQ::MonoidAlgebra, I::Vector{IndecInj})
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

@doc raw"""
    MonomialMatrix{T}

A map between direct sums of indecomposable injectives (`T = InjMod`) or of
their $Q$-graded parts (`T = IrrSum`) in the sense of [HM05](@cite). It consists of a scalar
matrix whose rows and columns are labelled by the summands of the source and
the target. The entry in row $r$ and column $c$ is the coefficient of the
monomial $x^{a_c - a_r}$ of the map between the corresponding summands.
"""
struct MonomialMatrix{T <: Union{InjMod, IrrSum}}
  matrix::MatElem   # scalar matrix over the coefficient field
  source::T
  target::T
  index::Int        # cohomological degree i of the differential d^i
end

matrix(mm::MonomialMatrix) = mm.matrix
domain(mm::MonomialMatrix) = mm.source
codomain(mm::MonomialMatrix) = mm.target
cohomological_degree(mm::MonomialMatrix) = mm.index

# scalar coefficients of a matrix with monomial entries over the monoid algebra
function _scalar_matrix(kQ::MonoidAlgebra, A::MatElem)
  k = coefficient_ring(kQ)
  return map(x -> evaluate(x, ones(elem_type(k), ngens(kQ))), A)
end

struct InjRes #ZZ^d-graded injective resolution
  mod::SubquoModule
  inj_mods::Vector{InjMod}
  augmentation_map::MatElem                     # scalar matrix of M -> J^0
  cochain_maps::Vector{MonomialMatrix{InjMod}}  # d^0, d^1, ..., d^{upto-1}
  upto::Int
  Q_graded_part::IrrRes
  shift::Vector{Int}
end

@doc raw"""
    indecomposable_injectives(J::InjMod)
    indecomposable_injectives(W::IrrSum)

Return the indecomposable injectives $k\{a_i + F_i - Q\}$ whose direct sum is `J`
(respectively whose $Q$-graded parts form the irreducible sum `W`).
"""
indecomposable_injectives(J::Union{InjMod,IrrSum}) = J.indec_injectives

@doc raw"""
    monoid_algebra(J::InjMod)
    monoid_algebra(W::IrrSum)

Return the monoid algebra over which `J` (respectively `W`) is defined.
"""
monoid_algebra(J::Union{InjMod,IrrSum}) = J.monoid_algebra

@doc raw"""
    injective_modules(res::InjRes)

Return the injective modules $J^0, J^1, \dots, J^i$ of the injective resolution `res`.
"""
injective_modules(res::InjRes) = res.inj_mods

@doc raw"""
    cochain_maps(res::InjRes)
    cochain_maps(res::IrrRes)

Return the differentials $d^0, d^1, \dots$ of the resolution `res`, that is,
the maps between consecutive terms. For an injective resolution these are
`MonomialMatrix` objects, for an irreducible resolution module homomorphisms.
The map from the resolved module into the first term is [`augmentation_map`](@ref).
"""
cochain_maps(res::InjRes) = res.cochain_maps

@doc raw"""
    augmentation_map(res::InjRes)
    augmentation_map(res::IrrRes)

Return the map $\epsilon$ from the resolved module $M$ into the first term of
the resolution `res`. For an injective resolution this is the scalar matrix
whose entry $(i, j)$ is the coefficient with which the $i$-th generator of $M$
maps into the $j$-th summand of $J^0$, for an irreducible resolution it is the
module homomorphism $M \to \overline{W}^0$.
"""
augmentation_map(res::InjRes) = res.augmentation_map
augmentation_map(res::IrrRes) = res.cochain_maps[1]

@doc raw"""
    Q_graded_part(J::InjMod)
    Q_graded_part(res::InjRes)

Return the $Q$-graded part of an injective module as a finitely generated
module, respectively the irreducible resolution of the shifted module
$M(-\alpha)$ from which the injective resolution `res` was computed, see
[`degree_shift`](@ref).
"""
Q_graded_part(res::InjRes) = res.Q_graded_part

@doc raw"""
    degree_shift(res::InjRes)

Return the degree $\alpha \in \mathbb{Z}^d$ such that all Bass numbers of
$M(-\alpha)$ up to the length of `res` lie in $Q$. The resolution `res` was
obtained by shifting an irreducible resolution of $M(-\alpha)$ back by $\alpha$.
"""
degree_shift(res::InjRes) = res.shift

@doc raw"""
    irreducible_sums(res::IrrRes)

Return the irreducible sums $\overline{W}^0, \overline{W}^1, \dots$ of the
irreducible resolution `res`.
"""
irreducible_sums(res::IrrRes) = res.irr_sums
cochain_maps(res::IrrRes) = res.cochain_maps[2:end]

@doc raw"""
    cochain_complex(res::IrrRes)

Return the irreducible resolution `res` as a cochain complex of modules.
"""
cochain_complex(res::IrrRes) = res.cochain_complex

@doc raw"""
    is_exact(res::IrrRes)

Check whether the cochain complex of the irreducible resolution `res` is exact.
"""
is_exact(res::IrrRes) = is_exact(res.cochain_complex)

function Base.show(io::IO, mm::MonomialMatrix{InjMod})
  print(io, "Monomial matrix for J^", mm.index, " -> J^", mm.index + 1)
end

function Base.show(io::IO, mm::MonomialMatrix{IrrSum})
  print(io, "Monomial matrix for W^", mm.index, " -> W^", mm.index + 1)
end

function _show_indec(io::IO, J::IndecInj, ::Type{InjMod})
  println(io, "  k{", J.vector, " + F - Q}, where p_F = ", J.face.prime)
end
function _show_indec(io::IO, J::IndecInj, ::Type{IrrSum})
  println(io, "  k{", J.vector, " + F - Q}_Q, where p_F = ", J.face.prime)
end

function Base.show(io::IO, ::MIME"text/plain", mm::MonomialMatrix{T}) where T
  src_label, tgt_label = T == InjMod ? ("J^$(mm.index)", "J^$(mm.index + 1)") : ("W^$(mm.index)", "W^$(mm.index + 1)")
  println(io, "Monomial matrix for $src_label -> $tgt_label")
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
    io, "Injective module given by direct sum of ", length(J.indec_injectives), " indecomposable injectives"
  )
end

function Base.show(io::IO, ::MIME"text/plain", J::InjMod)
  println(io, "Injective module given by direct sum of indecomposable injectives")
  for Ji in J.indec_injectives
    _show_indec(io, Ji, InjMod)
  end
  print(pretty(io), "over ", Lowercase(), J.monoid_algebra)
end

function Base.show(io::IO,W::IrrSum)
  print(io, "Irreducible sum given by direct sum of ", length(W.indec_injectives), " components")
end

function Base.show(io::IO, ::MIME"text/plain", W::IrrSum)
  println(io, "Irreducible sum given by direct sum of components")
  for Ji in W.indec_injectives
    _show_indec(io, Ji, IrrSum)
  end
  print(pretty(io), "over ", Lowercase(), W.monoid_algebra)
end

function Base.show(io::IO, ::MIME"text/plain", res::InjRes)
  println(io, "Injective resolution")
  println(io, "  ", join(["J^$i" for i in 0:res.upto], " -> "))
  println(io, "where")
  for i in eachindex(res.inj_mods)
    println(io, " J^$(i-1) = direct sum of")
    for Ji in res.inj_mods[i].indec_injectives
      print(io, "  "); _show_indec(io, Ji, InjMod)
    end
  end
  println(io, "of ", res.mod)
  print(pretty(io), "over ", Lowercase(), base_ring(res.mod))
end

function Base.show(io::IO, ::MIME"text/plain", res::IrrRes)
  println(io, "Irreducible resolution")
  println(io, "  ", join(["W^$i" for i in 0:(length(res.irr_sums) - 1)], " -> "))
  println(io, "where")
  for i in eachindex(res.irr_sums)
    println(io, " W^$(i-1) = direct sum of")
    for Ji in res.irr_sums[i].indec_injectives
      print(io, "  "); _show_indec(io, Ji, IrrSum)
    end
  end
  println(io, "of ", res.mod)
  print(pretty(io), "over ", Lowercase(), base_ring(res.mod))
end

@doc raw"""
    monomial_matrix(i::Int, res::IrrRes)
    monomial_matrix(i::Int, res::InjRes)

Return the differential $d^i$ of the (co)chain complex of `res` as a
`MonomialMatrix`.

Here `i` is the cohomological degree of the differential
$d^i \colon W^i \to W^{i+1}$ (resp. $J^i \to J^{i+1}$).
"""
function monomial_matrix(i::Int, res::IrrRes)
  n = length(res.irr_sums)
  @req 0 <= i <= n - 2 "cohomological degree i must be in 0:$(n - 2) for this resolution"
  A = _scalar_matrix(base_ring(res.mod), matrix(res.cochain_maps[i + 2]))
  return MonomialMatrix(A, res.irr_sums[i + 1], res.irr_sums[i + 2], i)
end

function monomial_matrix(i::Int, res::InjRes)
  n = length(res.cochain_maps)
  @req 0 <= i <= n - 1 "cohomological degree i must be in 0:$(n - 1) for this resolution"
  return res.cochain_maps[i + 1]
end

function Base.show(io::IO, ::MIME"text/plain", Ji::IndecInj)
  println(io, "Indecomposable injective")
  println(io, "  k{", Ji.vector, " + F - Q},")
  print(io, "where p_F = ", Ji.face.prime)
end

function Base.show(io::IO, Ji::IndecInj)
  print(
      io, "Indecomposable injective k{", Ji.vector, " + F - Q}, where p_F = ", Ji.face.prime
    )
end

@doc raw"""
    generators_W_H(kQ::MonoidAlgebra, H::HyperplaneQ, a::Vector{Int})

Return, for a hyperplane `H` bounding the cone $\mathbb{R}_{\geq 0}Q$ of the
monoid algebra `kQ` and a vector $a \in \mathbb{Z}^d$, a finite set $B$ such that

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
  rhs = H.A*a + H.b
  PaF = polyhedron(H.A, rhs) #a + RR H
  Z = zonotope(kQ)[1]
  for d in D
    I = intersect(PaF, d) #(a + RR H)\cap RR_+D
    if dim(I) >= 0
      for pt in lattice_points(I + Z)
        v = Vector{Int}(pt)
        # keep the points off the hyperplane a + RR H
        all(H.A*v .<= rhs) || push!(B, v)
      end
    end
  end
  return B
end

@doc raw"""
    degrees_of_bass_numbers(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)

Return the $\mathbb{Z}^d$-degrees of non-zero Bass numbers of `M` up to cohomological degree `i`.
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
    degrees_of_bass_numbers_bound(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)

Return a finite set $D$ of $\mathbb{Z}^d$-degrees such that every degree of a
non-zero Bass number of `M` up to cohomological degree `i` lies in $D + Q$.
Unlike `degrees_of_bass_numbers` this needs no $\mathrm{Ext}$ computation.
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
    # for normal Q the semigroup is the cone intersected with the lattice,
    # so membership is a check against the facet inequalities
    A_f, b_f = halfspace_matrix_pair(facets(cone(kQ)))
    A_m = Matrix(A_f)
    b_v = Vector(b_f)
    in_Q = b -> all(A_m*b .<= b_v)
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

Let `M` be finitely generated $\mathbb{Z}^d$-graded module over a monoid algebra $k[Q]$. This function computes $a\in \mathbb{Z}^d$
such that all $\mathbb{Z}^d$-degrees of non-zero Bass numbers of $M(-a)$ lie in $Q$.
"""
function compute_shift(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)
  #get all degrees of non-zero Bass numbers up to cohomological degree i
  return _shift_into_Q(base_ring(M), degrees_of_bass_numbers(M, i))
end

@doc raw"""
    compute_shift_bound(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)

Like `compute_shift`, but uses the cheap over-approximation
`degrees_of_bass_numbers_bound` instead of the exact Bass-number degrees.
"""
function compute_shift_bound(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)
  return _shift_into_Q(base_ring(M), degrees_of_bass_numbers_bound(M, i))
end

# MILP that minimises sum(a) subject to a + B[j] in Q for all j, where a in Q.
# Variables: [a (d); lambda_0 (n); lambda_1 (n); ...; lambda_k (n)], all integer.
function _shift_milp(kQ::MonoidAlgebra, degrees::Vector{Vector{Int}})
  isempty(degrees) && return zeros(Int, ambient_dimension(affine_semigroup(kQ)))
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
  # The multipliers are integers, so that a + b_j = A lambda_j lies in the
  # semigroup generated by the columns of A (rounding a real solution up
  # coordinatewise would not stay in the cone in general).
  int_vars = collect(d+1 : d+(k+1)*n)
  milp = mixed_integer_linear_program(P, c; integer_variables=int_vars, convention = :min)
  _, opt = solve_milp(milp)
  opt === nothing && error("MILP shift infeasible: no shift puts all Bass-number degrees into Q")
  return [Int(opt[i]) for i in 1:d]
end

@doc raw"""
    compute_shift_milp(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)

Like `compute_shift`, but chooses the shift $a$ with minimal coordinate sum
such that $a + b \in Q$ for every Bass-number degree $b$, by solving an
integer linear program.
"""
function compute_shift_milp(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)
  return _shift_milp(base_ring(M), degrees_of_bass_numbers(M, i))
end

@doc raw"""
    compute_shift_milp_bound(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)

Like `compute_shift_milp`, but uses
`degrees_of_bass_numbers_bound` for the Bass-number degrees.
"""
function compute_shift_milp_bound(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)
  return _shift_milp(base_ring(M), degrees_of_bass_numbers_bound(M, i))
end

@doc raw"""
    mod_quotient(M::SubquoModule, I::Ideal)

Return the submodule

$(0 :_M I) := \{m \in M \mid m\cdot I = 0\}.$

# Examples
```jldoctest
julia> R_Q, (x, y) = graded_polynomial_ring(QQ, [:x, :y]; weights = [[1, 0], [0, 1]])
(Graded multivariate polynomial ring in 2 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y])

julia> I = ideal(R_Q, [x^4, x^2*y^2, y^4])
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

julia> m = ideal(R_Q, [x, y])
Ideal generated by
  x
  y

julia> Oscar.InjectiveResolutions.mod_quotient(M, m)
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
    mod_saturate(M::SubquoModule, I::Ideal)

Compute the saturation

$(0 :_M I^\infty) := \{m \in M \mid m\cdot I^n = 0\text{ for some }n\in \NN_{>0}\}.$

# Examples
```jldoctest
julia> R_Q, (x, y) = graded_polynomial_ring(QQ, [:x, :y]; weights = [[1, 0], [0, 1]])
(Graded multivariate polynomial ring in 2 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y])

julia> I = ideal(R_Q, [x^4, x^2*y^2, y^4])
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

julia> m = ideal(R_Q, [x, y])
Ideal generated by
  x
  y

julia> Oscar.InjectiveResolutions.mod_saturate(M, m)
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
    return Bp, zero_matrix(kQ, 1, 1)
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
    @assert all(is_zero(_c_b[i]) || is_homogeneous(_c_b[i]) for i in 1:ngens(N)) "non-homogeneous coordinate, all-ones evaluation is wrong here"
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
    return Bp, zero_matrix(kQ, 1, 1)
  end
  lambda = _dual_basis_lambdas(k, ngens(N), p, degs, c_bs, C_rows)
  return Bp, matrix(kQ, map(kQ, hcat(lambda...)))
end

@doc raw"""
    irreducible_hull(M::SubquoModule{<:MonoidAlgebraElem}, j=0)

Return an irreducible hull of `M`, that is, an irreducible sum $\overline{W}$
together with a matrix $\Lambda$ defining an injective map $M \to \overline{W}$.
This is Algorithm 3.6 in [HM05](@cite) with the correction for socle
elements supported on several generators.
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
  return IrrSum(kQ,summands), hcat(lambda...)
end

@doc raw"""
    irreducible_decomposition(I::MonoidAlgebraIdeal)

Return an irreducible decomposition of `I`.

# Examples
```jldoctest
julia> kQ = monoid_algebra([[1, 0], [0, 1]], QQ)
Monoid algebra over rational field with cone of dimension 2

julia> x, y = gens(kQ)
2-element Vector{MonoidAlgebraElem{QQFieldElem, MonoidAlgebra{QQFieldElem, MPolyDecRing{QQFieldElem, QQMPolyRing}}}}:
 x_1
 x_2

julia> I = ideal(kQ, [x^4, x^2*y^2, y^4])
Ideal over monoid algebra over rational field with cone of dimension 2
generated by x_1^4, x_1^2*x_2^2, x_2^4

julia> W = irreducible_decomposition(I)
2-element Vector{MonoidAlgebraIdeal{MonoidAlgebraElem{QQFieldElem, MonoidAlgebra{QQFieldElem, MPolyDecRing{QQFieldElem, QQMPolyRing}}}}}:
 Ideal (x_1^2, x_1^2*x_2, x_2^4, x_1*x_2^4)
 Ideal (x_1^4, x_1^4*x_2, x_2^2, x_1*x_2^2)

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

Return the irreducible ideal $W \subseteq k[Q]$ with $J_Q = k[Q]/W$, the $Q$-graded
part of the indecomposable injective $J = k\{a + F - Q\}$ over the monoid algebra `kQ`.

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

      _G = filter(!is_zero,gens(sat_W_bar))

      #update W_bar
      W_bar,_ = quo(W_bar,_G)
    end
  end

  return ideal(kQ,[monomial_basis(kQ,b)[1] for b in _B])
end

# workaround for homomorphism between monoid algebras
function hom(kQ1::MonoidAlgebra, kQ2::MonoidAlgebra, V::Vector)
  return hom(kQ1.algebra,kQ2.algebra,V)
end

function monomial_basis(kQ::MonoidAlgebra, a::Vector{Int})
  return monomial_basis(kQ.algebra, a)
end

@doc raw"""
    irreducible_resolution(M::SubquoModule{<:MonoidAlgebraElem}, i::Union{Int,Nothing}=nothing; check::Bool=true)

Return an irreducible resolution of `M`. If `i` is specified then the resolution
is only computed up to cohomological degree `i`. With `check = false` the
internal verification that each map into an irreducible hull is well-defined
and injective is skipped.

# Examples
```jldoctest
julia> kQ = monoid_algebra([[1, 0], [0, 1]], QQ)
Monoid algebra over rational field with cone of dimension 2

julia> x, y = gens(kQ)
2-element Vector{MonoidAlgebraElem{QQFieldElem, MonoidAlgebra{QQFieldElem, MPolyDecRing{QQFieldElem, QQMPolyRing}}}}:
 x_1
 x_2

julia> I = ideal(kQ, [x^4, x^2*y^2, y^4])
Ideal over monoid algebra over rational field with cone of dimension 2
generated by x_1^4, x_1^2*x_2^2, x_2^4

julia> M = quotient_ring_as_module(I)
Graded subquotient of graded submodule of kQ^1 with 1 generator
  1: 1*e[1]
by graded submodule of kQ^1 with 3 generators
  1: x_1^4*e[1]
  2: x_1^2*x_2^2*e[1]
  3: x_2^4*e[1]

julia> irr_res = irreducible_resolution(M)
Irreducible resolution
  W^0 -> W^1
where
 W^0 = direct sum of
    k{[1, 3] + F - Q}_Q, where p_F = Ideal (x_1, x_2)
    k{[3, 1] + F - Q}_Q, where p_F = Ideal (x_1, x_2)
 W^1 = direct sum of
    k{[1, 1] + F - Q}_Q, where p_F = Ideal (x_1, x_2)
of Graded subquotient of graded submodule of kQ^1 with 1 generator
  1: 1*e[1]
by graded submodule of kQ^1 with 3 generators
  1: x_1^4*e[1]
  2: x_1^2*x_2^2*e[1]
  3: x_2^4*e[1]
over monoid algebra over rational field with cone of dimension 2
```
"""
function irreducible_resolution(M::SubquoModule{<:MonoidAlgebraElem}, i::Union{Int,Nothing}=nothing; check::Bool=true)
  kQ = base_ring(M)
  @req _generates_lattice(kQ) "the semigroup must generate ZZ^d"

  R_Q = kQ.algebra
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
    Wi = underlying_module(Ji)

    #multiply rows of lambda by degrees of generators of Mi
    m, n = size(_lambda)
    lambda = zero(_lambda)
    for ii in 1:m
      x_ii = monomial_basis(R_Q, degree(Mi[ii]))[1]
      for jj in 1:n
        lambda[ii, jj] = x_ii * _lambda[ii, jj]
      end
    end

    #define injective map Mi -> Wi
    fi = hom(Mi, Wi, matrix(lambda))
    @check is_welldefined(fi) "map into the irreducible hull is not well-defined"
    @check is_injective(fi) "map into the irreducible hull is not injective"

    #get boundary map W{i-1} -> Wi
    hi = gi*fi

    #compute cokernel, filtering zero image generators to avoid trivial relations
    nz_img_gens = filter(!is_zero, [hi(g) for g in gens(domain(hi))])
    if isempty(nz_img_gens)
      Mi, gi = Wi, identity_map(Wi)
    else
      Mi, gi = quo(Wi, sub(Wi, nz_img_gens)[1])
    end

    # replace the cokernel by a minimal presentation
    if !is_zero(Mi)
      Mi_min, phi = prune_with_map(Mi)
      gi = gi * inv(phi; check=false)
      Mi = Mi_min
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
    injective_resolution(M::SubquoModule{<:MonoidAlgebraElem}, i::Int; shift::Symbol=:bound, check::Bool=true)
    injective_resolution(I::MonoidAlgebraIdeal, i::Int; shift::Symbol=:bound, check::Bool=true)

Return an injective resolution of `M`, respectively of $k[Q]/I$, up to cohomological degree `i`.
With `check = false` the internal verification of the maps is skipped.

The module is first shifted so that the degrees of its Bass numbers at the maximal
ideal up to cohomological degree `i + d`, $d = \dim Q$, lie in $Q$, see Lemma 4.5 in [HM05](@cite).
The keyword `shift` selects how this shift is computed:
* `:bound` (default) uses `compute_shift_bound`, a cheap bound on the degrees of the Bass numbers.
* `:helm_miller` uses `compute_shift`, the shift described in [HM05](@cite), which needs the exact Bass numbers.
* `:milp_bound` and `:milp` choose the shift with minimal coordinate sum by an integer linear program, from the bound and from the exact Bass numbers respectively.

# Examples
```jldoctest
julia> kQ = monoid_algebra([[1, 0], [0, 1]], QQ)
Monoid algebra over rational field with cone of dimension 2

julia> x, y = gens(kQ)
2-element Vector{MonoidAlgebraElem{QQFieldElem, MonoidAlgebra{QQFieldElem, MPolyDecRing{QQFieldElem, QQMPolyRing}}}}:
 x_1
 x_2

julia> I = ideal(kQ, [x^4, x^2*y^2, y^4])
Ideal over monoid algebra over rational field with cone of dimension 2
generated by x_1^4, x_1^2*x_2^2, x_2^4

julia> M = quotient_ring_as_module(I)
Graded subquotient of graded submodule of kQ^1 with 1 generator
  1: 1*e[1]
by graded submodule of kQ^1 with 3 generators
  1: x_1^4*e[1]
  2: x_1^2*x_2^2*e[1]
  3: x_2^4*e[1]

julia> injective_resolution(M, 2)
Injective resolution
  J^0 -> J^1 -> J^2
where
 J^0 = direct sum of
    k{[1, 3] + F - Q}, where p_F = Ideal (x_1, x_2)
    k{[3, 1] + F - Q}, where p_F = Ideal (x_1, x_2)
 J^1 = direct sum of
    k{[-1, 3] + F - Q}, where p_F = Ideal (x_1, x_2)
    k{[1, 1] + F - Q}, where p_F = Ideal (x_1, x_2)
    k{[3, -1] + F - Q}, where p_F = Ideal (x_1, x_2)
 J^2 = direct sum of
    k{[-1, -1] + F - Q}, where p_F = Ideal (x_1, x_2)
of Graded subquotient of graded submodule of kQ^1 with 1 generator
  1: 1*e[1]
by graded submodule of kQ^1 with 3 generators
  1: x_1^4*e[1]
  2: x_1^2*x_2^2*e[1]
  3: x_2^4*e[1]
over monoid algebra over rational field with cone of dimension 2
```
"""
function injective_resolution(M::SubquoModule{<:MonoidAlgebraElem}, i::Int; shift::Symbol=:bound, check::Bool=true)
  kQ = base_ring(M)
  @req _generates_lattice(kQ) "the semigroup must generate ZZ^d"

  G = grading_group(kQ)

  # The shift must move the degrees of the Bass numbers at the maximal ideal
  # into Q up to cohomological degree i + d, d = dim Q: by [HM05, Lemma 4.5] a
  # summand k{a + F - Q} of J^j has non-zero Q-graded part once the summands
  # of Gamma_m J^{j + d - dim F} do, and dim F can be 0.
  depth = i + ambient_dimension(affine_semigroup(kQ))

  #compute irreducible resolution of shifted module
  if shift === :bound
    a_shift = compute_shift_bound(M, depth)
  elseif shift === :helm_miller
    a_shift = compute_shift(M, depth)
  elseif shift === :milp
    a_shift = compute_shift_milp(M, depth)
  elseif shift === :milp_bound
    a_shift = compute_shift_milp_bound(M, depth)
  else
    throw(ArgumentError("unknown shift strategy :$shift; expected :bound, :helm_miller, :milp, or :milp_bound"))
  end

  M_a = twist(M, -G(a_shift))
  irr_res = irreducible_resolution(M_a, i; check)

  #get injective modules up to cohomological degree i, i.e. J^0, J^1, ...,J^i
  inj_modules = Vector{InjMod}()
  for j in 1:min(i + 1, length(irr_res.irr_sums))
    shifted_comp = map(
      indec -> IndecInj(indec.face, indec.vector - a_shift), irr_res.irr_sums[j].indec_injectives
    )
    push!(inj_modules, InjMod(kQ, shifted_comp))
  end

  # the embedding M -> J^0 and the differentials d^k: J^k -> J^{k+1} as
  # monomial matrices (irr_res.cochain_maps[1] is the embedding, [k+2] is d^k)
  emb = _scalar_matrix(kQ, matrix(irr_res.cochain_maps[1]))
  cochain_maps = MonomialMatrix{InjMod}[]
  for k in 0:(length(inj_modules) - 2)
    A = _scalar_matrix(kQ, matrix(irr_res.cochain_maps[k + 2]))
    push!(cochain_maps, MonomialMatrix(A, inj_modules[k + 1], inj_modules[k + 2], k))
  end
  return InjRes(M, inj_modules, emb, cochain_maps, length(inj_modules) - 1, irr_res, a_shift)
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
    q = fdiv(bb[j], H[i, j])                   # floor division gives the canonical representative
    for k in 1:ncols(H)
      bb[k] -= q * H[i, k]
    end
  end
  return [Int(x) for x in bb]
end

@doc raw"""
    graded_bass_numbers(M::SubquoModule{<:MonoidAlgebraElem}, p::FaceQ, i::Int)

Return the graded Bass numbers $\mu^j(p_F, M)_a$ of `M` at the face `p` for
$0 \le j \le i$, which corresponds to the number of copies of $k\{a + F - Q\}$ in the $j$-th
term of a minimal injective resolution of `M`.
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

Check whether the injective resolution `res` is minimal by comparing the
multiplicities of its indecomposable injectives with the graded Bass numbers
of the resolved module, see `graded_bass_numbers`. With `verbose = true`
each mismatch is reported.
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

function injective_resolution(I::MonoidAlgebraIdeal, i::Int; shift::Symbol=:bound, check::Bool=true)
  return injective_resolution(quotient_ring_as_module(I), i; shift, check)
end

@doc raw"""
    injective_hull(M::SubquoModule{<:MonoidAlgebraElem})

Return the injective hull $E(M)$ of the finitely generated $\mathbb{Z}^d$-graded
module `M` as an `InjMod` together with the matrix defining the
embedding $M \hookrightarrow E(M)$.
"""
function injective_hull(M::SubquoModule{<:MonoidAlgebraElem})
  kQ = base_ring(M)
  @req _generates_lattice(kQ) "the semigroup must generate ZZ^d"

  R_Q = kQ.algebra
  G = grading_group(kQ)

  #compute irreducible hull of shifted module; the Bass numbers up to
  #cohomological degree d = dim Q control the summands of J^0, see
  #injective_resolution
  a_shift = compute_shift_bound(M, ambient_dimension(affine_semigroup(kQ)))
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

## Functions and types visible on the outside
export AffineSemigroup
export FaceQ
export IndecInj
export InjMod
export InjRes
export IrrRes
export IrrSum
export MonoidAlgebra
export MonoidAlgebraElem
export MonoidAlgebraIdeal
export MonomialMatrix
export affine_semigroup
export augmentation_map
export cochain_maps
export cohomological_degree
export degree_shift
export degrees_of_bass_numbers
export holes_module
export indecomposable_injectives
export injective_hull
export injective_modules
export injective_resolution
export irreducible_resolution
export irreducible_sums
export is_minimal
export local_cohomology
export local_cohomology_all
export monoid_algebra
export monoid_algebra_ideal
export monomial_matrix
export Q_graded_part
export sectors
export semigroup_generators
export zeroth_local_cohomology
