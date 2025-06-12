module InjectiveResolutions
using ..Oscar
using ..Oscar: IntegerUnion # add other things that are not exported from Oscar here
# functions with new methods
import ..Oscar:
  coefficient_ring,
  gens,
  ModuleGens,
  SubModuleOfFreeModule,
  is_zm_graded,
  grading_group,
  elem_type,
  singular_module, 
  singular_generators,
  cone,
  faces,
  hyperplanes,
  is_pointed,
  zonotope,
  degree,
  zero,
  one,
  inv,
  degree,
  _degree_fast,
  dim,
  is_subset,
  coordinates,
  intersect,
  primitive_generator,
  evaluate,
  coefficients,
  dim,
  is_normal,
  standard_basis,
  _build_sparse_row,
  normal_form,
  lift_std,
  sparse_row,
  syzygy_module,
  kernel,
  annihilator,
  twist,
  _saturation,
  free_resolution,
  _reduce, 
  oscar_assure,
  _graded_kernel,
  oscar_free_module,
  oscar_generators,
  singular_poly_ring,
  _simple_kernel,
  _extend_free_resolution


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

export irreducible_resolution
export irreducible_decomposition

export injective_resolution

export local_cohomology
export local_cohomology_all
export zeroth_local_cohomology

export MonoidAlgebra
export MonoidAlgebraIdeal
export MonoidAlgebraElem

#########################
# some composite types
#########################

include("MonoidAlgebra.jl")

struct IndecInj #indecomposable injective
  face::FaceQ
  vector::Vector{Int}
end
mutable struct InjMod #ZZ^d i-graded injective module over monoid algebra
  monoid_algebra::MonoidAlgebra
  indec_injectives::Vector{IndecInj}
  Q_graded_part::Union{SubquoModule,Nothing}

  function InjMod(A::MonoidAlgebra,I::Vector{IndecInj})
    return new(A,I,nothing)
  end
end

function Q_graded_part(I::InjMod)
  if I.Q_graded_part === nothing
    I.Q_graded_part = compute_Q_graded_part(I)
  end
  return I.Q_graded_part
end

function compute_Q_graded_part(I::InjMod)
  kQ = I.monoid_algebra
  irreducible_ideals = [_get_irreducible_ideal(kQ, J) for J in I.indec_injectives]
  irreducible_comp = [quotient_ring_as_module(Ji) for Ji in irreducible_ideals]
  return direct_sum(irreducible_comp...,task=:none)
end
struct IrrRes # irreducible resolution (including all computed data and the cochain complex)
  mod::SubquoModule
  irr_sums::Vector{InjMod}
  cochain_maps::Vector{SubQuoHom}
  cochain_complex::ComplexOfMorphisms # if sequence not exact return trivial cochain_complex (M0 -> M0)
end

#direct sum of two injective modules over the same monoid algebra
function +(I::InjMod, J::InjMod)
  @req I.monoid_algebra == J.monoid_algebra "monoid algebras not the same"
  return InjMod(I.monoid_algebra, vcat(I.indec_injectives, J.indec_injectives))
end

struct InjRes #ZZ^d-graded injective resolution
  mod::SubquoModule
  inj_mods::Vector{InjMod}
  cochain_maps::Vector{MatElem}
  upto::Int
  Q_graded_part::IrrRes
  shift::Vector{Int} #not needed
end

function Base.show(io::IO,J::InjMod)
  print(
    io, "injective module given by direct sum of ", length(J.indec_injectives)," indecomposable injectives"
  )
end

function Base.show(io::IO, ::MIME"text/plain", J::InjMod)
  println(io, "injective module given by direct sum of indecomposable injectives")
  for Ji in J.indec_injectives
    println(io, "  k{", Ji.vector, " + F - Q}, where p_F = ", Ji.face.prime)
  end
  print(
    io,
    "over ",
    J.monoid_algebra
  )
end

function Base.show(io::IO, ::MIME"text/plain", res::InjRes)
  println(io, "injective resolution ")
  println(io, "  ", join(["J^$i" for i in 0:res.upto], " -> "))
  println(io, "where ")
  for i in eachindex(res.inj_mods)
    j = i-1
    println(io, " J^$j = direct sum of")
    for Ji in res.inj_mods[i].indec_injectives
      println(io, "    k{", Ji.vector, " + F - Q}, where p_F = ", Ji.face.prime)
    end
  end
  println(io, "of ", res.mod)
  print(
    io,
    "over ",
    base_ring(res.mod)
  )
end

function Base.show(io::IO, ::MIME"text/plain", res::IrrRes)
  println(io, "irreducible resolution ")
  println(io, "  ", join(["W^$i" for i in 0:(length(res.irr_sums) - 1)], " -> "))
  println(io, "where ")
  for i in eachindex(res.irr_sums)
    j = i-1
    println(io, " W^$j = direct sum of")
    for Ji in res.irr_sums[i].indec_injectives
      println(io, "    k{", Ji.vector, " + F - Q}_Q, where p_F = ", Ji.face.prime)
    end
  end
  println(io, "of ", res.mod)
  print(
    io,
    "over ",
    base_ring(res.mod)
  )
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

  F = intersect(H.hyperplane, kQ.cone)

  #get faces of Q intersecting F only origin in Q
  D = Vector{Polyhedron}()
  for f in kQ.faces
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
        a in lattice_points(I + kQ.zonotope[1]) if (dim(intersect(convex_hull(a), PaF)) < 0)
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
function degrees_of_bass_numbers(M::SubquoModule{<:MonoidAlgebraElem}, i::Int) # now also for fin. gen. modules
  R_Q = base_ring(M)

  # residue field
  I_m = ideal(R_Q, gens(R_Q))
  k = quotient_ring_as_module(I_m)

  degrees = Vector{Vector{Int}}()
  for j in 0:i
    E = ext(k, M, j)
    for g in gens(E)
      push!(degrees, degree(Vector{Int}, g))
    end
  end
  return unique(degrees) # filter duplicates
end

@doc raw"""
    compute_shift(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)

Let $M$ be finitely generated $\mathbb{Z}^d$-graded module over a monoid algebra $k[Q]$. This function computes $a\in \mathbb{Z}^d$
such that all $\mathbb{Z}^d$-degrees of non-zero Bass numbers of $M(-a)$ lie in $Q$. 
"""
function compute_shift(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)
  # kQ = M.monoid_algebra
  kQ = base_ring(M)

  #get all degrees of non-zero Bass numbers up to cohomological degree i
  n_bass = degrees_of_bass_numbers(M, i)

  #sum of all primitive integer vectors along rays of Q
  c = zonotope(kQ)[2]

  j = 0
  while !all([is_subset(convex_hull(b), cone(kQ)) for b in n_bass]) #loop until all degrees of bass numbers lie in Q
    bass_ = [a_bass + c for a_bass in n_bass]
    n_bass = bass_
    j = j + 1
  end
  return j*c
end

@doc raw"""
    mod_quotient(M::SubquoModule,I::Ideal)

Computes the submodule

$(0 :_M I) := \{m \in M \mid m\cdot I = 0\}.$

# Examples
```jldoctest
julia> R_Q,(x,y) = graded_polynomial_ring(QQ,["x","y"];weights = [[1,0],[0,1]])
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
julia> R_Q,(x,y) = graded_polynomial_ring(QQ,["x","y"];weights = [[1,0],[0,1]])
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
function ZF_basis(M::SubquoModule{<:MonoidAlgebraElem}, p::FaceQ)
  kQ = base_ring(M)
  @assert kQ.algebra == base_ring(p.prime)

  T = elem_type(M)

  #compute quotient (0 :_M p)[\ZZ F]
  Np = mod_quotient(M, monoid_algebra_ideal(kQ,p.prime))[1]

  #initialize
  N = Np
  h_N = identity_map(N)
  B = Vector{T}() # empty vector of k[ZF]-basis 

  for g in gens(Np)
    N_g = sub(N, [h_N(g)])[1] #submodule of N =(0 :_M p_F)/(y0,...,yn) generated by g

    if annihilator(N_g) == MonoidAlgebraIdeal(kQ, p.prime) && !is_zero(N_g)
      push!(B, g)
    end
    M_B, _ = sub(Np, B)
    N, h_N = quo(Np, M_B) # update N
    if is_zero(N)
      break
    end
  end
  return filter(!is_zero, B)
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

@doc raw"""
    coefficients(N::SubquoModule, p_F::FaceQ)

Returns a subset Bp $\subseteq M$ and a $k$-matrix $\Lambda$ that defines an injective map

$(0 :_N p_F) \xrightarrow{\Lambda} \sum_{b\in Bp}k\{\deg(b) + F - Q\}.$

This fixes Algorithm 3.6. in [HM05](@cite).
"""
function coefficients(N::SubquoModule{T}, p_F::FaceQ) where {T <: MonoidAlgebraElem}
  kQ = base_ring(N)
  @assert base_ring(p_F.prime) == kQ.algebra

  #get the coefficient field
  k = coefficient_ring(kQ)

  # get socle degrees of indecomposable injectives kQ{a + F - Q}, i.e. compute a k[F]-basis of the localisation (0 :_M P_F)[ZZ F]
  Bp = ZF_basis(N, p_F)
  if is_empty(Bp)
    return [], zeros(kQ, 1, 1)
  end

  R = relations(N) #get all relations of N

  lambda = Vector{Vector{elem_type(k)}}()
  for b in Bp
    #get coefficient vector of w.r.t. generators of M
    b_amb = ambient_representative(b)
    _c_b = coordinates(N(b_amb)) #coordinates w.r.t. generators of M 
    c_b = [evaluate(_c_b[i], [1 for _ in 1:ngens(kQ)]) for i in 1:ngens(N)]
    # A possible alternative?
    # c_b = [only(AbstractAlgebra.coefficients(_c_b[i])) for i in 1:ngens(N)]

    #get all relevant generators of N, i.e., check (deg(b) + F) \cap (deg(g) + Q) ≠ ∅
    b_p = convex_hull(degree(Vector{Int}, b))
    G_b = Vector{SubquoModuleElem}()
    for g_N in filter(!is_zero, gens(N))
      g_p = convex_hull(degree(Vector{Int}, g_N))
      if dim(intersect(b_p + p_F.poly, g_p + kQ.cone)) >= 0
        push!(G_b, g_N)
      end
    end
    x_Gb = [monomial_basis(kQ, degree(g))[1] for g in G_b]
    _N = sub(ambient_free_module(N), [ambient_representative(g) for g in G_b])[1]

    #get all b-relevant relations w.r.t. F
    R_bF = Vector{FreeModElem{elem_type(kQ)}}()
    C_bF = Vector{Vector{elem_type(k)}}()
    for r in R
      #check (deg(b) + F)\cap (deg(r) + Q) ≠ ∅
      r_p = convex_hull(degree(Vector{Int}, r))
      if dim(intersect(b_p + p_F.poly, r_p + kQ.cone)) >= 0
        x_r = monomial_basis(kQ, degree(r))[1]
        a = lcm(x_Gb..., x_r)
        _r = (a//x_r).num*r #well-defined since a is lcm

        #can _r be written using generators in G_b
        if _r in _N #can _r be lifted?
          _c_r = coordinates(_N(_r))
          c_r = Vector{elem_type(k)}()
          for i in 1:ngens(N)
            j = findfirst(g -> g == N[i], G_b)
            if j !== nothing
              push!(c_r, evaluate(_c_r[j], [1 for _ in 1:ngens(kQ)]))
              # A possible alternative?
              # push!(c_r, only(AbstractAlgebra.coefficients(_c_r[j])))
            else
              push!(c_r, k())
            end
          end
          push!(C_bF, c_r)
          push!(R_bF, r)
        end
      end
    end

    #get kernel of coefficient matrix -> K
    if !is_empty(C_bF)
      _K = matrix(k, hcat(C_bF...))
      K = kernel(_K)
    else
      K = identity_matrix(k, ngens(N))
    end

    #get kernel of coefficient vector of b -> B
    _B = matrix(k, 1, length(c_b), c_b)
    B = kernel(transpose(_B))

    #get \lambda_b
    rows_K = [K[i, :] for i in 1:nrows(K)]
    l = rank(B)
    possible_rows = Vector{Vector{elem_type(k)}}()
    for c_K in rows_K
      B_K = vcat(B, matrix(k, 1, ngens(N), c_K))
      if rank(B_K)>l
        push!(possible_rows, c_K)
      end
    end
    push!(lambda, possible_rows[1])
  end
  return Bp, matrix(kQ,map(kQ, hcat(lambda...)))
end

@doc raw"""
    irreducible_hull(Mi::SubquoModule, kQ::MonoidAlgebra, j=0)

Return an irreducible hull of $M$.
"""
function irreducible_hull(Mi::SubquoModule{<:MonoidAlgebraElem}, j=0)
  kQ = base_ring(Mi)
  T = elem_type(kQ)

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
    M_sat = saturation((ideal(kQ, []) * Mi)[1], monoid_algebra_ideal(kQ,p.prime))
    if !is_zero(p.prime) && !is_zero(M_sat)
      N, _ = quo(Mi, M_sat)
    end
    if is_zero(N)
      break
    end
  end
  # return summands, hcat(lambda...)
  return InjMod(kQ,summands), hcat(lambda...)
end

@doc raw"""
    irreducible_decomposition(I::MonoidAlgebraIdeal)

Return an irreducible decomposition of $I$. 

!!! note
    The monoid algebra $k[Q]$ must be normal. 

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
  @req is_normal(kQ) "monoid algebra must be normal"

  J, _ = irreducible_hull(quotient_ring_as_module(I))
  return [_get_irreducible_ideal(kQ, I) for I in J.indec_injectives]
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

@doc raw"""
    irreducible_resolution(M::SubquoModule{<:MonoidAlgebraElem}, i::Int = 0)

Return an irreducible resolution of $M$.

!!! note
    The monoid algebra $k[Q]$ must be normal. 

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
function irreducible_resolution(M::SubquoModule{<:MonoidAlgebraElem}, i::Int=0)
  kQ = base_ring(M)
  @req is_normal(kQ) "monoid algebra must be normal"

  R_Q = kQ.algebra
  Mi = M # current module in resolution

  #initialize
  gi = identity_map(Mi)
  irreducible_sums = Vector{InjMod}()
  cochain_maps = Vector{SubQuoHom}()

  j = 1
  while !is_zero(Mi) #until cokernel Mi is zero
    #compute irreducible hull
    Ji, _lambda = irreducible_hull(Mi, j)
    
    #get Q-graded part
    Wi = Q_graded_part(Ji)

    #multiply rows of lambda by degrees of generators of Mi
    m, n = size(_lambda)
    lambda = zero(_lambda)
    for i in 1:m
      for j in 1:n
        lambda[i, j] = monomial_basis(R_Q, degree(Mi[i]))[1] * _lambda[i, j]
      end
    end

    #define injective map Mi -> Wi
    fi = hom(Mi, Wi, matrix(lambda))

    #get boundary map W{i-1} -> Wi
    hi = gi*fi

    #compute cokernel and then simplify
    Mi_, gi_ = quo(Wi, image(hi)[1]) #cokernel
    Mi, ji = prune_with_map(Mi_)
    gi = gi_*inv(ji)
    #Mi = Mi_
    #gi = gi_

    any(is_zero, relations(Mi)) && error("there must not be trivial relations")
    #fix modules with "zero" relations
    if length(filter(is_zero, relations(Mi))) > 0
      Mi, h = fix_module(Mi)
      gi = gi*h
    end

    push!(irreducible_sums, Ji)
    push!(cochain_maps, hi)

    # end at cohomological degree i
    if i > 0 && j == i
      break
    end
    j = j + 1
  end

  #get cochain complex
  C = cochain_complex(cochain_maps)

  return IrrRes(M, irreducible_sums, cochain_maps, C)
end

@doc raw"""
    injective_resolution(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)

Return an injective resolution of $M$ up to cohomological degree i.

!!! note
    The monoid algebra $k[Q]$ must be normal. 

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
function injective_resolution(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)
  kQ = base_ring(M)
  @req is_normal(kQ) "monoid algebra must be normal"

  G = grading_group(kQ)

  #compute irreducible resolution of shifted module
  a_shift = compute_shift(M, i+1)
  M_a = twist(M, -G(a_shift))
  irr_res = irreducible_resolution(M_a)

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
    matrix(irr_res.cochain_maps[k]) for k in eachindex(irr_res.cochain_maps) if 1 < k <= i+1
  ]
  return InjRes(M, inj_modules, cochain_maps, length(inj_modules)-1, irr_res, a_shift)
end

@doc raw"""
    injective_resolution(I::MonoidAlgebraIdeal,i::Int)
    
Return an injective resolution of $M = k[Q]/I$ up to cohomological degree i.  

!!! note
    The monoid algebra $k[Q]$ must be normal. 
"""
function injective_resolution(I::MonoidAlgebraIdeal, i::Int)
  return injective_resolution(quotient_ring_as_module(I), i)
end

#= to be enabled, once #861 is merged in Singular.jl
@attr Int function dim(M::SubquoModule{T}) where {CT<:FieldElem, T<:MPolyRingElem}
  F = ambient_free_module(M)

  if !all(repres(v) == F[i] for (i, v) in enumerate(gens(M)))
    MM, _ = present_as_cokernel(M)
    return dim(MM)
  end

  gb = groebner_basis(M.quo)
  return Singular.dimension(singular_generators(gb))
end
=#

#get Krull dimension of module
# TODO: For modules over polynomial rings this should make 
# use of the `dim` in Singular. But this does not seem 
# to be available as of yet.
@attr Union{Int,NegInf} function dim(M::ModuleFP)
  ann = annihilator(M)
  return dim(ann)
end


@doc raw"""
    is_normal(A::MonoidAlgebra{<:FieldElem, <:MPolyQuoRing})

Test if the given monoid algebra is normal by testing first the S2 and then the
R1 condition.

# Examples
```jldoctest
julia> A = monoid_algebra([[4,0],[3,1],[1,3],[0,4]],QQ)
monoid algebra over rational field with cone of dimension 2

julia> is_normal(A)
false
```
"""
@attr Bool function is_normal(A::MonoidAlgebra{<:FieldElem, <:MPolyQuoRing})
  # Implementation adapted from 
  #
  #    M2/Macaulay2/packages/IntegralClosure.m2,
  #
  # line 666 ff. of https://github.com/Macaulay2/M2/blob/2565455411d15a3386204aa62a00e20ee5c0e99f/M2/Macaulay2/packages/IntegralClosure.m2 
  # on Apr 25, 2025.
  R = A.algebra::MPolyQuoRing
  R_B = base_ring(R)
  I = modulus(R)
  M = quotient_ring_as_module(I)
  n = codim(I)

  # Check the S2 condition
  test_range = 0:(dim(R_B) - n - 2)

  for j in test_range
    # Check if codimension of Ext^{j+n+1} is at least j+n+3
    # get ext:
    E = ext(M,graded_free_module(R_B, 1),j + n + 1)
    # work around issue https://github.com/oscar-system/Oscar.jl/issues/4884
    if is_zero(E)
      d = -1
    else 
      d = dim(E)
    end
    cod = dim(R_B) - d
    if cod < j+n+3
        return false
    end               # S2 condition not satisfied
  end

  Jac = ideal(R, minors(map_entries(R, jacobian_matrix(gens(modulus(R)))), n))  # Compute minors of the Jacobian
  d = dim(Jac)                   # Get dimension of the Jacobian
  d < 0 && (d = -Inf)            # Handle negative dimensions
  return (dim(R) - d >= 2)       # Check the condition
end

is_normal(A::MonoidAlgebra{<:FieldElem, <:MPolyRing}) = true


# import local cohomology functions
include("LocalCohomology.jl")
include("ModuleFunctionality.jl")

end # module InjectiveResolutions

using .InjectiveResolutions

## Functions visible on the outside
export monoid_algebra
export faces
export hyperplanes

export irreducible_resolution
export irreducible_decomposition

export injective_resolution

export local_cohomology
export local_cohomology_all
export zeroth_local_cohomology

export MonoidAlgebra
export MonoidAlgebraIdeal
export MonoidAlgebraElem

