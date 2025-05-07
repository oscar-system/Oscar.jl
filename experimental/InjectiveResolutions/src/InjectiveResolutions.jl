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
  singular_assure,
  oscar_assure, 
  _graded_kernel,
  oscar_generators,
  singular_poly_ring,
  _simple_kernel,
  _extend_free_resolution, 
  free_show


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
export get_monoid_algebra
export monoid_algebra_from_lattice
export monoid_algebra
export monoid_algebra_ideal
export monoid_algebra_module

export irreducible_res
export irreducible_dec

export injective_res

#########################
# some composite types
#########################

struct FaceQ # face of semigroup  
  prime::Union{MPolyIdeal,MPolyQuoIdeal} #homogeneous prime corresponding to face
  poly::Polyhedron  #face as polyhedron
end

struct HyperplaneQ # a hyperplane bounding the cone RR_{\geq 0}Q
  hyperplane::Polyhedron
  A::Matrix{Int}
  b::Vector{Int}
end

@attributes mutable struct MonoidAlgebra{CoeffType,AlgebraType} <: Ring # monoid algebra with associated data
  algebra::AlgebraType
  pointed::Union{Nothing,Bool} #cone pointed
  polyhedral_cone::Cone{QQFieldElem}
  cone::Polyhedron
  faces::Vector{FaceQ}
  hyperplanes::Vector{HyperplaneQ}
  zonotope::Tuple{Polyhedron,Vector{Int}}

  function MonoidAlgebra(
    A::AlgebraType; check::Bool=true
  ) where {AlgebraType<:Union{MPolyRing,MPolyQuoRing}}
    #check if monoid algebra
    @check is_zm_graded(A) "given algebra is not ZZ^d-graded"
    gg_Q = grading_group(A)
    @check is_free(gg_Q) && is_abelian(gg_Q) "given algebra is not a monoid algebra"
    kk = coefficient_ring(A)
    result = new{elem_type(kk),AlgebraType}(A, nothing)
  end
end

is_graded(A::MonoidAlgebra) = true
is_zm_graded(A::MonoidAlgebra) = true

function polyhedral_cone(A::MonoidAlgebra)
  if !isdefined(A, :polyhedral_cone)
    A.polyhedral_cone = get_polyhedral_cone(A.algebra)
  end
  return A.polyhedral_cone
end

function cone(A::MonoidAlgebra)
  if !isdefined(A, :cone)
    A.cone = polyhedron(polyhedral_cone(A))
  end
  return A.cone
end

function faces(A::MonoidAlgebra)
  if !isdefined(A, :faces)
    A.faces = get_faces_of_polyhedral_cone(A.algebra, zonotope(A)[1], cone(A))
  end
  return A.faces
end

function hyperplanes(A::MonoidAlgebra)
  if !isdefined(A, :hyperplanes)
    A.hyperplanes = get_bounding_hyperplanes(cone(A))
  end
  return A.hyperplanes
end

function is_pointed(A::MonoidAlgebra)
  if isnothing(A.pointed)
    A.pointed = is_pointed(polyhedral_cone(A))
  end
  return A.pointed::Bool
end

function zonotope(A::MonoidAlgebra)
  if !isdefined(A, :zonotope)
    A.zonotope = get_zonotope(cone(A))
  end
  return A.zonotope
end

coefficient_ring(A::MonoidAlgebra) = coefficient_ring(A.algebra)
number_of_variables(A::MonoidAlgebra) = ngens(A.algebra)

### Elements of MonoidAlgebras
mutable struct MonoidAlgebraElem{CoeffType,ParentType} <: RingElem
  parent::ParentType
  elem::RingElem
  # TODO: Do we want to store additional information on the elements here?

  function MonoidAlgebraElem(
    A::ParentType
  ) where {CoeffType,ParentType<:MonoidAlgebra{CoeffType}}
    return new{CoeffType,ParentType}(A)
  end

  function MonoidAlgebraElem(
    A::ParentType,
    a::RingElem;
    check::Bool=true,
  ) where {CoeffType,ParentType<:MonoidAlgebra{CoeffType}}
    @check parent(a) === A.algebra
    return new{CoeffType,ParentType}(A, a)
  end
end

function ==(A::MonoidAlgebra, B::MonoidAlgebra)
  return A.algebra == B.algebra
end

parent(a::MonoidAlgebraElem) = a.parent

elem_type(::Type{T}) where {CoeffType,T<:MonoidAlgebra{CoeffType}} = MonoidAlgebraElem{
  CoeffType,T
}

function degree(
    a::MonoidAlgebraElem{<:FieldElem, PT};
    check::Bool=true
  ) where {RT <: MPolyRing, PT <: MonoidAlgebra{<:FieldElem, RT}}
  !isdefined(a, :elem) && return zero(grading_group(parent(a)))
  return degree(underlying_element(a); check)
end

grading_group(A::MonoidAlgebra) = grading_group(A.algebra)

function underlying_element(a::MonoidAlgebraElem)
  if !isdefined(a, :elem)
    a.elem = zero(parent(a).algebra)
  end
  return a.elem::elem_type(parent(a).algebra)
end

### implementation of some common ring functionality
function (A::MonoidAlgebra)()
  return MonoidAlgebraElem(A)
end

function zero(A::MonoidAlgebra)
  return A()
end

function one(A::MonoidAlgebra)
  return MonoidAlgebraElem(A, one(A.algebra); check=false)
end

function +(a::MonoidAlgebraElem, b::MonoidAlgebraElem)
  A = parent(a)
  @assert A === parent(b)
  return MonoidAlgebraElem(A, underlying_element(a) + underlying_element(b); check=false)
end

function -(a::MonoidAlgebraElem, b::MonoidAlgebraElem)
  A = parent(a)
  @assert A === parent(b)
  return MonoidAlgebraElem(A, underlying_element(a) - underlying_element(b); check=false)
end

function ==(a::MonoidAlgebraElem, b::MonoidAlgebraElem)
  A = parent(a)
  @assert A === parent(b)
  return underlying_element(a) == underlying_element(b)
end

function *(a::MonoidAlgebraElem, b::MonoidAlgebraElem)
  A = parent(a)
  @assert A === parent(b)
  return MonoidAlgebraElem(A, underlying_element(a) * underlying_element(b); check=false)
end

function *(c::CoeffType, b::MonoidAlgebraElem{CoeffType}) where {CoeffType}
  A = parent(b)
  @assert coefficient_ring(A) === parent(c)
  return MonoidAlgebraElem(A, c*underlying_element(b); check=false)
end

function *(c::IntegerUnion, b::MonoidAlgebraElem)
  A = parent(b)
  kk = coefficient_ring(A)
  return MonoidAlgebraElem(A, kk(c)*underlying_element(b); check=false)
end

#=
# TODO: should not be necessary! remove?
function *(m::MPolyDecRingElem{CoeffType},a::MonoidAlgebraElem{CoeffType}) where {CoeffType}
  kQ = parent(a)
  return kQ(m*underlying_element(a))
end
=#

function -(b::MonoidAlgebraElem)
  A = parent(b)
  return MonoidAlgebraElem(A, -underlying_element(b); check=false)
end

function (A::MonoidAlgebra{CoeffType})(c::CoeffType) where {CoeffType}
  return MonoidAlgebraElem(A, A.algebra(c))
end

function (A::MonoidAlgebra)(c::Any)
  return MonoidAlgebraElem(A, A.algebra(c))
end

function (A::MonoidAlgebra)(a::MonoidAlgebraElem)
  @assert parent(a) === A
  return a
end

function deepcopy_internal(a::MonoidAlgebraElem, dict::IdDict)
  return parent(a)(deepcopy_internal(underlying_element(a), dict))
end

#=
# TODO: should not be necessary! remove?
function Base.in(a::MonoidAlgebraElem,M::SubquoModule{<:MonoidAlgebraElem})
  @assert parent(a) == base_ring(M)
  return underlying_element(a) in M
end
=#

is_unit(a::MonoidAlgebraElem) = is_unit(underlying_element(a))

function inv(a::MonoidAlgebraElem)
  return MonoidAlgebraElem(parent(a), inv(underlying_element(a)))
end

monomial_basis(A::MonoidAlgebra, g::FinGenAbGroupElem) = monomial_basis(A.algebra, g)
#elem_type(A)[A(x) for x in monomial_basis(A.algebra,g)]

evaluate(a::MonoidAlgebraElem, vals::Vector) = evaluate(underlying_element(a),vals) 

is_homogeneous(a::MonoidAlgebraElem) = is_homogeneous(underlying_element(a))

function degree(
        a::MonoidAlgebraElem{<:FieldElem, PT};
        check::Bool=true
    ) where {RT <: MPolyQuoRing, PT <: MonoidAlgebra{<:FieldElem, RT}}
  !isdefined(a, :elem) && return zero(grading_group(parent(a)))
  _a = underlying_element(a)
  simplify(_a)
  @req !iszero(_a) "Element must be non-zero"
  return degree(_a.f; check)
end

function _degree_fast(a::MonoidAlgebraElem)
  return _degree_fast(underlying_element(a))
end

function degree(::Type{Vector{Int}}, a::MonoidAlgebraElem; check::Bool=true)
  _a = underlying_element(a)
  @assert is_zm_graded((base_ring(parent(_a))))
  d = degree(_a; check)
  return Int[d[i] for i=1:ngens(parent(d))]
end

function dim(A::MonoidAlgebra)
  return dim(A.algebra)
end

#=
# TODO: Do we want this?
function AbstractAlgebra.coefficients(a::MonoidAlgebraElem)
  return AbstractAlgebra.coefficients(underlying_element(a))
end

function AbstractAlgebra.coefficients(a::MPolyQuoRingElem)
  return AbstractAlgebra.coefficients(a.f)
end
=#


# TODO: Fix up promote rules?
AbstractAlgebra.promote_rule(
                             ::Type{CoeffType}, ::Type{T}
                            ) where {CoeffType,T<:MonoidAlgebraElem{CoeffType}} = T

# TODO: We would like to use the parametrization of the `MonoidAlgebra` 
# with the type of its underlying ring directly. But this parameter only
# provides the type of the ring, not its elements. So we have to work 
# as if `MonoidAlgebra` was not parametrized.
AbstractAlgebra.promote_rule(
                             ::Type{MPolyDecRingElem{CoeffType, T}}, ::Type{U}
                            ) where {CoeffType, T, U<:MonoidAlgebraElem{CoeffType}} = U

AbstractAlgebra.promote_rule(
                             ::Type{MPolyQuoRingElem{PolyType}}, ::Type{T}
                            ) where {CoeffType, PolyType<:MPolyDecRingElem{CoeffType}, T<:MonoidAlgebraElem{CoeffType}} = T

AbstractAlgebra.promote_rule(::Type{Int}, ::Type{T}) where {T<:MonoidAlgebraElem} = T
AbstractAlgebra.promote_rule(::Type{ZZRingElem}, ::Type{T}) where {T<:MonoidAlgebraElem} = T

gens(A::MonoidAlgebra) = [MonoidAlgebraElem(A, x) for x in gens(A.algebra)]
number_of_generators(A::MonoidAlgebra) = ngens(A.algebra)
getindex(A::MonoidAlgebra, i::Int) = MonoidAlgebraElem(A, A.algebra[i])

function Base.show(io::IO, a::MonoidAlgebraElem)
  print(io, underlying_element(a))
end

parent_type(
  ::Type{ElemType}
) where {ParentType,CoeffType,ElemType<:MonoidAlgebraElem{CoeffType,ParentType}} =
  ParentType

# TODO: Finish implementation of the ring interface! 

### Ideals over `MonoidAlgebra`s
mutable struct MonoidAlgebraIdeal{ElemType} <: Ideal{ElemType}
  monoid_algebra::MonoidAlgebra
  gens::Vector{ElemType}
  ideal::Ideal

  function MonoidAlgebraIdeal(A::MonoidAlgebra, v::Vector{T}) where {T<:MonoidAlgebraElem}
    @assert all(parent(x) === A for x in v)
    return new{T}(A, v)
  end

  # constructor from an `underlying_ideal`
  function MonoidAlgebraIdeal(A::MonoidAlgebra, I::Ideal)
    @assert base_ring(I) === A.algebra
    return new{elem_type(A)}(A, elem_type(A)[A(x) for x in gens(I)], I)
  end
end

function base_ring(I::MonoidAlgebraIdeal{ElemType}) where {ElemType}
  return I.monoid_algebra::parent_type(ElemType)
end

function gens(I::MonoidAlgebraIdeal{ElemType}) where {ElemType}
  return I.gens::Vector{ElemType}
end

function underlying_ideal(I::MonoidAlgebraIdeal)
  if !isdefined(I, :ideal)
    I.ideal = ideal(base_ring(I).algebra, [underlying_element(x) for x in gens(I)])
  end
  return I.ideal::ideal_type(base_ring(I).algebra)
end

# A sample for how to extend functionality via deflection to the underlying ideal
# and wrapping the result.
function radical(I::MonoidAlgebraIdeal)
  return MonoidAlgebraIdeal(base_ring(I), radical(underlying_ideal(I)))
end

dim(I::MonoidAlgebraIdeal) = dim(underlying_ideal(I))

# some generic functionality which should probably be elsewhere
function is_subset(I::Ideal, J::Ideal)
  return all(x in J for x in gens(I))
end

function Base.:(==)(I::Ideal, J::Ideal)
  return is_subset(I, J) && is_subset(J, I)
end

function Base.:*(I::Ideal, J::Ideal)
  @assert base_ring(I) === base_ring(J)
  return ideal(base_ring(I), [x*y for x in gens(I) for y in gens(J)])
end

# user facing constructor
ideal(A::MonoidAlgebra, v::Vector) = MonoidAlgebraIdeal(A, elem_type(A)[A(x) for x in v])

function Base.in(a::MonoidAlgebraElem, I::MonoidAlgebraIdeal)
  return underlying_element(a) in underlying_ideal(I)
end

function coordinates(a::MonoidAlgebraElem, I::MonoidAlgebraIdeal)
  return coordinates(underlying_element(a), underlying_ideal(I))
end

function intersect(a::MonoidAlgebraIdeal, b::MonoidAlgebraIdeal...)
  kQ = base_ring(a)
  for g in b
    @req base_ring(g) == kQ "base rings must match"
  end
  as = Singular.intersection(singular_generators(a), [singular_generators(g) for g in b]...)
  return MonoidAlgebraIdeal(kQ, [kQ(x) for x in gens(as)])
end

function Base.intersect(V::Vector{T}) where {T <: MonoidAlgebraIdeal}
  @assert length(V) != 0
  length(V) == 1 && return V[1]

  return intersect(V[1], V[2:end]...)
end

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
  surjections::Vector{SubQuoHom}
  inclusions::Vector{SubQuoHom}
  cokernels::Vector{SubquoModule}
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

function Base.show(io::IO, ::MIME"text/plain", kQ::MonoidAlgebra)
  println(io, "monoid algebra over ", lowercase(string(coefficient_ring(kQ))),
    " with cone of dimension $(dim(cone(kQ)))",
  )
  #=
  println(io,"Monoid algebra k[Q] over ", lowercase(string(coefficient_ring(kQ.algebra))), " generated by ", join(gens(kQ.algebra), ", "), " quotient by ideal (", join(gens(modulus(kQ.algebra)),", "),"), ")
  print(io, "where Q is a subset of ZZ^", ambient_dim(kQ.cone)," and the cone RR_{>= 0}Q has dimension ", dim(kQ.cone), ".")
  =#
end

function Base.show(io::IO, kQ::MonoidAlgebra)
  println(io, "monoid algebra over ", lowercase(string(coefficient_ring(kQ))),
    " with cone of dimension $(dim(cone(kQ)))",
  )
end

function Base.show(io::IO, ::MIME"text/plain", I::MonoidAlgebraIdeal)
  println(
    io, "ideal over $(base_ring(I)) generated by "*join(["$(x)" for x in gens(I)], ", ")
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
  println(io, "with maps")
  println(io, " J^{j-1} -> J^j for 1 <= j <= ", res.upto)
  println(io, "of ", res.mod)
  println(
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
  println(io, "with maps")
  println(io, " W^{i-1} -> W^i for 1 <= i <= ", length(res.irr_sums)-1)
  println(io, "of ", res.mod)
  println(
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
  if is_terse(io)
    print(io, "k{", Ji.vector, " + F - Q}, where p_F = ", Ji.face.prime)
  else
    print(
      io, "indecomposable injective k{", Ji.vector, " + F - Q}, where p_F = ", Ji.face.prime
    )
  end
end

function monoid_algebra(A::Union{MPolyRing,MPolyQuoRing})
  return MonoidAlgebra(A)
end

@doc raw"""
    monoid_algebra_from_lattice(B::Matrix{Int},k::Field)

Given a integer matrix $M_Q\in \mathbb{Z}^{d\times n}$ return the monoid algebra
\[k[x_1,\dots,x_n]/I_L,\]
where $I_L = (\bold{x}^u - \bold{x}^v \mid u,v \in \mathbb{N}^d \text{ with } u -v \in L)$ is the lattice ideal of the lattice $L$ generated by the columns $v_1,\dots,v_n$ of $M_Q$.

# Examples
```jldoctest
julia> M_Q = [1 0; 0 1]
2×2 Matrix{Int64}:
 1  0
 0  1

julia> monoid_algebra_from_lattice(M_Q,QQ)
monoid algebra over rational field with cone of dimension 2
```
"""
function monoid_algebra_from_lattice(M_Q::Matrix{Int}, k::Field)
  d = size(M_Q, 1)

  # construct k[t_1,...,t_d]
  t_vars = [Symbol("t_$i") for i in 1:d]
  T, t = graded_polynomial_ring(k, t_vars)

  # construct k[x_1,...,x_n] where n is the number of columns/generators
  x_vars = [Symbol("x_$i") for i in 1:size(M_Q, 2)]
  R, _ = graded_polynomial_ring(
    k, x_vars; weights=[Vector(row) for row in eachcol(M_Q)]
  )

  # construct map x_i \to t^(a_i)
  targ = [prod(t[j]^M_Q[j, i] for j in 1:d) for i in 1:size(M_Q, 2)]
  map_T_R = hom(R, T, targ)

  # return monoid algebra 
  return MonoidAlgebra(quo(R, ideal(kernel(map_T_R)[1]))[1])
end

@doc raw"""
    monoid_algebra_from_lattice(V_Q::Vector{Vector{Int}},k::Field)

Given a finite number of vectors $v_1,\dots,v_n$ in $\mathbb{Z}^d$ return the monoid algebra
\[k[x_1,\dots,x_n]/I_L,\]
where $I_L = (\bold{x}^u - \bold{x}^v \mid u,v \in \mathbb{N}^d \text{ with } u -v \in L)$ is the lattice ideal of the lattice $L$ generated by $v_1,\dots,v_n$.

# Examples
```jldoctest
julia> monoid_algebra_from_lattice([[0,1],[1,1],[2,1]],QQ)
monoid algebra over rational field with cone of dimension 2
```
"""
function monoid_algebra_from_lattice(V_Q::Vector{Vector{Int}}, k::Field)
  return monoid_algebra_from_lattice(Matrix{Int}(transpose(matrix(V_Q))), k)
end

function monoid_algebra_ideal(kQ::MonoidAlgebra, I::Ideal)
  @req base_ring(I) == kQ.algebra "base rings do not match"
  return MonoidAlgebraIdeal(kQ, I)
end

function Oscar.quotient_ring_as_module(I::MonoidAlgebraIdeal)
  R = base_ring(I)
  F = graded_free_module(R,1)
  e1 = F[1]
  return quo_object(F, [x * e1 for x = gens(I)]) 
end

@doc raw"""
  prime_to_face(kQ::Union{MPolyRing,MPolyQuoRing}, zonotope::Polyhedron, F::Polyhedron)

Let kQ be a monoid algebra over some semigroup $Q$. Given a face $F$ of the cone $C = \RR_{\geq 0}Q$ of the monoid algebra kQ,
return the corresponding homogeneous prime ideal
\[P_F = k\{Q\setminus F\}.\]
The zonotope is the Minkowski sum of all primitive integer vectors along rays of $Q$. 
"""
function prime_to_face(
  kQ::Union{MPolyRing,MPolyQuoRing}, zonotope::Polyhedron, F::Polyhedron
)
  gens = Vector{MPolyDecRingElem}()
  for lp in lattice_points(zonotope)
    if !(lp in F) #check if lattice point is in F
      a_v = Vector{ZZRingElem}()
      for a_p in lp # PointVector -> Vector
        push!(a_v, a_p)
      end
      push!(gens, monomial_basis(kQ, a_v)[1])
    end
  end
  return ideal(kQ, gens)
end

# given a monoid algebra, this function returns the corresponding polyhedral cone
# INPUT:    monoid algebra k[Q]
# OUTPUT:   polyhedral cone \RR_{\geq 0}Q
function get_polyhedral_cone(R::Union{MPolyDecRing,MPolyQuoRing})
  D = [degree(Vector{Int}, g) for g in gens(R)]
  return positive_hull(D)
end

# (works more general for polyhedra)
# given a polyhedral cone, this function returns its faces
# INPUT:    polyhedral cone C
# OUTPUT:   set of faces of C
function get_faces_of_polyhedral_cone(
  kQ::Union{MPolyRing,MPolyQuoRing}, zonotope::Polyhedron, P::Polyhedron
)
  P_faces = Vector{Polyhedron}()
  for i in 0:dim(P)
    in
    append!(P_faces, faces(P, i))
  end
  return [FaceQ(prime_to_face(kQ, zonotope, F), F) for F in P_faces]
end

# given a polyhedral cone, this function returns the hyperplanes bounding it
# INPUT:    polyhedral cone C
# OUTPUT:   set of hyperplanes bounding C
function get_bounding_hyperplanes(P::Polyhedron)
  hyperplanes = Vector{HyperplaneQ}()
  for f in facets(Polyhedron, P)
    hyperplane = polyhedron(affine_hull(f)[1])
    A, b = get_hyperplane_H_presentation(hyperplane)
    push!(hyperplanes, HyperplaneQ(hyperplane, A, b))
  end
  return hyperplanes
end

# TODO: Is this an admissible signature for a method of this function?
# Should this be moved to the PolyhedralGeometry section?
# This returns the primitive generator as a `Vector{Int}`.
function primitive_generator(::Type{Int}, r::AbstractVector{T}) where {T<:RationalUnion}
  return Vector{Int}(first(primitive_generator_with_scaling_factor(Int, r)))
end

# This method returns a triple `(v, num, den)` where `v` is a `Vector{Int}` 
# and `num` and `den` are both `Int`s so that `num//den` is the scaling factor.
function primitive_generator_with_scaling_factor(
    ::Type{Int},
    r::AbstractVector{T}
  ) where {T<:RationalUnion}
  @req !is_zero(r) "input must not be a zero vector"
  first_scaling_factor = lcm(denominator.(r))
  result = Int[Int(numerator(a)*divexact(first_scaling_factor, denominator(a))) for a in r]
  g = gcd(result)
  result = Int[divexact(a, g) for a in result]
  return Tuple{Vector{Int}, Int, Int}((result, first_scaling_factor, g))
end

# given a hyperplane, return the H-presentation of it
# INPUT:    hyperplane h
# OUTPUT:   matrix A, vector b corresponding to Ax \leq b which defines h
function get_hyperplane_H_presentation(h::Polyhedron)
  aff_hull = affine_hull(h).Obj.pm_polytope.AFFINE_HULL
  _M = Matrix{Rational}(aff_hull)
  M = hcat(map(row -> reshape(primitive_generator(Int, row), 1, :), eachrow(_M))...)
  A = [M[:, 2:n_columns(M)]; -M[:, 2:n_columns(M)]]
  b = [M[:, 1]; -M[:1]]
  return A, b
end

# given a polyhedral cone, return the zonotope as in Lemma 3.10 in [HM05]
# INPUT:    polyhedral cone C  
# OUTPUT:   zonotope, sum of primitive integer vector ong rays of C 
function get_zonotope(P::Polyhedron)
  d = ambient_dim(P)
  P_rays = [primitive_generator(Int, Vector(r)) for r in rays(P)]

  c = zeros(Int, d)
  zonotope = convex_hull(zeros(Int, d))
  for r in P_rays
    zonotope = zonotope + convex_hull([zeros(Int, 1, d); reshape(r, 1, length(r))])
    c = c + r
  end
  return zonotope, c
end

@doc raw"""
  generators_W_H(kQ::MonoidAlgebra, H::HyperplaneQ, a::Vector{Int})

Given a monoid algebra  kQ = $k[Q]$, a hyperplane $H$ that bounds the polyhedral cone $\RR_{\geq 0}Q$ and a vector
$a \in \mathbb{Z}^d$, return a finite set $B$ such that
\[(x^b \mid b \in B) \cong k\{(a + H_+^\circ)\cap Q\}.\]
This is Algorithm 3.11. in~\cite{HM05}.
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

Computes the submodule \[(0 :_M I) := \{m \in M \mid m\cdot I = 0\}.\]

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

julia> mod_quotient(M,m)
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

Compute the saturation \[(0 :_M I^\infty) := \{m \in M \mid m\cdot I^n = 0\text{ for some }n\in \NN_{>0}\}.\]

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

julia> mod_saturate(M,m)
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

Let $p = k\{Q\setminus F\}$ for some face $F$. This functions computes a $k[\mathbb{Z}F]$-basis of the quotient $(0 :_M p)[\mathhbb{Z}F] = \{m \in M \mid m\cdot p = 0\}[\mathbb{Z}F]$.
"""
function ZF_basis(M::SubquoModule{<:MonoidAlgebraElem}, p::FaceQ)
  kQ = base_ring(M)
  @assert kQ.algebra == base_ring(p.prime)

  T = elem_type(M)

  #compute quotient (0 :_M p)[\ZZ F]
  Np = mod_quotient(M, monoid_algebra_ideal(kQ,p.prime))[1]

  #initilalize
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
  coefficients(N::SubquoModule, p_F::FaceQ, kQ::MonoidAlgebra)

Returns a subset Bp $\subseteq M$ and a $k$-matrix $\Lambda$ that defines an injective map
\[(0 :_N p_F) \xrightarrow{\Lambda} \sum_{b\in Bp}k\{\deg(b) + F - Q\}.\]
This fixes Algorithm 3.6. in~\cite{HM05}.
"""
function coefficients(N::SubquoModule{<:MonoidAlgebraElem}, p_F::FaceQ)
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
  return Bp, map(kQ, hcat(lambda...))
end

@doc raw"""
  irreducible_hull(Mi::SubquoModule, kQ::MonoidAlgebra, j=0)

Returns an irreducible hull of a $\mathbb{Z}^d$-graded $k[Q]$-module M. It consists of indecomposable injectives $J_1,...,J_k$ and a $k$-matrix $\Lamdba$
such that $Mi \xhookrightarrow{\Lambda} \sum{i=1}^k J_i$. 
"""
function irreducible_hull(Mi::SubquoModule{<:MonoidAlgebraElem}, j=0)
  kQ = base_ring(Mi)
  T = elem_type(kQ)

  #initilalize
  N = Mi
  summands = Vector{IndecInj}()
  lambda = Vector{Matrix{T}}()

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
    irreducible_dec(I::MonoidAlgebraIdeal)

Let $R$ be a monoid algebra and $I\subset R$ an ideal. This function computes an irreducible decomposition of $I$, i.e.
\[I = W_1 \cap \dots \cap W_r,\]
where $W_i\subset R$ are irreducible ideals for all $i$.   

# Examples
```jldoctest
julia> R_Q,(x,y) = graded_polynomial_ring(QQ,["x","y"];weights = [[1,0],[0,1]])
(Graded multivariate polynomial ring in 2 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y])

julia> kQ = monoid_algebra(R_Q)
monoid algebra over rational field with cone of dimension 2


julia> I = ideal(kQ,[x^4,x^2*y^2,y^4])
ideal over monoid algebra over rational field with cone of dimension 2
 generated by x^4, x^2*y^2, y^4


julia> irreducible_dec(I)
2-element Vector{MPolyIdeal{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 Ideal (x^2, x^2*y, y^4, x*y^4)
 Ideal (x^4, x^4*y, y^2, x*y^2)
```
"""
function irreducible_dec(I::MonoidAlgebraIdeal)
  kQ = base_ring(I)

  J, _ = irreducible_hull(quotient_ring_as_module(I))
  return [_get_irreducible_ideal(kQ, I) for I in J.indec_injectives]
end

@doc raw"""
  _get_irreducible_ideal(kQ::MonoidAlgebra, J::IndecInj)

Given a monoid algebra $k[Q]$ and an indecomposable injective $J = k\{a + F - Q\}$ return the irreducible ideal $W\subseteq k[Q]$ such that $J_Q = k[Q]/W$ ($Q$-graded part of $J$). 
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
    irreducible_res(M::SubquoModule{<:MonoidAlgebraElem}, i::Int = 0)

Let $k[Q]$ be a monoid algebra and let $M$ be a $\mathbb{Z}^d$-graded module over $k[Q]$. This function computes a minimal irreducible resolution
\[0 \to M \xrightarrow{ϵ} \overline{W}^0 \xrightarrow{d^0} \overline{W}^1 \dots \xrightarrow{d^{r-1} \overline{W}^r,\]
where $\overline{W}^i = \sum_{j=1}^{n_i} \overline{W_{i_j}} = \sum_{j=1}^{n_i} k[Q]/W_{i_j}$ for irreducible ideals $W_{i_j}$. 

# Examples 
```jldoctest
julia> R_Q,(x,y) = graded_polynomial_ring(QQ,["x","y"];weights = [[1,0],[0,1]])
(Graded multivariate polynomial ring in 2 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y])

julia> kQ = monoid_algebra(R_Q)
monoid algebra over rational field with cone of dimension 2


julia> I = ideal(kQ,[x^4,x^2*y^2,y^4])
ideal over monoid algebra over rational field with cone of dimension 2
 generated by x^4, x^2*y^2, y^4


julia> M = quotient_ring_as_module(I)
Graded subquotient of graded submodule of R_Q^1 with 1 generator
  1: e[1]
by graded submodule of R_Q^1 with 3 generators
  1: x^4*e[1]
  2: x^2*y^2*e[1]
  3: y^4*e[1], 
over monoid algebra generated by x, y quotient by ideal ().

julia> irr_res = irreducible_res(M)
irreducible resolution 
  W^0 -> W^1
where 
 W^0 = direct sum of
    k{[1, 3] + F - Q}_Q, where p_F = Ideal (y, x, x*y)
    k{[3, 1] + F - Q}_Q, where p_F = Ideal (y, x, x*y)
 W^1 = direct sum of
    k{[1, 1] + F - Q}_Q, where p_F = Ideal (y, x, x*y)
with maps
 W^{i-1} -> W^i for 1 <= i <= 1
of Graded subquotient of graded submodule of R_Q^1 with 1 generator
  1: e[1]
by graded submodule of R_Q^1 with 3 generators
  1: x^4*e[1]
  2: x^2*y^2*e[1]
  3: y^4*e[1]
over monoid algebra generated by x, y quotient by ideal ()
```
"""
function irreducible_res(M::SubquoModule{<:MonoidAlgebraElem}, i::Int=0)
  kQ = base_ring(M)
  R_Q = kQ.algebra
  Mi = M # current module in resolution

  #initialize
  gi = identity_map(Mi)
  res_Wi = Vector{InjMod}()
  res_Mi = [Mi] #cokernels 
  res_hi = Vector{SubQuoHom}()
  res_fi = Vector{SubQuoHom}()
  res_gi = [gi] #quotient maps

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

    #fix modules with "zero" relations
    if length(filter(is_zero, relations(Mi))) > 0
      Mi, h = fix_module(Mi)
      gi = gi*h
    end

    push!(res_Mi, Mi)
    push!(res_gi, gi)
    push!(res_Wi, Ji)
    # TODO: Remove once bugs are gone!
    if !isempty(res_hi) && !is_zero(compose(last(res_hi), hi))
      @show length(res_hi)
      @show last(res_hi)
      @show hi
      error()
    end
    push!(res_hi, hi)
    push!(res_fi, fi)

    # end at cohomological degree i
    if i > 0 && j == i
      break
    end
    j = j + 1
  end

  #get cochain complex
  C = cochain_complex(res_hi)

  return IrrRes(M, res_Wi, res_hi, res_gi, res_fi, res_Mi, C)
end

@doc raw"""
    injective_res(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)

Let $k[Q]$ be a monoid algebra and $M$ a finitely generated $Q$-graded module over $k[Q]$. This function computes a minimal injective
resolution 
\[0 \to M \xrightarrow{ϵ} J^0 \xrightarrow{d^0} J^1 \xrightarrow{d^1} \dots \xrightarrow{d^{i-1}} J^i.\]
The maps $d^j$ are given by monomial matrices. This is an implementation of the algorithms in~\cite{HM05}.

# Examples
```jldoctest
julia> R_Q,(x,y) = graded_polynomial_ring(QQ,["x","y"];weights = [[1,0],[0,1]])
(Graded multivariate polynomial ring in 2 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y])

julia> kQ = monoid_algebra(R_Q)
monoid algebra over rational field with cone of dimension 2


julia> I = ideal(kQ,[x^4,x^2*y^2,y^4])
ideal over monoid algebra over rational field with cone of dimension 2
 generated by x^4, x^2*y^2, y^4


julia> M = quotient_ring_as_module(I)
Graded subquotient of graded submodule of R_Q^1 with 1 generator
  1: e[1]
by graded submodule of R_Q^1 with 3 generators
  1: x^4*e[1]
  2: x^2*y^2*e[1]
  3: y^4*e[1], 
over monoid algebra generated by x, y quotient by ideal ().

julia> injective_res(M,2)
injective resolution 
  J^0 -> J^1 -> J^2
where 
 J^0 = direct sum of
    k{[1, 3] + F - Q}, where p_F = Ideal (y, x, x*y)
    k{[3, 1] + F - Q}, where p_F = Ideal (y, x, x*y)
 J^1 = direct sum of
    k{[-1, 3] + F - Q}, where p_F = Ideal (y, x, x*y)
    k{[1, 1] + F - Q}, where p_F = Ideal (y, x, x*y)
    k{[3, -1] + F - Q}, where p_F = Ideal (y, x, x*y)
 J^2 = direct sum of
    k{[-1, -1] + F - Q}, where p_F = Ideal (y, x, x*y)
with maps
 J^{j-1} -> J^j for 1 <= j <= 2
of Graded subquotient of graded submodule of R_Q^1 with 1 generator
  1: e[1]
by graded submodule of R_Q^1 with 3 generators
  1: x^4*e[1]
  2: x^2*y^2*e[1]
  3: y^4*e[1]
over monoid algebra generated by x, y quotient by ideal ()
```

```jldoctest
julia> S,(x,y,z) = graded_polynomial_ring(QQ,["x","y","z"]; weights = [[0,1],[1,1],[2,1]])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> J = ideal(S,[x*z-y^2])
Ideal generated by
  x*z - y^2

julia> kQ = monoid_algebra(quo(S,J)[1])
monoid algebra over rational field with cone of dimension 2


julia> R = kQ.algebra
Quotient
  of multivariate polynomial ring in 3 variables over QQ graded by
    x -> [0 1]
    y -> [1 1]
    z -> [2 1]
  by ideal (x*z - y^2)

julia> F = graded_free_module(R,2)
Graded free module R^2([0 0]) of rank 2 over R

julia> a = R[y y;0 x^2]
[y     y]
[0   x^2]

julia> b = R[x^2*z 0; x^4*y 0; 0 x^5*y; 0 z^3]
[x^2*z       0]
[x^4*y       0]
[    0   x^5*y]
[    0     z^3]

julia> M = monoid_algebra_module(kQ,SubquoModule(F,a,b))
Graded subquotient of graded submodule of F with 2 generators
  1: y*e[1] + y*e[2]
  2: x^2*e[2]
by graded submodule of F with 4 generators
  1: x^2*z*e[1]
  2: x^4*y*e[1]
  3: x^5*y*e[2]
  4: z^3*e[2], 
over monoid algebra generated by x, y, z quotient by ideal (x*z - y^2).

julia> injective_res(M,2)
injective resolution 
  J^0 -> J^1 -> J^2
where 
 J^0 = direct sum of
    k{[1, 4] + F - Q}, where p_F = Ideal (x, y, z, x*z)
    k{[5, 7] + F - Q}, where p_F = Ideal (x, y, z, x*z)
    k{[4, 7] + F - Q}, where p_F = Ideal (x, y, z, x*z)
    k{[0, 2] + F - Q}, where p_F = Ideal (y, z, x*z)
    k{[1, 2] + F - Q}, where p_F = Ideal (x, y, x*z)
 J^1 = direct sum of
    k{[1, 2] + F - Q}, where p_F = Ideal (x, y, z, x*z)
    k{[0, 3] + F - Q}, where p_F = Ideal (x, y, z, x*z)
    k{[-1, 3] + F - Q}, where p_F = Ideal (x, y, z, x*z)
    k{[5, 4] + F - Q}, where p_F = Ideal (x, y, z, x*z)
    k{[0, 5] + F - Q}, where p_F = Ideal (x, y, z, x*z)
    k{[4, 6] + F - Q}, where p_F = Ideal (x, y, z, x*z)
    k{[3, 6] + F - Q}, where p_F = Ideal (x, y, z, x*z)
    k{[-1, -2] + F - Q}, where p_F = Ideal (y, z, x*z)
    k{[-4, -2] + F - Q}, where p_F = Ideal (x, y, x*z)
 J^2 = direct sum of
    k{[0, 0] + F - Q}, where p_F = Ideal (x, y, z, x*z)
    k{[-1, 1] + F - Q}, where p_F = Ideal (x, y, z, x*z)
    k{[-1, 2] + F - Q}, where p_F = Ideal (x, y, z, x*z)
    k{[-2, 2] + F - Q}, where p_F = Ideal (x, y, z, x*z)
    k{[3, 5] + F - Q}, where p_F = Ideal (x, y, z, x*z)
    k{[2, 5] + F - Q}, where p_F = Ideal (x, y, z, x*z)
with maps
 J^{j-1} -> J^j for 1 <= j <= 2
of Graded subquotient of graded submodule of F with 2 generators
  1: y*e[1] + y*e[2]
  2: x^2*e[2]
by graded submodule of F with 4 generators
  1: x^2*z*e[1]
  2: x^4*y*e[1]
  3: x^5*y*e[2]
  4: z^3*e[2]
over monoid algebra generated by x, y, z quotient by ideal (x*z - y^2)
```
"""
function injective_res(M::SubquoModule{<:MonoidAlgebraElem}, i::Int)
  kQ = base_ring(M)
  G = grading_group(kQ)

  #compute irreducible resolution of shifted module
  a_shift = compute_shift(M, i+1)
  M_a = twist(M, -G(a_shift))
  irr_res = irreducible_res(M_a)

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
    injective_res(I::MonoidAlgebraIdeal,i::Int)
Let $k[Q]$ be a monoid algebra and $I\subset k[Q]$ a $\mathbb{Z}^d$-graded ideal. This function computes an injective resolution of $M = k[Q]/I$.  
"""
function injective_res(I::MonoidAlgebraIdeal, i::Int)
  return injective_res(quotient_ring_as_module(I), i)
end

##fix SubquoModule with "zero"-relations
function fix_module(M::SubquoModule)
  R = base_ring(M)
  M_F = ambient_free_module(M)
  M_gens = gens(M)
  if length(M_gens) > 0
    M_rels = filter(!is_zero, relations(M))
    d = length(M_gens)
    M_fixed, _ = quo(sub(M_F, M_gens)[1], M_rels)
    return M_fixed, hom(M, M_fixed, identity_matrix(R, d))
  else
    return M, identity_map(M)
  end
end


#get Krull dimension of module
# TODO: For modules over polynomial rings this should make 
# use of the `dim` in Singular. But this does not seem 
# to be available as of yet.
@attr Int function dim(M::ModuleFP)
  ann = annihilator(M)
  return dim(ann)
end

#=
function dim(F::FreeMod)
  return dim(base_ring(F))
end
=#

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


@attr Bool function is_normal(A::MonoidAlgebra{<:FieldElem, <:MPolyQuoRing})
  # Implementation adapted from 
  #
  #    M2/Macaulay2/packages/IntegralClosure.m2,
  #
  # line 666 ff. of https://github.com/Macaulay2/M2/blob/2565455411d15a3386204aa62a00e20ee5c0e99f/M2/Macaulay2/packages/IntegralClosure.m2 
  # on Apr 25, 2025.
  R = A.algebra::MPolyQuoRing
  I = modulus(R)
  M = quotient_ring_as_module(I)
  n = codim(I)                # Calculate the codimension of the ideal

  # Check the S2 condition
  test_range = 0:(dim(R) - n - 2)

  all(dim(R) - dim(ext(M, graded_free_module(R_B, 1), j + n + 1)) >= (j + n + 3) for j in test_range) || return false # S2 condition not satisfied

  Jac = ideal(R, minors(map_entries(R, jacobian_matrix(gens(modulus(R)))), n))  # Compute minors of the Jacobian
  d = dim(Jac)                   # Get dimension of the Jacobian
  d < 0 && (d = -Inf)            # Handle negative dimensions
  return (dim(R) - d >= 2)       # Check the condition
end


# import local cohomology functions
include("LocalCohomology.jl")
include("ModuleFunctionality.jl")

end # module InjectiveResolutions

using .InjectiveResolutions

## Functions visible on the outside
#export get_monoid_algebra
export monoid_algebra_from_lattice
export monoid_algebra
#export monoid_algebra_ideal
#export monoid_algebra_module

export irreducible_res
export irreducible_dec

export injective_res

#export mod_quotient
#export mod_saturate
#export ZF_basis # should be renamed ?
export InjMod
export IndecInj
#export underlying_element
#export underlying_ideal
export irreducible_hull
#export _get_irreducible_ideal
#export compute_shift
#export Q_graded_part # should be renamed?
#export compute_Q_graded_part

export MonoidAlgebra
export MonoidAlgebraElem
