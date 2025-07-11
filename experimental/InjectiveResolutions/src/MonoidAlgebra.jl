########################################################################
# MonoidAlgebra.jl
#
# This is a wrapper type for `MPolyRing` or `MPolyQuoRing` which is 
# graded and available at the field `.algebra`. This datum is enhanced 
# by various combinatorial stuff. 
#
# In order to facilitate functionality for ideals, we introduce a 
# wrapper type `MonoidAlgebraIdeal` for `MPolyIdeal` (resp. 
# `MPolyQuoIdeal`). For finitely generated modules, on the other hand, 
# we implement the translation to the singular side within the 
# `ModuleGens` directly; see `ModuleFunctionality.jl` for that. 
########################################################################

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

@doc raw"""
    cone(A::MonoidAlgebra)

Given a monoid algebra with underlying monoid $Q$, this function returns the polyhedral cone $\mathbb{R}_{\geq 0}Q$ as a polyhedron. 

"""
function cone(A::MonoidAlgebra)
  if !isdefined(A, :cone)
    A.cone = polyhedron(polyhedral_cone(A))
  end
  return A.cone
end


@doc raw"""
    faces(A::MonoidAlgebra)

Given a monoid algebra with underlying monoid $Q$, this function a list of all faces of the polyhedral cone $\mathbb{R}_{\geq 0}Q$ with their corresponding homogeneous prime ideals. 
"""
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

is_unit(a::MonoidAlgebraElem) = is_unit(underlying_element(a))

function inv(a::MonoidAlgebraElem)
  return MonoidAlgebraElem(parent(a), inv(underlying_element(a)))
end

monomial_basis(A::MonoidAlgebra, g::FinGenAbGroupElem) = monomial_basis(A.algebra, g)

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
  return krull_dim(A.algebra)
end


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

@doc raw"""
    prime_to_face(kQ::Union{MPolyRing,MPolyQuoRing}, zonotope::Polyhedron, F::Polyhedron)

Let kQ be a monoid algebra over some semigroup $Q$. Given a face $F$ of the cone $C = \RR_{\geq 0}Q$ of the monoid algebra kQ,
return the corresponding homogeneous prime ideal

$P_F = k\{Q\setminus F\}.$

The zonotope is the Minkowski sum of all primitive integer vectors along rays of $Q$. 
"""
function prime_to_face(
  kQ::Union{MPolyRing,MPolyQuoRing}, zonotope::Polyhedron, F::Polyhedron
)
  gens = Vector{MPolyDecRingElem}()
  for lp in lattice_points(zonotope)
    if !(lp in F) #check if lattice point is in F
      a_v = [a_p for a_p in lp]
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
    monoid_algebra(B::Matrix{Int},k::Field)

Return the monoid algebra generated by monomials $x^{v_1},\dots,x^{v_n}\in k[x_1,\dots,x_d]$, where $v_1,\dots,v_n\in \mathbb{Z}^d$ are the columns of $M_Q$. 

# Examples
```jldoctest
julia> M_Q = [1 0; 0 1]
2×2 Matrix{Int64}:
 1  0
 0  1

julia> monoid_algebra(M_Q,QQ)
monoid algebra over rational field with cone of dimension 2
```
"""
function monoid_algebra(M_Q::Matrix{Int}, k::Field)
  d = size(M_Q, 1)

  # construct k[t_1,...,t_d]
  t_vars = [Symbol("t_$i") for i in 1:d]
  T, t = graded_polynomial_ring(k, t_vars; cached = false)

  # construct k[x_1,...,x_n] where n is the number of columns/generators
  x_vars = [Symbol("x_$i") for i in 1:size(M_Q, 2)]
  R, _ = graded_polynomial_ring(
    k, x_vars; weights=[Vector(row) for row in eachcol(M_Q)], cached = false
  )

  # construct map x_i \to t^(a_i)
  targ = [prod(t[j]^M_Q[j, i] for j in 1:d) for i in 1:size(M_Q, 2)]
  map_T_R = hom(R, T, targ)

  # return monoid algebra
   if is_zero(ideal(gens(kernel(map_T_R))))
    return MonoidAlgebra(R)
   else
    return  MonoidAlgebra(quo(R, ideal(gens(kernel(map_T_R))))[1])
   end
end

@doc raw"""
    monoid_algebra(V_Q::Vector{Vector{Int}},k::Field)

Return the monoid algebra generated by monomials $x^{v_1},\dots,x^{v_n}\in k[x_1,\dots,x_n]$, where `V_Q` $= v_1,\dots,v_n \in \mathbb{Z}^d$.

# Examples
```jldoctest
julia> kQ = monoid_algebra([[0,1],[1,1],[2,1]],QQ)
monoid algebra over rational field with cone of dimension 2

julia> kQ.algebra
Quotient
  of multivariate polynomial ring in 3 variables over QQ graded by
    x_1 -> [0 1]
    x_2 -> [1 1]
    x_3 -> [2 1]
  by ideal (-x_1*x_3 + x_2^2)
```
"""
function monoid_algebra(V_Q::Vector{Vector{Int}}, k::Field)
  return monoid_algebra(Matrix{Int}(transpose(matrix(V_Q))), k)
end

function Base.show(io::IO,F::FaceQ)
  print(io,"face corresponding to homogeneous prime ",lowercase(string(F.prime)))
end

function Base.show(io::IO, ::MIME"text/plain", kQ::MonoidAlgebra)
  print(io, "monoid algebra over ", lowercase(string(coefficient_ring(kQ))),
    " with cone of dimension $(dim(cone(kQ)))",
  )
end

function Base.show(io::IO, kQ::MonoidAlgebra)
  print(io, "monoid algebra over ", lowercase(string(coefficient_ring(kQ))),
    " with cone of dimension $(dim(cone(kQ)))",
  )
end

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

dim(I::MonoidAlgebraIdeal) = krull_dim(underlying_ideal(I))

# some generic functionality which should probably be elsewhere
function is_subset(I::T, J::T) where {T<:Ideal}
  return all(x in J for x in gens(I))
end

function Base.:(==)(I::T, J::T) where {T <: Ideal}
  return is_subset(I, J) && is_subset(J, I)
end

function Base.:*(I::T, J::T) where {T<:Ideal}
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
  @req all(base_ring(g) === kQ for g in b) "base rings must match"
  as = Singular.intersection(singular_generators(a), [singular_generators(g) for g in b]...)
  return MonoidAlgebraIdeal(kQ, [kQ(x) for x in gens(as)])
end

function Base.intersect(V::Vector{T}) where {T <: MonoidAlgebraIdeal}
  @assert length(V) != 0
  length(V) == 1 && return V[1]

  return intersect(V[1], V[2:end]...)
end

# random elements for testing 
rand(A::MonoidAlgebra, v...) = A(rand(A.algebra, v...))

### Additional functionality for ring conformance tests
Base.hash(a::MonoidAlgebraElem, h::UInt) = hash(underlying_element, h)
characteristic(A::MonoidAlgebra) = characteristic(A.algebra)
divexact(a::MonoidAlgebraElem, b::MonoidAlgebraElem; check::Bool=true) = parent(a)(divexact(underlying_element(a), underlying_element(b); check))

function divides(a::MonoidAlgebraElem, b::MonoidAlgebraElem)
  success, q = divides(underlying_element(a), underlying_element(b))
  return success, parent(a)(q)
end

is_nilpotent(a::MonoidAlgebraElem) = is_nilpotent(underlying_element(a))
canonical_unit(a::MonoidAlgebraElem) = canonical_unit(underlying_element(a))
is_domain_type(::Type{MonoidAlgebraElem{CT, PT}}) where {CT, AT, PT <: MonoidAlgebra{CT, AT}} = is_domain_type(elem_type(AT))
# dummy method required by the conformance test
#divrem(a::MonoidAlgebraElem, b::MonoidAlgebraElem; check::Bool=true) = zero(a), a

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

function Base.show(io::IO, ::MIME"text/plain", I::MonoidAlgebraIdeal)
  print(
    io, "ideal over $(base_ring(I)) generated by "*join(["$(x)" for x in gens(I)], ", ")
  )
end

function Base.show(io::IO,I::MonoidAlgebraIdeal)
  print(io,underlying_ideal(I))
end

