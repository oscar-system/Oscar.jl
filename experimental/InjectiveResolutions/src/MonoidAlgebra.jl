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
  return dim(A.algebra)
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

