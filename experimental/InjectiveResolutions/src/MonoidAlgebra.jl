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
  A::Union{Matrix{Int},Nothing} #semigroup generators of the face as columns of a matrix
end

struct HyperplaneQ # a hyperplane bounding the cone RR_{\geq 0}Q
  hyperplane::Polyhedron
  A::Matrix{Int}
  b::Vector{Int}
end
@attributes mutable struct AffineSemigroup # affine semigroup
  generators::Matrix{Int}  # semigroup generators as columns of the matrix

  function AffineSemigroup(generators::Matrix{Int})
    return new(generators)
  end
end

function Base.show(io::IO, S::AffineSemigroup)
  print(io, "affine semigroup in ambient dimension ", ambient_dimension(S))
end

function Base.show(io::IO, ::MIME"text/plain", S::AffineSemigroup)
  print(io, "affine semigroup with generators: ")
  gens_list = gens(S)
  for (i, s) in enumerate(gens_list)
    print(io, s)
    if i < length(gens_list)
      print(io, ", ")
    end
  end
end

function gens(Q::AffineSemigroup)
  return [Q.generators[:, i] for i in 1:size(Q.generators, 2)]
end

function rank(Q::AffineSemigroup)
  return rank(matrix(ZZ,Q.generators))
end

function ambient_dimension(Q::AffineSemigroup)
  return size(Q.generators, 1)
end

@doc raw"""
    affine_semigroup(M::Matrix{Int})

Create an affine semigroup from a matrix of integer generators (as columns).
"""
function affine_semigroup(M::Matrix{Int})
  return AffineSemigroup(M)
end

@doc raw"""
    affine_semigroup(V_Q::Vector{Vector{Int}})

Create an affine semigroup.
"""
function affine_semigroup(V_Q::Vector{Vector{Int}})
  return affine_semigroup(hcat(V_Q...))
end

@doc raw"""
    is_pointed(S::AffineSemigroup)

Given an affine semigroup, check if its cone is pointed.
"""
@attr Bool function is_pointed(S::AffineSemigroup)
  return is_pointed(polyhedral_cone(S))
end

@attr Cone{QQFieldElem} function polyhedral_cone(S::AffineSemigroup)
  G = [S.generators[:, i] for i in 1:size(S.generators, 2)]
  return positive_hull(G)
end

@attr Polyhedron function cone(S::AffineSemigroup)
  return polyhedron(polyhedral_cone(S))
end

@attr Vector{HyperplaneQ} function hyperplanes(S::AffineSemigroup)
  return get_bounding_hyperplanes(cone(S))
end

@attr Tuple{Polyhedron,Vector{Int}} function zonotope(S::AffineSemigroup)
  return get_zonotope(cone(S))
end

@attributes mutable struct MonoidAlgebra{CoeffType,AlgebraType} <: Ring # monoid algebra with associated data
  algebra::AlgebraType
  affine_semigroup::AffineSemigroup

  function MonoidAlgebra(
    A::AlgebraType, Q::AffineSemigroup; check::Bool=true
  ) where {AlgebraType<:Union{MPolyRing,MPolyQuoRing}}
    @check is_zm_graded(A) "given algebra is not ZZ^d-graded"
    gg_Q = grading_group(A)
    @check is_free(gg_Q) && is_abelian(gg_Q) "given algebra is not a monoid algebra"
    kk = coefficient_ring(A)
    return new{elem_type(kk),AlgebraType}(A, Q)
  end

  function MonoidAlgebra(
    A::AlgebraType; check::Bool=true
  ) where {AlgebraType<:Union{MPolyRing,MPolyQuoRing}}
    @check is_zm_graded(A) "given algebra is not ZZ^d-graded"
    gg_Q = grading_group(A)
    @check is_free(gg_Q) && is_abelian(gg_Q) "given algebra is not a monoid algebra"
    D = [degree(Vector{Int}, g) for g in gens(A)]
    Q = AffineSemigroup(hcat(D...))
    kk = coefficient_ring(A)
    return new{elem_type(kk),AlgebraType}(A, Q)
  end
end

is_graded(A::MonoidAlgebra) = true
is_zm_graded(A::MonoidAlgebra) = true

@doc raw"""
    affine_semigroup(A::MonoidAlgebra)

Given a monoid algebra kQ, this function returns the underlying affine semigroup Q.
"""
function affine_semigroup(A::MonoidAlgebra)
  return A.affine_semigroup
end

@doc raw"""
    semigroup_generators(A::MonoidAlgebra)

Given a monoid algebra kQ, this function returns the matrix of semigroup generators as columns.
"""
function semigroup_generators(A::MonoidAlgebra)
  return affine_semigroup(A).generators
end

function polyhedral_cone(A::MonoidAlgebra)
  return polyhedral_cone(affine_semigroup(A))
end

@doc raw"""
    cone(A::MonoidAlgebra)

Given a monoid algebra with underlying monoid $Q$, this function returns the polyhedral cone $\mathbb{R}_{\geq 0}Q$ as a polyhedron. 

"""
function cone(A::MonoidAlgebra)
  return cone(affine_semigroup(A))
end


@doc raw"""
    faces(A::MonoidAlgebra)

Given a monoid algebra with underlying monoid $Q$, this function a list of all faces of the polyhedral cone $\mathbb{R}_{\geq 0}Q$ with their corresponding homogeneous prime ideals. 
"""
@attr Vector{FaceQ} function faces(A::MonoidAlgebra)
  # Compute all faces of the cone
  P = cone(A)
  P_faces = Vector{Polyhedron}()
  for i in 0:dim(P)
    append!(P_faces, Oscar.faces(P, i))
  end

  # Get semigroup generators
  D = semigroup_generators(A)

  # For each face, compute its semigroup generators and prime ideal
  _faces = Vector{FaceQ}()
  for F in P_faces
    # Determine which generators lie on this face
    A_face = [D[:,i] for i in 1:size(D,2) if is_subset(convex_hull(D[:,i]), F)]
    if is_empty(A_face)
      face_gens = nothing
    else
      face_gens = hcat(A_face...)
    end

    # Compute the corresponding homogeneous prime ideal
    prime = prime_of_face(A.algebra, F)
    push!(_faces, FaceQ(prime, F, face_gens))
  end

  return _faces
end

function hyperplanes(A::MonoidAlgebra)
  return hyperplanes(affine_semigroup(A))
end

function is_pointed(A::MonoidAlgebra)
  semigroup = affine_semigroup(A)
  return is_pointed(semigroup)
end

function zonotope(A::MonoidAlgebra)
  return zonotope(affine_semigroup(A))
end

coefficient_ring(A::MonoidAlgebra) = coefficient_ring(A.algebra)
number_of_variables(A::MonoidAlgebra) = ngens(A.algebra)

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
  test_range = 0:(krull_dim(R_B) - n - 2)

  for j in test_range
    # Check if codimension of Ext^{j+n+1} is at least j+n+3
    E = ext(M, graded_free_module(R_B, 1), j + n + 1)
    # work around issue https://github.com/oscar-system/Oscar.jl/issues/4884
    if is_zero(E)
      d = -1
    else
      d = krull_dim(E)
    end
    cod = krull_dim(R_B) - d
    if cod < j + n + 3
      return false
    end
  end

  Jac = ideal(R, minors(map_entries(R, jacobian_matrix(gens(modulus(R)))), n))
  d = krull_dim(Jac)
  d < 0 && return true
  return krull_dim(R) - d >= 2
end

is_normal(A::MonoidAlgebra{<:FieldElem, <:MPolyRing}) = true

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

function degree(
    ::Type{Vector{Int}},
    a::MonoidAlgebraElem{<:FieldElem, PT};
    check::Bool=true
  ) where {RT <: MPolyRing, PT <: MonoidAlgebra{<:FieldElem, RT}}
  !isdefined(a, :elem) && return zero(grading_group(parent(a)))
  return degree(Vector{Int},underlying_element(a))
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

krull_dim(A::MonoidAlgebra) = krull_dim(A.algebra)

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
    prime_of_face(kQ::Union{MPolyRing,MPolyQuoRing}, F::Polyhedron)

Let kQ be a monoid algebra over some semigroup $Q$. Given a face $F$ of the cone $C = \RR_{\geq 0}Q$ of the monoid algebra kQ,
return the corresponding homogeneous prime ideal

$P_F = k\{Q\setminus F\}.$ 
"""
function prime_of_face(kQ::Union{MPolyRing,MPolyQuoRing},F::Polyhedron)
  #get degrees of generators of kQ
  G = [degree(Vector{Int},g) for g in Oscar.gens(kQ)]

  #generators of P_F
  gens_PF = [g for g in G if !is_subset(convex_hull(g),F)]
  return ideal(kQ,[monomial_basis(kQ,g)[1] for g in gens_PF])
end

# given a monoid algebra, this function returns the corresponding polyhedral cone
# INPUT:    monoid algebra k[Q]
# OUTPUT:   polyhedral cone \RR_{\geq 0}Q
function get_polyhedral_cone(R::Union{MPolyDecRing,MPolyQuoRing})
  D = [degree(Vector{Int}, g) for g in gens(R)]
  return positive_hull(D)
end


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
  b = [M[:, 1]; -M[:, 1]]
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

  t_vars = [Symbol("t_$i") for i in 1:d]
  if all(M_Q .>= 0)
    T, t = graded_polynomial_ring(k, t_vars; cached = false)
  else
    T, t = laurent_polynomial_ring(k, t_vars) 
  end
  # construct k[x_1,...,x_n] where n is the number of columns/generators
  x_vars = [Symbol("x_$i") for i in 1:size(M_Q, 2)]
  R, _ = graded_polynomial_ring(
    k, x_vars; weights=[Vector(row) for row in eachcol(M_Q)], cached = false
  )

  # construct map x_i \to t^(a_i)
  targ = [prod(t[j]^M_Q[j, i] for j in 1:d) for i in 1:size(M_Q, 2)]
  map_T_R = hom(R, T, targ)

  # return monoid algebra
  Q = AffineSemigroup(M_Q)
  if is_zero(ideal(gens(kernel(map_T_R))))
    kQ = MonoidAlgebra(R, Q)
  else
    kQ = MonoidAlgebra(quo(R, ideal(gens(kernel(map_T_R))))[1], Q)
  end
  return kQ
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

@doc raw"""
    monoid_algebra(Q::AffineSemigroup, k::Field)

Return the monoid algebra over affine semigroup.
"""
function monoid_algebra(Q::AffineSemigroup, k::Field)
  @assert is_pointed(Q) "the semigroup must be pointed"
  return monoid_algebra(Q.generators,k)
end

# compute the saturation of a monoid algebra
@attr MonoidAlgebra function saturation(kQ::MonoidAlgebra)
  C = polyhedral_cone(kQ)
  k = coefficient_ring(kQ)
  Csat = matrix(ZZ, hilbert_basis(C))
  Csat_int = Int.(Csat)
  Csat_gens = [Csat_int[i,:] for i in 1:size(Csat_int,1)]
  return monoid_algebra(Csat_gens, k)
end

@doc raw"""
    saturation_map(kQ::MonoidAlgebra)

Return the inclusion map $k[Q] \hookrightarrow k[\overline{Q}]$ from a monoid algebra to its saturation.
"""
function saturation_map(kQ::MonoidAlgebra)
  kQsat = saturation(kQ)
  im_phi = [monomial_basis(kQsat, degree(g))[1] for g in gens(kQ.algebra)]
  return Oscar.hom(kQ, kQsat, im_phi)
end

@doc raw"""
    saturation_ideal(kQ::MonoidAlgebra)

Return the ideal in $k[\overline{Q}]$ generated by the image of the generators of $k[Q]$,
i.e., the ideal defining $k[Q]$ as a $k[\overline{Q}]$-module.
"""
function saturation_ideal(kQ::MonoidAlgebra)
  kQsat = saturation(kQ)
  phi = saturation_map(kQ)
  return ideal(kQsat, [phi(underlying_element(g)) for g in gens(kQ)])
end

@doc raw"""
    holes_module(kQ::MonoidAlgebra)

Return the $k[Q]$-module $k[\overline{Q}]/k[Q]$, where $\overline{Q}$ is the saturation of $Q$.
This module encodes the "holes" of the semigroup, i.e., the elements in $\overline{Q} \setminus Q$.

Returns a `SubquoModule` over `kQ`.
"""
function holes_module(kQ::MonoidAlgebra)
  if is_normal(kQ)
    F = graded_free_module(kQ, 0)
    return quo(F, [zero(F)])[1]
  end

  # present kQsat as a kQ-module: kQ^r →^PM kQ^s → kQsat → 0
  # the last generator in gs is always 1, corresponding to the copy of kQ
  # dropping it directly gives kQsat/kQ
  phi = saturation_map(kQ)
  gs, PM, sect = present_finite_extension_ring(phi)

  s = length(gs)
  r = size(PM, 1)

  # build graded free module over kQ with shifts from the non-trivial generators
  kQsat = saturation(kQ)
  G = grading_group(kQ)
  degrees = [G(degree(Vector{Int}, kQsat(g))) for g in gs[1:s-1]]
  F = graded_free_module(kQ, degrees)

  # relations: first s-1 columns of PM (since e_s = 0)
  rels = [sum(kQ(PM[i,j]) * F[j] for j in 1:s-1) for i in 1:r]
  filter!(!is_zero, rels)

  return quo(F, rels)[1]
end

@doc raw"""
    is_Q_graded(M::SubquoModule{<:MonoidAlgebraElem})

Check if all generators of $M$ have degrees in the semigroup $Q$.
"""
function is_Q_graded(M::SubquoModule{<:MonoidAlgebraElem})
  kQ = base_ring(M)
  return all(g -> is_zero(g) || is_in_semigroup(kQ, degree(Vector{Int}, g)), gens(M))
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

function number_of_generators(I::MonoidAlgebraIdeal{ElemType}) where {ElemType}
  return length(I.gens)
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

function minimal_generating_set(I::MonoidAlgebraIdeal)
  kQ = base_ring(I)
  return [kQ(g) for g in minimal_generating_set(underlying_ideal(I))]
end

dim(I::MonoidAlgebraIdeal) = krull_dim(underlying_ideal(I))
krull_dim(I::MonoidAlgebraIdeal) = krull_dim(underlying_ideal(I))

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
Base.hash(a::MonoidAlgebraElem, h::UInt) = hash(underlying_element(a), h)
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

# Let p in \ZZ^d, a in \ZZ^d and F a face of a semigroup Q. This function checks if p is in a + ZF.
function is_in_aZF(a::Vector{Int},F::FaceQ,p::Vector{Int})
    if F.A === nothing
        return a == p
    else
        A = vcat(F.A,-F.A)
        #we want to check if p is in a + ZF, where F is a face of the semigroup generated by the columns of A
        #this is equivalent to checking if p is in the semigroup generated by the columns of A and -A
        b = vcat(p .- a,-p .+ a)
        P = polyhedron(A,b)
        milp = mixed_integer_linear_program(P,zeros(Int64,size(A,2)))
        return solve_milp(milp) != (nothing,nothing)
    end
end

# return integer identity matrix
function I_n(n)
    I = zeros(Int64, n, n)
    for i in 1:n
        I[i, i] = 1
    end
    return I
end

@attr Dict{Vector{Int},Bool} function _is_in_semigroup_cache(kQ::MonoidAlgebra)
    return Dict{Vector{Int},Bool}()
end

#check if a point in ZZ^d is in the semigroup generated by the degrees of the generators of a monoid algebra
function is_in_semigroup(kQ::MonoidAlgebra, p::Vector{Int})
    cache = _is_in_semigroup_cache(kQ)
    return get!(cache, p) do
        D = [degree(Vector{Int},g) for g in gens(kQ)]
        _is_in_semigroup(hcat(D...),p)
    end
end

#check if a point in ZZ^d is in the semigroup generated by the columns of a matrix A
function _is_in_semigroup(A,p)
    @assert size(A,1) == length(p)
    n = size(A,2)
    #only positive integer combinations
    _In = - I_n(n)
    _A = vcat(_In,A,-A)
    _b = vcat(zeros(Int64,n),p,-p)
    P = polyhedron(_A,_b)
    milp = mixed_integer_linear_program(P,zeros(Int64,n))
    return solve_milp(milp) != (nothing,nothing)
end

# this function checks if the semigroup generated Z^d
function generates_Zd(kQ::MonoidAlgebra)
  _M = [degree(Vector{Int},g) for g in gens(kQ)]
  M = matrix(ZZ,hcat(_M...))
  S = snf(M)
  diag_entries = [S[i,i] for i in 1:min(nrows(S), ncols(S))]
  return rank(M) == rank(grading_group(kQ)) && all(diag_entries .== 1)
end

# check if the intersection (g + Q) \cap (a + ZF) is non-empty
function in_intersection(kQ::MonoidAlgebra,g_vec::Vector{Int},a_vec::Vector{Int},F::FaceQ)
    # cheap cone-level pre-filter: (g + cone(Q)) ∩ (a + ZF) = ∅ implies (g + Q) ∩ (a + ZF) = ∅
    g_p = convex_hull(g_vec)
    a_p = convex_hull(a_vec)
    if dim(intersect(g_p + cone(kQ), a_p + F.poly + (-1)*F.poly)) < 0
        return false
    end

    n = ngens(kQ)
    _D = [degree(Vector{Int},g) for g in gens(kQ)]
    D = hcat(_D...)
    if F.A === nothing
        _A = vcat(D,-D)
        A = vcat(_A,-I_n(n))
    else
        _A = vcat(hcat(-F.A,D),hcat(F.A,-D))
        A = vcat(_A,hcat(zeros(Int64,n,size(F.A,2)),-I_n(n)))
    end
    b = vcat(a_vec.-g_vec,g_vec .- a_vec, zeros(Int64,n))
    P = polyhedron(A,b)
    milp = mixed_integer_linear_program(P,zeros(Int64,size(A,2)))
    return solve_milp(milp) != (nothing,nothing)
end
