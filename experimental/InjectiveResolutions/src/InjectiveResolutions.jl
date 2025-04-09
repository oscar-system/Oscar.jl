## Functions visible on the outside
export get_monoid_algebra
export monoid_algebra_from_lattice
export get_monoid_algebra_module
export get_monoid_algebra_ideal

export irreducible_res
export irreducible_dec

export injective_res

export mod_quotient
export mod_saturate

#########################
# some composite types
#########################

struct FaceQ # face of semigroup  
  prime::Union{MPolyIdeal,MPolyQuoIdeal}
  poly::Polyhedron  #face of Q corresponding to prime
end

struct HyperplaneQ # a hyperplane bounding the cone RR_{\geq 0}Q
  hyperplane::Polyhedron
  A::Matrix{Int}
  b::Vector{Int}
end

struct IndecInj # indecomposable injective
  face::FaceQ
  vector::Vector{Int}
end

struct IrrSum # irreducible sum
  mod::SubquoModule
  components::Vector{IndecInj}
end

mutable struct MonoidAlgebra{CoeffType, AlgebraType} <: Ring # monoid algebra with associated data
    algebra::AlgebraType
    pointed::Union{Nothing, Bool} #cone pointed
    polyhedral_cone::Cone{QQFieldElem}
    cone::Polyhedron
    faces::Vector{FaceQ}
    hyperplanes::Vector{HyperplaneQ}
    zonotope::Tuple{Polyhedron,Vector{Int}}

    function MonoidAlgebra(A::AlgebraType; check::Bool=true) where {AlgebraType <: Union{MPolyRing, MPolyQuoRing}}
      #check if monoid algebra
      @check is_zm_graded(A) "given algebra is not ZZ^d-graded"
      gg_Q = grading_group(A)
      @check is_free(gg_Q) && is_abelian(gg_Q) "given algebra is not a monoid algebra"
      kk = coefficient_ring(A)
      result = new{elem_type(kk), AlgebraType}(A, nothing)
    end
end

function polyhedral_cone(A)
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

### Elements of MonoidAlgebras
mutable struct MonoidAlgebraElem{CoeffType, ParentType} <: RingElem
  parent::ParentType
  elem::RingElem
  # TODO: Do we want to store additional information on the elements here?

  function MonoidAlgebraElem(
      A::ParentType
    ) where {CoeffType, ParentType <: MonoidAlgebra{CoeffType}}
    return new{CoeffType, ParentType}(parent)
  end
  
  function MonoidAlgebraElem(
      A::ParentType,
      a::RingElem;
      check::Bool=true
    ) where {CoeffType, ParentType <: MonoidAlgebra{CoeffType}}
    @check parent(a) === A.algebra
    return new{CoeffType, ParentType}(A, a)
  end
end

parent(a::MonoidAlgebraElem) = a.parent

elem_type(::Type{T}) where {CoeffType, T <: MonoidAlgebra{CoeffType}} = MonoidAlgebraElem{CoeffType, T}

function degree(a::MonoidAlgebraElem)
  !isdefined(a, :elem) && return zero(grading_group(parent(a)))
  return degree(underlying_element(a))
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

# TODO: Fix up promote rules?
AbstractAlgebra.promote_rule(::Type{CoeffType}, ::Type{T}) where {CoeffType, T<:MonoidAlgebraElem{CoeffType}} = T

AbstractAlgebra.promote_rule(::Type{Int}, ::Type{T}) where {T<:MonoidAlgebraElem} = T
AbstractAlgebra.promote_rule(::Type{ZZRingElem}, ::Type{T}) where {T<:MonoidAlgebraElem} = T


gens(A::MonoidAlgebra) = [MonoidAlgebraElem(A, x) for x in gens(A.algebra)]
getindex(A::MonoidAlgebra, i::Int) = MonoidAlgebraElem(A, A.algebra[i])

function Base.show(io::IO, a::MonoidAlgebraElem)
  print(io, underlying_element(a))
end

parent_type(::Type{ElemType}) where {ParentType, CoeffType, ElemType <: MonoidAlgebraElem{CoeffType, ParentType}} = ParentType

# TODO: Finish implementation of the ring interface! 

### Ideals over `MonoidAlgebra`s
mutable struct MonoidAlgebraIdeal{ElemType} <: Ideal{ElemType}
  monoidAlgebra::MonoidAlgebra
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
  return I.monoidAlgebra::parent_type(ElemType)
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

# some generic functionality which should probably be elsewhere
function is_subset(I::Ideal, J::Ideal)
  return all(x in J for x in gens(I))
end

# user facing constructor
ideal(A::MonoidAlgebra, v::Vector) = MonoidAlgebraIdeal(A, elem_type(A)[A(x) for x in v])

function Base.in(a::MonoidAlgebraElem, I::MonoidAlgebraIdeal)
  return underlying_element(a) in underlying_ideal(I)
end

function coordinates(a::MonoidAlgebraElem, I::MonoidAlgebraIdeal)
  return coordinates(underlying_element(a), underlying_ideal(I))
end

struct MonoidAlgebraModule
  monoidAlgebra::MonoidAlgebra
  mod::SubquoModule
end

struct InjMod
  monoidAlgebra::MonoidAlgebra
  indecInjectives::Vector{IndecInj}
end
struct IrrRes # irreducible resolution (including all computed data and the cochain complex)
  mod::MonoidAlgebraModule
  irrSums::Vector{IrrSum}
  cochainMaps::Vector{SubQuoHom}
  surjections::Vector{SubQuoHom}
  inclusions::Vector{SubQuoHom}
  cokernels::Vector{SubquoModule}
  cochainComplex::ComplexOfMorphisms # if sequence not exact return trivial cochain_complex (M0 -> M0)
end
struct InjRes
  mod::MonoidAlgebraModule
  injMods::Vector{InjMod}
  cochainMaps::Vector{MatElem}
  upto::Int
  irrRes::IrrRes
  shift::Vector{Int} #not needed
end

function Base.show(io::IO, ::MIME"text/plain", M::MonoidAlgebraModule)
  println(io, M.mod, ", ")
  print(
    io,
    "over monoid algebra generated by ",
    join(gens(M.monoidAlgebra.algebra), ", "),
    " quotient by ideal (",
    join(gens(modulus(M.monoidAlgebra.algebra))),
    ").",
  )
end

function Base.show(io::IO,::MIME"text/plain",kQ::MonoidAlgebra)
    println(io,"monoid algebra over ", lowercase(string(coefficient_ring(kQ))), 
            " with cone of dimension $(dim(cone(kQ)))"
           )
    #=
    println(io,"Monoid algebra k[Q] over ", lowercase(string(coefficient_ring(kQ.algebra))), " generated by ", join(gens(kQ.algebra), ", "), " quotient by ideal (", join(gens(modulus(kQ.algebra)),", "),"), ")
    print(io, "where Q is a subset of ZZ^", ambient_dim(kQ.cone)," and the cone RR_{>= 0}Q has dimension ", dim(kQ.cone), ".")
    =#
end

function Base.show(io::IO, kQ::MonoidAlgebra)
    println(io,"monoid algebra over ", lowercase(string(coefficient_ring(kQ))), 
            " with cone of dimension $(dim(cone(kQ)))"
           )
end

function Base.show(io::IO,::MIME"text/plain",I::MonoidAlgebraIdeal)
  println(io,"ideal over $(base_ring(I)) generated by "*join(["$(x)" for x in gens(I)], ", "))
end

function Base.show(io::IO, ::MIME"text/plain", J::InjMod)
  println(io, "Injective module given by direct sum of indecomposable injectives")
  for Ji in J.indecInjectives
    println(io, "  k{", Ji.vector, " + F - Q}, where p_F = ", Ji.face.prime)
  end
  print(
    io,
    "over monoid algebra generated by ",
    join(gens(J.monoidAlgebra.algebra), ", "),
    " quotient by ideal (",
    join(gens(modulus(J.monoidAlgebra.algebra)), ", "),
    ").",
  )
end

function Base.show(io::IO, ::MIME"text/plain", res::InjRes)
  println(io, "Injective resolution ")
  println(io, "  ", join(["J^$i" for i in 0:res.upto], " -> "))
  println(io, "where ")
  for i in eachindex(res.injMods)
    j = i-1
    println(io, " J^$j = direct sum of")
    for Ji in res.injMods[i].indecInjectives
      println(io, "    k{", Ji.vector, " + F - Q}, where p_F = ", Ji.face.prime)
    end
  end
  println(io, "with maps")
  println(io, " J^{j-1} -> J^j for 1 <= j <= ", res.upto)
  println(io, "of ", res.mod.mod)
  println(
    io,
    "over monoid algebra generated by ",
    join(gens(res.mod.monoidAlgebra.algebra), ", "),
    " quotient by ideal (",
    join(gens(modulus(res.mod.monoidAlgebra.algebra)), ", "),
    ").",
  )
end

function Base.show(io::IO, ::MIME"text/plain", res::IrrRes)
  println(io, "Irreducible resolution ")
  println(io, "  ", join(["W^$i" for i in 0:(length(res.irrSums) - 1)], " -> "))
  println(io, "where ")
  for i in eachindex(res.irrSums)
    j = i-1
    println(io, " W^$j = direct sum of")
    for Ji in res.irrSums[i].components
      println(io, "    k{", Ji.vector, " + F - Q}_Q, where p_F = ", Ji.face.prime)
    end
  end
  println(io, "with maps")
  println(io, " W^{i-1} -> W^i for 1 <= i <= ", length(res.irrSums)-1)
  println(io, "of ", res.mod.mod)
  println(
    io,
    "over monoid algebra generated by ",
    join(gens(res.mod.monoidAlgebra.algebra), ", "),
    " quotient by ideal (",
    join(gens(modulus(res.mod.monoidAlgebra.algebra)), ", "),
    ").",
  )
end

function Base.show(io::IO, ::MIME"text/plain", Ji::IndecInj)
  println(io, "Indecomposable injective")
  println(io, "  k{", Ji.vector, " + F - Q},")
  print(io, "where p_F = ", Ji.face.prime, ".")
end

function ddirect_sum(I::InjMod...)
  kQ = I[1].MonoidAlgebra
  for J in I
    @assert J[1].monoidAlgebra == kQ "Monoid algebra not the same!"
  end
  return [J[2] for J in I]
end

#get Krull dimension of module
function dim(M::ModuleFP)
  ann = annihilator(M)
  return dim(ann)
end

#get (Krull) codimension of module
function codim(M::ModuleFP)
  R = base_ring(M)
  return dim(R) - dim(M)
end

@doc raw"""
    get_monoid_algebra(k_Q::Union{MPolyRing, MPolyQuoRing})

Computes polyhedral data associated to k_Q that is a monoid algebra over some affine semigroup Q: 
    - the polyhedral cone $C = \RR_{\geq 0}Q$
    - the faces $F$ of $C$ and the corresponding prime ideals $p_F= k\{Q\setminus F\}$
    - the hyperplanes bounding $C$
    - the zonotope $G_Q$ that is the Minkowski sum of all primitive integer vectors along rays of $Q$

# Examples
julia> S,(x,y,z) = graded_polynomial_ring(QQ,["x","y","z"]; weights = [[0,1],[1,1],[2,1]])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> J = ideal(S,[x*z-y^2])
Ideal generated by
  x*z - y^2

julia> R_Q,_ = quo(S,J)
(Quotient of multivariate polynomial ring by ideal (x*z - y^2), Map: S -> R_Q)

julia> get_monoid_algebra(R_Q)
Monoid algebra k[Q] over rational field generated by x, y, z quotient by ideal (x*z - y^2), 
where Q is a subset of ZZ^2 and the cone RR_{>= 0}Q has dimension 2.
"""
function get_monoid_algebra(k_Q::Union{MPolyRing,MPolyQuoRing})
  #check if monoid algebra
  @req is_zm_graded(k_Q) "Not ZZ^d-graded."
  gg_Q = grading_group(k_Q)
  @req is_free(gg_Q) && is_abelian(gg_Q) "Not a monoid algebra."

  # polyhedral cone RR_{\geq 0}Q
  C_Q = get_polyhedral_cone(k_Q)
  P_Q = polyhedron(C_Q) # C_Q as polyhedron
  G_Q, c = get_zonotope(P_Q)
  return MonoidAlgebra(
    k_Q,
    P_Q,
    get_faces_of_polyhedral_cone(k_Q, G_Q, P_Q),
    get_bounding_hyperplanes(P_Q),
    is_pointed(C_Q),
    (G_Q, c),
  )
end

@doc raw"""
    monoid_algebra_from_lattice(B::Union{Matrix{Int}, Vector{Vector{Int}}},k::Field)

Given a finite number of vectors $v_1,\dots,v_n$ in $\ZZ^d$ return the monoid algebra
\[k[x_1,\dots,x_n]/I_L,\]
where $I_L = (\bold{x}^u - \bold{x}^v \mid u,v \in \NN^d \text{ with } u -v \in L)$ is the lattice ideal of the lattice $L$ generated by $v_1,\dots,v_n$.

# Examples
```jldoctest
julia> monoid_algebra_from_lattice([[0,1],[1,1],[2,1]],QQ)

Monoid algebra k[Q] over rational field generated by x_1, x_2, x_3 quotient by ideal (-x_1*x_3 + x_2^2), 
where Q is a subset of ZZ^2 and the cone RR_{>= 0}Q has dimension 2.
```
"""
function monoid_algebra_from_lattice(B::Union{Matrix{Int},Vector{Vector{Int}}}, k::Field)
  if B isa Vector{Vector{Int}}
    genMatrix = Matrix{Int}(transpose(matrix(B)))
  else
    genMatrix = B
  end
  d = size(genMatrix, 1)

  # construct k[t_1,...,t_d]
  t_vars = [Symbol("t_$i") for i in 1:d]
  T, t = graded_polynomial_ring(k, t_vars)

  # construct k[x_1,...,x_n] where n is the number of columns/generators
  x_vars = [Symbol("x_$i") for i in 1:size(genMatrix, 2)]
  R, _ = graded_polynomial_ring(
    k, x_vars; weights=[Vector(row) for row in eachcol(genMatrix)]
  )

  # construct map x_i \to t^(a_i)
  targ = [prod(t[j]^genMatrix[j, i] for j in 1:d) for i in 1:size(genMatrix, 2)]
  map_T_R = hom(R, T, targ)

  # return monoid algebra 
  return get_monoid_algebra(quo(R, ideal(kernel(map_T_R)[1]))[1])
end

@doc raw"""
    get_monoid_algebra_module(kQ::MonoidAlgebra, M::SubquoModule)

Given a monoid algebra $k[Q]$ and a $k[Q]$-module M output M as a MonoidAlgebraModule. 
"""
function get_monoid_algebra_module(kQ::MonoidAlgebra, M::SubquoModule)
  @req base_ring(M) == kQ.algebra "Base rings do not match."
  return MonoidAlgebraModule(kQ, M)
end

function get_monoid_algebra_ideal(kQ::MonoidAlgebra, I::Ideal)
  @req base_ring(I) == kQ.algebra "Base rings do not match."
  return MonoidAlgebraIdeal(kQ, I)
end

function base_ring(M::MonoidAlgebraModule)
  return M.monoidAlgebra
end

function quotient_ring_as_module(I::MonoidAlgebraIdeal)
  return MonoidAlgebraModule(base_ring(I),quotient_ring_as_module(underlying_ideal(I)))
end

# given a face F of a the cone C = \RR_{\geq 0}Q of a monoid algebra, return the prime ideal k{Q\F}
# G_Q is the zonotope of Q, i.e. the Minkowski sum of all primitive integer vectors along rays of Q
function _get_prime_of_face(
  k_Q::Union{MPolyRing,MPolyQuoRing}, G_Q::Polyhedron, F::Polyhedron
)
  gens_PF = []
  for lp in lattice_points(G_Q)
    # if dim(intersect(convex_hull(lp),F))== -1 
    if !(lp in F) #check if lattice point is in F
      a_v = Vector{ZZRingElem}()
      for a_p in lp # PointVector -> Vector
        push!(a_v, a_p)
      end
      push!(gens_PF, monomial_basis(k_Q, a_v)[1])
    end
  end
  return ideal(k_Q, gens_PF)
end

# given a QQ^d vector v, this function returns a vector w that lies on the ray through v and the origin
# INPUT:    rational d-vector 
# OUTPUT:   integer d-vector
function rational_to_integer_vector(
  v::Union{Vector{QQFieldElem},AbstractVector{Rational},Vector{Rational}}
)
  denominators = [denominator(x) for x in v]
  lcm_denominators = lcm(denominators)

  return [Int(numerator(x*lcm_denominators)) for x in v]
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
  k_Q::Union{MPolyRing,MPolyQuoRing}, G_Q::Polyhedron, P::Polyhedron
)
  P_faces = Vector{Polyhedron}()
  for i in 0:dim(P)
    in
    append!(P_faces, faces(P, i))
  end
  return [FaceQ(_get_prime_of_face(k_Q, G_Q, F), F) for F in P_faces]
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

# given a hyperplane, return the H-presentation of it
# INPUT:    hyperplane h
# OUTPUT:   matrix A, vector b corresponding to Ax \leq b which defines h
function get_hyperplane_H_presentation(h::Polyhedron)
  aff_hull = affine_hull(h).Obj.pm_polytope.AFFINE_HULL
  _M = Matrix{Rational}(aff_hull)
  M = hcat(map(row -> reshape(rational_to_integer_vector(row), 1, :), eachrow(_M))...)
  A = [M[:, 2:n_columns(M)]; -M[:, 2:n_columns(M)]]
  b = [M[:, 1]; -M[:1]]
  return A, b
end

# given a polyhedral cone, return the zonotope as in Lemma 3.10 in [HM2004]
# INPUT:    polyhedral cone C  
# OUTPUT:   zonotope, sum of primitive integer vector ong rays of C 
function get_zonotope(P::Polyhedron)
  d = ambient_dim(P)
  P_rays = [rational_to_integer_vector(Vector(r)) for r in rays(P)]

  c = zeros(Int, d)
  zonotope = convex_hull(zeros(Int, d))
  for r in P_rays
    zonotope = zonotope + convex_hull([zeros(Int, 1, d); reshape(r, 1, length(r))])
    c = c + r
  end
  return zonotope, c
end

# given a $\ZZ^d$ graded module M and a vector $a \in ZZ^d$, return the module $M(-a)$ shifted by $a$
# INPUT:    $\ZZ^d$ graded module M,  vector $a \in ZZ^d$
# OUTPUT:   $M(-a)$ the module M shifted by a
function shifted_module(M::MonoidAlgebraModule, shift::Vector)
  S = M.monoidAlgebra.algebra
  m_shift = monomial_basis(S, shift)[1]

  A_shift = [m_shift*ambient_representative(m) for m in gens(M.mod)]
  B_shift = [m_shift*r for r in relations(M.mod)]

  return MonoidAlgebraModule(
    M.monoidAlgebra, SubquoModule(ambient_free_module(M.mod), A_shift, B_shift)
  )
end

# computes the generators of k{(a + H_+^°)\cap Q}
# Algorithm 3.11 in HM2004
# INPUT:    MonoidAlgebra kQ
#           HyperplaneQ H that bounding the polyhedral cone RR_{\geq 0}kQ
#           vector a in \ZZ^d
# OUTPUT:   finite set B such that (x^b : b \in B) equals k{(a + H_+^°)\cap Q}
function generators_W_H(kQ::MonoidAlgebra, H::HyperplaneQ, a::Vector{Int})
  F = intersect(H.hyperplane, kQ.cone)

  #get faces of Q intersecting F only at 0 in Q
  D = []
  for f in kQ.faces
    if dim(intersect(f.poly, F)) == 0
      push!(D, f.poly)
    end
  end

  B = []
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

@doc raw"
    compute_bass_numbers(M::SubquoModule,i::Int)

Let R = base_ring(M). Computes degrees of non-zero Bass numbers of $M$ up to some given cohomological degree i. 

INPUT:  monomial ideal I
        integer i
OUTPUT: list of ZZ^d graded degrees of Bass numbers up to cohomological degree i
"
function compute_bass_numbers(M::SubquoModule, i::Int) # now also for fin. gen. modules
  R_Q = base_ring(M)

  # residue_field
  I_m = ideal(R_Q, gens(R_Q))
  R_k = quotient_ring_as_module(I_m)

  D = Vector{Vector{Int}}()
  for j in 0:i
    in
    E = Nothing
    try
      E = ext(R_k, M, j)
    catch e
      if isa(e, AssertionError)
        E = Nothing
      end else
        for g in gens(E)
          push!(D, degree(Vector{Int}, g))
        end
    end
  end
  return unique(D) # filter duplicates
end

##check if all points lie in P_Q
##INPUT:    list of points in ZZ^d
##          polyhedron P_Q
##OUTPUT:   true if all points lie in P_Q
function points_in_Q(points::Vector{Vector{Int}}, P_Q::Polyhedron)
  for b in points
    if !is_subset(convex_hull(b), P_Q)
      return false
    end
  end
  return true
end

##computes a in ZZ^d such that all Bass numbers of M(a) lie in Q
##INPUT:    monomial ideal I
##          integer i up to which cohomological degree Bass numbers are considered
##OUTPUT:   ZZ^d degree of shift
# function compute_shift(I::MonoidAlgebraIdeal,i::Int)
function compute_shift(M::MonoidAlgebraModule, i::Int)
  n_bass = compute_bass_numbers(M.mod, i)
  c = M.monoidAlgebra.zonotope[2] #sum of all primitive integer vectors along rays of Q

  j = 0
  while !points_in_Q(n_bass, M.monoidAlgebra.cone) #loop until all degrees of bass numbers lie in Q
    bass_ = [a_bass + c for a_bass in n_bass]
    n_bass = bass_
    j = j + 1
  end
  return j*c
end

#computes the quotient of a module M by an ideal I, i.e. (0 :_M I)
#INPUT:     SubquoModule M
#           Ideal I
#OUTPUT:    SubquoModule (0 :_M I)
@doc raw"""
    mod_quotient(M::SubquoModule,I::Ideal)

This function computes the quotient \[(0 :_M I) := \{m \in M \mid m\cdot I = 0\}.\]

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

#compute saturation of a module M by an ideal, i.e. (0 :_M I^{infty})
#INPUT:     Submodule M
#           Ideal I
#OUTPUT:    SubquoModule (0 :_M I^{infty})
@doc raw"""
    mod_saturate(M::SubquoModule,I::Ideal)

This function computes the saturation \[(0 :_M I^\infty) := \{m \in M \mid m\cdot I^n = 0\text{ for some }n\in \NN_{>0}\}.\]

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

#Compute k[ZF]-basis as in Algorithm 3.3 of HM2004
#INPUT:     SubquoModule M = (0 :_M P_F) computed with mod_quotient(_,_)
#           prime ideal P_F = k{Q\F}
#OUTPUT:    k[ZF]-basis of localization (0 :_M P_F)[ZF]
function ZF_basis(M::SubquoModule, PF::Ideal)
  T = elem_type(M)
  #initilalize
  N = M
  h_N = identity_map(N)
  B = Vector{T}() # empty vector of k[ZF]-basis 

  for g in gens(M)
    N_g = sub(N, [h_N(g)])[1] #submodule of N =(0 :_M P_F)/(y0,...,yn) generated by g
    if annihilator(N_g) == PF && !is_zero(N_g)
      # if is_subset(PF,annihilator(N_g))&& !is_zero(N_g)
      push!(B, g)
    end
    M_B, _ = sub(M, B)
    N, h_N = quo(M, M_B) # update N
    if is_zero(N)
      break
    end
  end
  return filter(!is_zero, B)
end

#get the coefficient of a monomial 
function _f_k(m::Union{MPolyDecRingElem,MPolyQuoRingElem})
  R = parent(m)
  k = coefficient_ring(R)
  f_k = hom(R, k, ones(elem_type(k), ngens(R)))
  return f_k(m)
end

# this is an implementation of Algorithm 2 (fixed version of Algorithm 3.6. in HM2004)
function coefficients_fixed(N::SubquoModule, p::FaceQ, kQ::MonoidAlgebra)
  R_N = base_ring(N)
  @req R_N == base_ring(p.prime) "Base rings of module and ideal do not match."

  k = coefficient_ring(R_N) #get the field
  Np = mod_quotient(N, p.prime)[1]

  if is_zero(Np) # if Np is zero there is nothing to compute
    return [], zeros(R_N, 1, 1)
  end

  # get socle degrees of indecomposable injectives kQ{a + F - Q}, i.e. compute a k[F]-basis of the localisation (0 :_M P_F)[ZZ F]
  Bp = ZF_basis(Np, p.prime)

  R = relations(N) #get all relations of N
  # n = ngens(N)
  n = ngens(ambient_free_module(N))

  lambda = []
  for b in Bp
    #get coefficient vector of w.r.t. generators of M
    b_amb = ambient_representative(b)
    _c_b = coordinates(b_amb) #coordinates w.r.t. ambient free module
    _c_b = coordinates(N(b_amb)) #coordinates w.r.t. generators of M 
    c_b = [_f_k(_c_b[i]) for i in 1:ngens(N)]

    #get all relevant generators of N, i.e., check (deg(b) + F) \cap (deg(g) + Q) ≠ ∅
    b_p = convex_hull(degree(Vector{Int}, b))
    G_b = Vector{SubquoModuleElem}()
    for g_N in filter(!is_zero, gens(N))
      g_p = convex_hull(degree(Vector{Int}, g_N))
      if dim(intersect(b_p + p.poly, g_p + kQ.cone)) >= 0
        push!(G_b, g_N)
      end
    end
    x_Gb = [monomial_basis(R_N, degree(g))[1] for g in G_b]
    _N = sub(ambient_free_module(N), [ambient_representative(g) for g in G_b])[1]

    #get all b-relevant relations w.r.t. F
    R_bF = []
    C_bF = []
    for r in R
      #check (deg(b) + F)\cap (deg(r) + Q) ≠ ∅
      r_p = convex_hull(degree(Vector{Int}, r))
      if dim(intersect(b_p + p.poly, r_p + kQ.cone)) >= 0
        x_r = monomial_basis(R_N, degree(r))[1]
        a = lcm(x_Gb..., x_r)
        _r = (a//x_r).num*r #well-defined since a is lcm
        try #can _r be written using generators in G_b
          _c_r = coordinates(_N(_r))
          c_r = Vector{elem_type(k)}()
          for i in 1:ngens(N)
            j = findfirst(g -> g == N[i], G_b)
            if j !== nothing
              push!(c_r, _f_k(_c_r[j]))
            else
              push!(c_r, k())
            end
          end
          # c_r = [_f_k(_c_r[i]) for i = 1:ngens(N)]
          push!(C_bF, c_r)
          push!(R_bF, r)
        catch
        end
      end
    end

    #get kernel of coefficient matrix -> K
    if C_bF != []
      _K = matrix(QQ, hcat(C_bF...))
      K = kernel(_K)
    else
      K = identity_matrix(QQ, ngens(N))
    end

    #get kernel of coefficient vector of b -> B
    _B = matrix(QQ, 1, length(c_b), c_b)
    B = kernel(transpose(_B))

    #get \lambda_b
    rows_K = [K[i, :] for i in 1:nrows(K)]
    l = rank(B)
    possible_rows = []
    for c_K in rows_K
      B_K = vcat(B, matrix(QQ, 1, ngens(N), c_K))
      if rank(B_K)>l
        push!(possible_rows, c_K)
      end
    end
    push!(lambda, possible_rows[1])
  end
  return Bp, map(x -> x*one(R_N), hcat(lambda...))
end

@doc raw"""
    irreducible_hull(M::SubquoModule,P::Vector{FaceQ})

Computes an irreducible hull of a $\ZZ^d$-graded $k[Q]$-module M. 

INPUT:  SubquoModule M over semigroup ring k[Q]
        List of prime ideals corresponding to the faces of Q
OUTPUT: effective irreducible hull W_bar and effective vector set lambda
"""
function irreducible_hull(Mi::SubquoModule, P::Vector{FaceQ}, kQ::MonoidAlgebra, j=0)
  N = Mi
  summands = Vector{IndecInj}()
  lambda = []

  # for p in filter(p -> !is_zero(p.prime),P) 
  for p in P
    # Bp,lambda_p = _coefficients(N,p,kQ) #compute k[ZZF]-basis of (0 :_M P_F) and the scalar matrix that ensures these basis-elements are mapped to something non-zero in W^i (that is not completely constructed)
    Bp, lambda_p = coefficients_fixed(N, p, kQ)
    for b in Bp
      push!(summands, IndecInj(p, degree(Vector{Int}, b)))
    end

    if length(Bp) > 0 # we don't want to add zero vectors to lambda...
      push!(lambda, lambda_p)
    end
    # M_sat = mod_saturate(Mi,p.prime)
    M_sat = saturation((ideal(kQ.algebra, []) * Mi)[1], p.prime)
    # if !is_zero(saturation((ideal(kQ.algebra,[])*Mi)[1],p.prime))
    if !is_zero(p.prime) && !is_zero(M_sat)
      # N,_ = quo(Mi,saturation((ideal(kQ.algebra,[])*Mi)[1],p.prime)) # use new saturation method
      N, _ = quo(Mi, M_sat)
    end
    # N,_ = quo(Mi,mod_saturate(Mi,p.prime)) #update N
    if is_zero(N)
      break
    end
  end
  # return summands,foldl(hcat,lambda), degLambda
  return summands, hcat(lambda...)
end

#get irreducible decomposition of ideal over monoid algebra
@doc raw"""
    irreducible_dec(I::MonoidAlgebraIdeal)

Let $R$ be a monoid algebra and $I\subset R$ an ideal. This function computes an irreducible decomposition of $I$, i.e.
\[I = W_1 \cap \dots \cap W_r,\]
where $W_i\subset R$ are irreducible ideals for all $i$.   

# Examples
```jldoctest
julia> R_Q,(x,y) = graded_polynomial_ring(QQ,["x","y"];weights = [[1,0],[0,1]])
(Graded multivariate polynomial ring in 2 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y])

julia> kQ = get_monoid_algebra(R_Q)
Monoid algebra k[Q] over rational field generated by x, y quotient by ideal (), 
where Q is a subset of ZZ^2 and the cone RR_{>= 0}Q has dimension 2.

julia> I = ideal(kQ,[x^4,x^2*y^2,y^4])
Ideal generated by
 x^4
 x^2*y^2
 y^4
over monoid algebra generated by x, y quotient by ideal ().


julia> irreducible_dec(I)
2-element Vector{MPolyIdeal{MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}}}:
 Ideal (x^2, x^2*y, y^4, x*y^4)
 Ideal (x^4, x^4*y, y^2, x*y^2)
```
"""
function irreducible_dec(I::MonoidAlgebraIdeal)
  kQ = base_ring(I)

  indec_injectives, _ = irreducible_hull(quotient_ring_as_module(I.ideal), kQ.faces, kQ)
  return [_get_irreducible_ideal(kQ, I) for I in indec_injectives]
end

#given a monoid algebra k[Q] and an indecomposable injective J = k{a + F - Q} compute the irreducible ideal W such that J_Q = k[Q]/W 
function _get_irreducible_ideal(kQ::MonoidAlgebra, J::IndecInj)
  B_i = []

  for h in hyperplanes(kQ)
      if is_subset(J.face.poly,h.hyperplane)
          B_h = generators_W_H(kQ,h,J.vector)
          push!(B_i,B_h)
      end
  end
  end

  G_W = Vector{MPolyDecRingElem{QQFieldElem,QQMPolyRingElem}}()
  for b in B_i
    for bb in b
      a_v = Vector{ZZRingElem}()
      for a in bb
        push!(a_v, a)
      end
      push!(G_W, monomial_basis(kQ.algebra, a_v)[1])
    end
  end
  return ideal(kQ.algebra, G_W)
end

@doc raw"""
    irreducible_res(M::MonoidAlgebraModule, i::Int = 0)

Let $k[Q]$ be a monoid algebra and let $M$ be a $\ZZ^d$-graded module over $k[Q]$. This function computes a minimal irreducible resolution
\[0 \to M \xrightarrow{ϵ} \overline{W}^0 \xrightarrow{d^0} \overline{W}^1 \dots \xrightarrow{d^{r-1} \overline{W}^r,\]
where $\overline{W}^i = \sum_{j=1}^{n_i} \overline{W_{i_j}} = \sum_{j=1}^{n_i} k[Q]/W_{i_j}$ for irreducible ideals $W_{i_j}$. 

# Examples 
```jldoctest
julia> R_Q,(x,y) = graded_polynomial_ring(QQ,["x","y"];weights = [[1,0],[0,1]])
(Graded multivariate polynomial ring in 2 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y])

julia> kQ = get_monoid_algebra(R_Q)
Monoid algebra k[Q] over rational field generated by x, y quotient by ideal (), 
where Q is a subset of ZZ^2 and the cone RR_{>= 0}Q has dimension 2.

julia> I = ideal(kQ,[x^4,x^2*y^2,y^4])
Ideal generated by
 x^4
 x^2*y^2
 y^4
over monoid algebra generated by x, y quotient by ideal ().


julia> M = quotient_ring_as_module(I)
Graded subquotient of graded submodule of R_Q^1 with 1 generator
  1: e[1]
by graded submodule of R_Q^1 with 3 generators
  1: x^4*e[1]
  2: x^2*y^2*e[1]
  3: y^4*e[1], 
over monoid algebra generated by x, y quotient by ideal ().

julia> irr_res = irreducible_res(M)
Irreducible resolution 
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
over monoid algebra generated by x, y quotient by ideal ().
```
"""
function irreducible_res(M::MonoidAlgebraModule, i::Int = 0)
    kQ = M.monoidAlgebra
    R_Q = kQ.algebra
    Mi = M.mod # current module in resolution
    gi = identity_map(Mi) #initilalize
    res_Wi = Vector{IrrSum}()
    res_Mi = [Mi] #cokernels 
    res_hi = Vector{SubQuoHom}()
    res_fi = Vector{SubQuoHom}()
    res_gi = [gi] #quotient maps

    j = 1
    while !is_zero(Mi) #until cokernel Mi is zero
        #compute irreducible hull
        indec_injectives, _lambda = irreducible_hull(Mi,faces(kQ),kQ,j)

        #get the corresponding irreducible ideal for each indecomposable injective
        irreducible_ideals = [_get_irreducible_ideal(kQ,J) for J in indec_injectives]

        #get irreducible sum Wi, i.e. the direct sum of all quotients k[Q]/K_i, where K_i is an irreducible ideal
        irreducible_comp = map(I -> quotient_ring_as_module(I),irreducible_ideals)
        d_sum(x,y) = direct_sum(x,y,task=:none)
        Wi = foldl(d_sum,irreducible_comp)

        #multiply rows of lambda by degrees of generators of Mi
        m,n = size(_lambda)
        lambda = zero(_lambda)
        for i in 1:m
            for j in 1:n
                lambda[i, j] = _lambda[i, j] * monomial_basis(R_Q,degree(Mi[i]))[1]
            end
        end

        #define injective map Mi -> Wi
        fi = hom(Mi,Wi,matrix(lambda))

        #get boundary map W{i-1} -> Wi
        hi = gi*fi
        
        #compute cokernel and then simplify
        Mi_,gi_ = quo(Wi,image(hi)[1]) #cokernel
        Mi,ji = prune_with_map(Mi_)
        gi = gi_*inv(ji)
        
        if length(filter(is_zero,relations(Mi))) > 0 # fix modules with "zero" relations
            Mi,h = fix_module(Mi)
            gi = gi*h
        end
        push!(res_Mi,Mi)
        push!(res_gi,gi)
        push!(res_Wi,IrrSum(Wi,indec_injectives))
        push!(res_hi,hi)
        push!(res_fi,fi)

        # end at cohomological degree i
        if i > 0 && j == i
            break
        end
        j = j + 1
    end

    #define injective map Mi -> Wi
    fi = hom(Mi, Wi, matrix(lambda))

    #get boundary map W{i-1} -> Wi
    hi = gi*fi

    #compute cokernel and then simplify
    Mi_, gi_ = quo(Wi, image(hi)[1]) #cokernel
    Mi, ji = prune_with_map(Mi_)
    gi = gi_*inv(ji)

    if length(filter(is_zero, relations(Mi))) > 0 # fix modules with "zero" relations
      Mi, h = fix_module(Mi)
      gi = gi*h
    end
    push!(res_Mi, Mi)
    push!(res_gi, gi)
    push!(res_Wi, IrrSum(Wi, indec_injectives))
    push!(res_hi, hi)
    push!(res_fi, fi)

    # end at cohomological degree i
    if i > 0 && j == i
      break
    end
    j = j + 1
  end

  #get cochain complex
  C = cochain_complex([identity_map(M.mod)]) # default value if sequence not exact
  try
    C = cochain_complex(res_hi)
  catch
    ;
  end

  return IrrRes(M, res_Wi, res_hi, res_gi, res_fi, res_Mi, C)
end

# compute a minimal injective resolution up to cohomological degree i
# INPUT:    MonoidAlgebraIdeal I
#           positive integer i
# OUTPUT:   injective resolution of R/I up to cohomological degree I
@doc raw"""
    injective_res(M::MonoidAlgebraModule, i::Int)

Let $k[Q]$ be a monoid algebra and $M$ a finitely generated $Q$-graded module over $k[Q]$. This function computes a minimal injective
resolution 
\[0 \to M \xrightarrow{ϵ} J^0 \xrightarrow{d^0} J^1 \xrightarrow{d^1} \dots \xrightarrow{d^{i-1}} J^i.\]
The maps $d^j$ are given by monomial matrices. This is an implementation of the algorithms in~\cite{HM2004}.

# Examples
```jldoctest
julia> R_Q,(x,y) = graded_polynomial_ring(QQ,["x","y"];weights = [[1,0],[0,1]])
(Graded multivariate polynomial ring in 2 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y])

julia> kQ = get_monoid_algebra(R_Q)
Monoid algebra k[Q] over rational field generated by x, y quotient by ideal (), 
where Q is a subset of ZZ^2 and the cone RR_{>= 0}Q has dimension 2.

julia> I = ideal(kQ,[x^4,x^2*y^2,y^4])
Ideal generated by
 x^4
 x^2*y^2
 y^4
over monoid algebra generated by x, y quotient by ideal ().


julia> M = quotient_ring_as_module(I)
Graded subquotient of graded submodule of R_Q^1 with 1 generator
  1: e[1]
by graded submodule of R_Q^1 with 3 generators
  1: x^4*e[1]
  2: x^2*y^2*e[1]
  3: y^4*e[1], 
over monoid algebra generated by x, y quotient by ideal ().

julia> injective_res(M,2)
Injective resolution 
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
over monoid algebra generated by x, y quotient by ideal ().
```

```jldoctest
julia> S,(x,y,z) = graded_polynomial_ring(QQ,["x","y","z"]; weights = [[0,1],[1,1],[2,1]])
(Graded multivariate polynomial ring in 3 variables over QQ, MPolyDecRingElem{QQFieldElem, QQMPolyRingElem}[x, y, z])

julia> J = ideal(S,[x*z-y^2])
Ideal generated by
  x*z - y^2

julia> kQ = get_monoid_algebra(quo(S,J)[1])
Monoid algebra k[Q] over rational field generated by x, y, z quotient by ideal (x*z - y^2), 
where Q is a subset of ZZ^2 and the cone RR_{>= 0}Q has dimension 2.

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

julia> M = get_monoid_algebra_module(kQ,SubquoModule(F,a,b))
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
Injective resolution 
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
over monoid algebra generated by x, y, z quotient by ideal (x*z - y^2).
```
"""
function injective_res(M::MonoidAlgebraModule, i::Int)
  kQ = M.monoidAlgebra

  #compute irreducible resolution of shifted module
  a_shift = compute_shift(M, i+1)
  irrRes = irreducible_res(shifted_module(M, a_shift))

  #get injective modules up to cohomological degree i, i.e. J^0, J^1, ...,J^i
  inj_modules = Vector{InjMod}()
  for j in 1:min(i + 1, length(irrRes.irrSums))
    shifted_comp = map(
      indec -> IndecInj(indec.face, indec.vector - a_shift), irrRes.irrSums[j].components
    )
    push!(inj_modules, InjMod(kQ, shifted_comp))
  end

  #get all needed maps (as k-matrix or k[Q]-matrix?)
  cochain_maps = [
    matrix(irrRes.cochainMaps[k]) for k in eachindex(irrRes.cochainMaps) if 1 < k <= i+1
  ]
  return InjRes(M, inj_modules, cochain_maps, length(inj_modules)-1, irrRes, a_shift)
end

@doc raw"""
    injective_res(I::MonoidAlgebraIdeal,i::Int)
Let $k[Q]$ be a monoid algebra and $I\subset k[Q]$ a $\ZZ^d$-graded ideal. This function computes an injective resolution of $M = k[Q]/I$.  
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

function jacobian(I::Ideal)
  R = base_ring(I)
  gens_I = gens(I)

  if is_zero(I)
    return Matrix{elem_type(R)}(undef, 0, 0)
  else
    jacobian = Matrix{elem_type(R)}(undef, length(gens_I), length(gens(R)))

    for i in eachindex(gens_I)
      for j in eachindex(gens(R))
        jacobian[i, j] = derivative(gens_I[i], j)
      end
    end

    return jacobian
  end
end

function jacobian(R::Union{MPolyRing,MPolyQuoRing})
  if R isa MPolyQuoRing
    return jacobian(modulus(R))
  else
    return Matrix{elem_type(R)}(undef, 0, 0)
  end
end

function is_normal(R::Union{MPolyRing,MPolyQuoRing})
  # 1 argument: A ring - usually a quotient ring.
  # Return: A boolean value, true if the ring is normal and false otherwise.
  if R isa MPolyQuoRing
    I = R.I
    R_B = base_ring(R)
  else
    I = ideal(R_Q, []) #zero ideal
    R_B = R
  end
  # Get the ideal associated with the ring
  # M = cokernel(generators(I)) # Compute the cokernel of the generators of the ideal
  M = quotient_ring_as_module(I)
  n = codim(I)                # Calculate the codimension of the ideal

  # Check the S2 condition
  test = 0:(dim(R) - n - 2)

  # test = [dim(ring(I)) - n - 1 - i for i in 0:(dim(ring(I)) - n - 1)]

  if all([
    (codim(ext(M, graded_free_module(R_B, 1), j + n + 1)) >= (j + n + 3)) for j in test
  ])
    Jac = minors(n, jacobian(R))  # Compute minors of the Jacobian
    d = dim(Jac)                   # Get dimension of the Jacobian
    d < 0 && (d = -Inf)            # Handle negative dimensions
    return (dim(R) - d >= 2)       # Check the condition
  else
    return false                    # S2 condition not satisfied
  end
end

# # Algorithm 3.6. in HM2004 that is not working for the general case
# function coefficients(N::SubquoModule, p::FaceQ)
#     R_N = base_ring(N)
#     @req R_N == base_ring(p.prime) "Base rings of module and ideal do not match."

#     T = elem_type(R_N)
#     k = coefficient_ring(R_N)
#     Np,incl_Np = mod_quotient(N,p.prime) # quotient (0 :_N p.prime)
#     if is_zero(Np)
#         return Np,zeros(R_N,1,1)
#     end

#     # get socle degrees of indecomposable injectives kQ{a + F - Q}
#     Bp = ZF_basis(Np,p.prime)

#     #### (maybe) new better version
#     _lambda = Vector{Vector}()
#     m_g = one(R_N)
#     for b in Bp
#         lambda_b = []
#         for i=1:ngens(N)
#             m_g = monomial_basis(R_N,degree(N[i]))[1]*one(R_N)
#             m_bg = monomial_basis(R_N,degree(b)-degree(N[i]))[1]
#             if !is_zero(m_bg*N[i]) 
#                 b_rest = zero(N)
#                 for j = 1:ngens(N)
#                     if length(coordinates(incl_Np(b)-b_rest)) == 1
#                         break
#                     end 
#                     if j == i || is_zero(coordinates(ambient_representative(b))[j])
#                         continue
#                     end
#                     mu_j = coordinates(ambient_representative(b))[j]
#                     if !is_zero(mu_j)
#                         b_rest = b_rest + mu_j*N[j]
#                     end
#                 end
#                 Ni,incl_i = sub(N,[N[i]])
#                 _b = Nothing
#                 try
#                     _b = preimage(incl_i, incl_Np(b) - b_rest)
#                 catch
#                 end
#                 if _b != Nothing
#                     push!(lambda_b,leading_coefficient(lift(coordinates(_b)[1])))
#                 else
#                     push!(lambda_b,k())
#                 end
#             else
#                 push!(lambda_b,k())
#             end
#         end
#         push!(_lambda,lambda_b)
#     end
#     return Bp, transpose(matrix(k,_lambda)), m_g
# end

# # given a fin. gen. module and a prime ideal p, compute a k[ZZF]-basis Bp of (0 :_M p) and for each generator b in Bp compute the scalar matrix lambda that defines a well-defined injective map 
# function _coefficients(N::SubquoModule, p::FaceQ, kQ::MonoidAlgebra)
#     R_N = base_ring(N)
#     @req R_N == base_ring(p.prime) "Base rings of module and ideal do not match."

#     k = coefficient_ring(R_N) #get the field
#     Np = mod_quotient(N,p.prime)[1]

#     if is_zero(Np) # if Np is zero there is nothing to compute
#         return [],zeros(R_N,1,1)
#     end

#     # get socle degrees of indecomposable injectives kQ{a + F - Q}
#     Bp = ZF_basis(Np,p.prime)

#     R = relations(N) #get all relations of N
#     # n = ngens(N)
#     n = ngens(ambient_free_module(N))

#     lambda = []
#     for b in Bp # for every basis vector of (0 :_M P_F) we compute a scalar vector \lambda_b 
#         # m_b = monomial_basis(R_N,degree(b))[1]
#         lambda_b = []
#         b_c = dense_row(coordinates(ambient_representative(b)),n)[1,:]

#         # check in which positions (i.e. coordinates of free presentation) b is non-zero
#         m = ngens(N)
#         s = [] #positions of non-zero terms of b
#         for i=1:n
#             # if is_zero(b_c[i]*N[i])
#             # if is_zero(b_c[i]) || all((degree(Vector{Int},b_c[i]) - degree(Vector{Int},N[i])) .>= degree(Vector{Int},one(R_N))) && is_zero(monomial_basis(R_N,degree(Vector{Int},b_c[i])-degree(Vector{Int},N[i]))[1]*N[i]) # check if the coefficient of b[i] is already zero or if the term b[i] is zero (b[i] = x^(deg(b_c[i]-deg(N[i]))*N[i]) 
#             set_zero = true
#             for j = 1:m
#                 # if is_zero(coordinates(b)[j]) ## is this working???
#                 #     continue
#                 # end
#                 if is_zero(b_c[i]) || (points_in_Q([degree(Vector{Int},b) - degree(Vector{Int},N[j])],kQ.cone) && is_zero(monomial_basis(R_N,degree(Vector{Int},b)-degree(Vector{Int},N[j]))[1]*N[j])) # check if the coefficient of b[i] is already zero or if the term b[i] is zero (b[i] = x^(deg(b_c[i]-deg(N[i]))*N[i]) 
#                     # b_c[i] = zero(R_N)
#                 else #term b[i] is non-zero
#                     push!(s,i)
#                     set_zero = false
#                     break
#                 end 
#             end
#             if set_zero
#                 b_c[i] = zero(R_N)
#             end
#         end

#         # get all relations that apply in deg(b)
#         rels_b = [r for r in R if points_in_Q([degree(Vector{Int},b) - degree(Vector{Int},r)],kQ.cone)]
#         r = [] # we are only interested in non-trivial relations, i.e. not of the form x^a*N[i] = 0
#         # also we want a reduced presentation of these relations, i.e. the coefficient of r[i] is zero if r[i] = 0  
#         for r_b in rels_b 
#             _r = dense_row(coordinates(ambient_representative(r_b)),n)[1,:] # get all coefficients of the relation
#             s_r = [] #positions of non-zero terms of r_b 
#             for i=1:n 
#                 #check if entry is already non-zero, i.e. not relevant in the relation... 
#                 # if !is_zero(m_b*N[i])
#                 if !is_zero(_r[i]) ## maybe here something goes wrong!!!
#                     if !is_zero(monomial_basis(R_N,degree(Vector{Int},b)-degree(Vector{Int},N[i]))[1]*N[i]) # check if term is non-zero
#                         push!(s_r,i)
#                     else
#                         _r[i] = zero(R_N)
#                     end 
#                 end
#             end
#             if length(intersect(s,s_r)) > 0 && count(!is_zero,_r) > 1  #avoid simple relations x^a*e_i = 0 and we only want relations that are relevant for b, i.e. that have a common non-zero term
#                 push!(r,(r_b,_r))
#             end
#         end

#         # compute coefficients in \lambda_b
#         if r != []
#             rel = r[argmin([degree(Vector{Int},m[1]) for m in r])][2] #get the relevant relation of lowest degree    
#             last_non_zero = findlast(!is_zero,rel) #index of last non-zero term

#             sum_of_coeffs = k() #needed to compute the coefficient of last non-zero entry
#             # f_k = hom(R_N,k,ones(elem_type(k),ngens(R_N))) #get coefficient
#             for i=1:n #set all (except the last) entry to (one if the tern r[i] is non-zero) or (zero if the term is zero)  
#                 if i == last_non_zero #set the last entry such that rel*lambda_b = 0 (dot product of coefficients of rel with lambda_b)
#                     coeff = -(sum_of_coeffs//_f_k(rel[i]))
#                     push!(lambda_b,coeff)
#                 else
#                     if !is_zero(rel[i])
#                         push!(lambda_b,one(R_N))
#                         sum_of_coeffs = sum_of_coeffs + _f_k(rel[i]) #add scalar coefficient of rel[i]
#                     else
#                         push!(lambda_b,zero(R_N))
#                     end 
#                 end
#             end
#         else # no non-trivial relations in deg(b)
#             l = findfirst(!is_zero,b_c) # first non-zero term of b 
#             for i=1:n 
#                 if i == l # alternatively check if i in s
#                     push!(lambda_b,one(R_N)) # if suffices to set the first non-zero entry to one(R_N)
#                 else
#                     push!(lambda_b,zero(R_N))
#                 end 
#             end
#         end
#             push!(lambda,lambda_b)
#         end
#     # end
#     if length(Bp) == 0
#         return [],[]
#     end
#     return Bp, matrix(R_N,hcat(lambda...))
# end

# ## not used???
# # checks if two ring elements are equal modulo ZF
# # INPUT:     RingElem a
# #            RingElem b
# #            prime ideal p_F 
# # OUTPUT:    true if a is congruent to b modulo ZF
# function equal_mod_ZF(a::RingElem,b::RingElem,p_F::Ideal)
#     # @assert parent(a) == parent(b) && parent(a) == base_ring(I) "base ring and parents must match"
#     if degree(a) == degree(b)
#         return true 
#     elseif !(div(a,b) == 0)
#         return !ideal_membership(div(a,b),p_F)
#     elseif !(div(b,a) == 0)
#         return !ideal_membership(div(b,a),p_F)
#     else
#         return false
#     end
# end

# import local cohomology functions
include("LocalCohomology.jl")
