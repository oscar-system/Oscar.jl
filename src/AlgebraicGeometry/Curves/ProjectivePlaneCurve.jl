
################################################################################
@doc raw"""
    ProjectivePlaneCurve <: AbsProjectiveCurve

A reduced curve in the projective plane.

# Examples
```jldoctest
julia> R, (x,y,z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> C = plane_curve(y^3*x^6 - y^6*x^2*z)
Projective plane curve
  defined by 0 = x^5*y - x*y^4*z

```
"""
@attributes mutable struct ProjectivePlaneCurve{BaseRingType<:Field, RingType<:Ring} <: AbsProjectiveCurve{BaseRingType, RingType}
  X::ProjectiveAlgebraicSet{BaseRingType, RingType}
  defining_equation::MPolyDecRingElem

  function ProjectivePlaneCurve(X::ProjectiveAlgebraicSet{S,T}, check::Bool=true) where {S,T}
    @check begin
      dim(X) == 1 || error("not of dimension one")
      dim(ambient_space(X)) == 2 || error("not a plane curve")
    end
    new{S,T}(X)
  end

  function ProjectivePlaneCurve(eqn::MPolyDecRingElem; is_radical::Bool=false)
    ngens(parent(eqn)) == 3 || error("not a plane curve")
    if !is_radical
      eqn = prod([i[1] for i in factor(eqn)], init=one(parent(eqn)))
    end
    X = ProjectivePlaneCurve(algebraic_set(eqn; is_radical=true))
    X.defining_equation = eqn
    return X
  end
end



ProjectivePlaneCurve(I::MPolyIdeal{<:MPolyDecRingElem}; kwargs...) = ProjectivePlaneCurve(algebraic_set(eq; kwargs...))

plane_curve(I::MPolyIdeal{<:MPolyDecRingElem};kwargs...) = ProjectivePlaneCurve(I; kwargs...)
plane_curve(eq::MPolyDecRingElem{<:FieldElem}; kwargs...) = ProjectivePlaneCurve(eq; kwargs...)

underlying_scheme(X::ProjectivePlaneCurve) = X.X
fat_scheme(X::ProjectivePlaneCurve) = fat_scheme(underlying_scheme(X))


function Base.show(io::IO, ::MIME"text/plain", C::ProjectivePlaneCurve)
  io = pretty(io)
  println(io, "Projective plane curve")
  print(io, Indent(), "defined by 0 = ", defining_equation(C), Dedent())
end


################################################################################
# Plane Curves related functions.
################################################################################

@doc raw"""
    defining_equation(C::ProjectivePlaneCurve)

Return the defining equation of the (reduced) plane curve `C`.
"""
function defining_equation(C::ProjectivePlaneCurve{S,MPolyQuoRing{T}}) where {S,T}
  if isdefined(C, :defining_equation)
    return C.defining_equation::T
  end
  m = minimal_generating_set(vanishing_ideal(C))
  @assert length(m) == 1
  return m[1]::T
end

################################################################################
# hash function

function Base.hash(C::ProjectivePlaneCurve, h::UInt)
  F = 1//AbstractAlgebra.leading_coefficient(defining_equation(C))*defining_equation(C)
  return hash(F, h)
end

################################################################################

@doc raw"""
    degree(C::ProjectivePlaneCurve)

Return the degree of the defining polynomial of `C`.
"""
@attr Int function degree(C::ProjectivePlaneCurve)
  return total_degree(defining_equation(C))
end

union(X::T,Y::T) where {T<:ProjectivePlaneCurve} = ProjectivePlaneCurve(union(underlying_scheme(X),underlying_scheme(Y)))

@doc raw"""
    common_components(C::S, D::S) where {S<:ProjectivePlaneCurve}

Return the projective plane curve consisting of the common components of `C`
and `D`, or an empty vector if they do not have a common component.
"""
function common_components(C::S, D::S) where {S<:ProjectivePlaneCurve}
  G = gcd(defining_equation(C), defining_equation(D))
  if isone(G)
    return Vector{S}()
  else
    return S[plane_curve(G)]
  end
end

function irreducible_components(C::ProjectivePlaneCurve)
  return [ProjectivePlaneCurve(i[1]) for i in factor(defining_equation(C))]
end





################################################################################
# multiplicity

@doc raw"""
     multiplicity(C::ProjectivePlaneCurve{S}, P::AbsProjectiveRationalPoint)

Return the multiplicity of `C` at `P`.
"""
function multiplicity(C::ProjectivePlaneCurve, P::AbsProjectiveRationalPoint)
  P in C || return 0
  P = C(coordinates(P))
  S = standard_covering(C)
  i = findfirst(!iszero, coordinates(P))
  Ca = S[i]
  Q = Ca(P)
  return multiplicity(Ca, Q)
end


################################################################################
# homogeneization for lines

function help_homogene_line(R::MPolyRing, r::MPolyRing, F::MPolyRingElem, i::Int)
  total_degree(F) == 1 || error("This is not a degree one polynomial")
  V = gens(R)
  W = gens(r)
  v = V[i]
  deleteat!(V, i)
  phi = hom(r, R, V)
  G = phi(F)
  G = G - evaluate(G, [0,0,0])*(1-v)
  return G
end

################################################################################
# tangent lines

 @doc raw"""
      tangent_lines(C::ProjectivePlaneCurve{S}, P::Geometry.ProjSpcElem{S}) where S <: FieldElem

Return the tangent lines at `P` to `C` with their multiplicity.
"""
function tangent_lines(C::ProjectivePlaneCurve, P::AbsProjectiveRationalPoint)
  P in C || error("The point is not on the curve.")
  P = C(P)
  R = ambient_coordinate_ring(C)
  S = standard_covering(C)
  i = findfirst(!iszero, coordinates(P))
  Ca = S[i]
  Q = Ca(P)
  L = tangent_lines(Ca, Q)
  D = Dict{ProjectivePlaneCurve, Int}()
  if !isempty(L)
    D = Dict(ProjectivePlaneCurve(help_homogene_line(R, ambient_coordinate_ring(Ca), defining_equation(x), i)) => L[x] for x in keys(L))
  end
  return D
end

@doc raw"""
     intersection_multiplicity(C::S, D::S, P::AbsProjectiveRationalPoint) where S <: ProjectivePlaneCurve

Return the intersection multiplicity of `C` and `D` at `P`.
"""
function intersection_multiplicity(C::S, D::S, P::AbsProjectiveRationalPoint) where S <: ProjectivePlaneCurve
  P in ambient_space(C) || error("The point needs to be in a projective two dimensional space")
  SC = standard_covering(C)
  SD = standard_covering(D)
  i = findfirst(!iszero, coordinates(P))
  return intersection_multiplicity(SC[i], SD[i], SC[i](dehomogenization(P,i)))
end

################################################################################


@doc raw"""
     is_transverse_intersection(C::S, D::S, P::AbsProjectiveRationalPoint) where S <: ProjectivePlaneCurve

Return `true` if `C` and `D` intersect transversally at `P` and `false` otherwise.
"""
function is_transverse_intersection(C::S, D::S, P::AbsProjectiveRationalPoint) where S <: ProjectivePlaneCurve
  P in C && P in D || return false
  any(P in i for i in common_components(C,D)) && return false
  intersection_multiplicity(C, D, P) == 1
end

################################################################################

@doc raw"""
    arithmetic_genus(C::ProjectivePlaneCurve)

Return the arithmetic genus of `C`.

# Examples
```jldoctest
julia> T, (x, y, z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> C = plane_curve(y^2 * z - x^3 - x * z^2)
Projective plane curve
  defined by 0 = x^3 + x*z^2 - y^2*z

julia> arithmetic_genus(C)
1
```
"""
function arithmetic_genus(C::ProjectivePlaneCurve)
  H = hilbert_polynomial(C)
  return -ZZ(coeff(H, 0)) + 1
end

################################################################################

@doc raw"""
    geometric_genus(C::ProjectivePlaneCurve; check::Bool=true)

Return the geometric genus of `C`.

If `C` is singular this is defined as the geometric genus of any smooth birational model of `C`.

If `check` is true, checks that `C` is an irreducible curve.

# Examples
```jldoctest
julia> R, (x,y,z) = graded_polynomial_ring(QQ, ["x", "y", "z"]);

julia> C = plane_curve(z*x^2-y^3)
Projective plane curve
  defined by 0 = x^2*z - y^3

julia> geometric_genus(C)
0

```
"""
geometric_genus(C::AbsProjectiveCurve)

function geometric_genus(C::ProjectivePlaneCurve{<:QQField}; check::Bool=true)
  @check is_irreducible(C) "the curve must be irreducible"
  I = vanishing_ideal(C)
  singular_assure(I)
  return ZZ(Singular.LibNormal.genus(I.gens.S)::Int)
end

function geometric_genus(C::ProjectivePlaneCurve)
  # buggy
  A = homogeneous_coordinate_ring(C)
  L = normalization(A)
  m = length(L)
  pa = zero(ZZ)
  for i in 1:m
    J = L[i][1].I
    T, _ = grade(parent(J[1]))
    V = T.(gens(J))
    JJ = ideal(T, V)
    B = quo(T, JJ)
    H = hilbert_polynomial(B[1])
    pa = pa - ZZ(coeff(H, 0))
  end
  return pa + 1
end
