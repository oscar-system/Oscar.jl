
################################################################################
@doc raw"""
    ProjectivePlaneCurve <: AbsProjectiveCurve

A reduced curve in the projective plane.

# Examples
```jldoctest
julia> R, (x,y,z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

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
      is_equidimensional(X) || error("not equidimensional")
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



ProjectivePlaneCurve(I::MPolyIdeal{<:MPolyDecRingElem}; kwargs...) = ProjectivePlaneCurve(algebraic_set(I; kwargs...))

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
  if !isdefined(C, :defining_equation)
    m = minimal_generating_set(vanishing_ideal(C))
    C.defining_equation = only(m)
  end
  return C.defining_equation::T
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
    tangent_lines(C::ProjectivePlaneCurve{S}, P::AbsProjectiveRationalPoint) where S <: FieldElem

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
