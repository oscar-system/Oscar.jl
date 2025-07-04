
@doc raw"""
    AffinePlaneCurve{BaseField<:Field, RingType<:Ring} <: AbsAffineCurve{BaseField, RingType}

Type for reduced affine plane curves.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> F = y^3*x^6 - y^6*x^2;

julia> C = plane_curve(F)
Affine plane curve
  defined by 0 = x^5*y - x*y^4

```
"""
@attributes mutable struct AffinePlaneCurve{BaseField<:Field, RingType<:Ring} <: AbsAffineCurve{BaseField, RingType}
  X::AbsAffineAlgebraicSet{BaseField, RingType}
  defining_equation::RingElem  # reduced defining equation

  function AffinePlaneCurve(X::AbsAffineAlgebraicSet{S,T}; check::Bool=true) where {S,T}
    @check begin
      ngens(ambient_coordinate_ring(X)) == 2 || error("$(X) is not contained in the affine plane")
      dim(X) == 1 || error("wrong dimension, not a curve")
      is_equidimensional(X) || error("not a curve")
    end
    new{S,T}(X)
  end
  function AffinePlaneCurve(eqn::E; check::Bool=true) where {E<:MPolyRingElem}
    @check begin
      ngens(parent(eqn)) == 2 || error("number of variables must be 2")
      (isone(eqn) || iszero(eqn)) && error("the equation must not be trivial")
      true
    end
    eqn = prod([i[1] for i in factor_squarefree(eqn)], init=one(parent(eqn)))
    X = algebraic_set(eqn; is_radical=true, check=check)
    C = new{base_ring_type(X),ring_type(X)}(X)
    C.defining_equation = eqn
    return C
  end
end

AffinePlaneCurve(X::AbsAffineScheme; args...) = AffinePlaneCurve(AffineAlgebraicSet(X;args...);args...)

# Functions to comply with the AbsAlgebraicSet interface

underlying_scheme(C::AffinePlaneCurve) = C.X

fat_scheme(C::AffinePlaneCurve) = fat_scheme(underlying_scheme(C))

dim(::AbsAffineCurve) = 1

plane_curve(eq::MPolyRingElem) = AffinePlaneCurve(algebraic_set(eq))

@doc raw"""
    defining_equation(C::AffinePlaneCurve)

Return the defining equation of `C`.
"""
function defining_equation(C::AffinePlaneCurve{S, MPolyQuoRing{E}}) where {S, E}
  if isdefined(C,:defining_equation)
    return C.defining_equation::E
  end
  J = vanishing_ideal(C)
  g = small_generating_set(J)
  @assert length(g) == 1
  C.defining_equation = g[1]
  return g[1]::E
end

function Base.hash(C::AffinePlaneCurve, h::UInt)
  lc = leading_coefficient(defining_equation(C))
  F = inv(lc)*defining_equation(C)
  return hash(F, h)
end


function Base.show(io::IO, ::MIME"text/plain", C::AffinePlaneCurve)
  io = pretty(io)
  println(io, "Affine plane curve")
  print(io, Indent(), "defined by 0 = ", defining_equation(C), Dedent())
end

# one line printing is inherited from AbsAlgebraicSet

function union(X::AffinePlaneCurve, Y::AffinePlaneCurve)
  return AffinePlaneCurve(union(underlying_scheme(X),underlying_scheme(Y)))
end

################################################################################
# gives the common components of two affine plane curves

@doc raw"""
    common_components(C::AffinePlaneCurve, D::AffinePlaneCurve)

Return the affine plane curve consisting of the common components of `C` and `D`,
or an empty vector if they do not have a common component. This
component can be reducible.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> C = plane_curve(x*(x+y)*(x^2 + x + 1));

julia> D = plane_curve(x*(x+y)*(x-y));

julia> common_components(C, D)
1-element Vector{AffinePlaneCurve{QQField, MPolyQuoRing{QQMPolyRingElem}}}:
 scheme(x^2 + x*y)

```
"""
function common_components(C::S, D::S) where {S<:AffinePlaneCurve}
  G = gcd(defining_equation(C), defining_equation(D))
  if isone(G)
     return Vector{S}()
  else
     return S[plane_curve(G)]
  end
end

function irreducible_components(C::AffinePlaneCurve)
  return [AffinePlaneCurve(i[1]) for i in factor(defining_equation(C))]
end

################################################################################
# Helping function:
# change of variable to send a point at the origin.

function curve_map_point_origin(C::AffinePlaneCurve, P::AbsAffineRationalPoint)
  F = defining_equation(C)
  R = parent(F)
  V = gens(R)
  G = evaluate(F, [V[1] + P[1], V[2] + P[2]])
  return AffinePlaneCurve(G)
end

################################################################################
# Helping function:
# used to order elements to find the multiplicity.

function _sort_helper_multiplicity(a::FinGenAbGroupElem)
  return a.coeff[1, 1]
end

################################################################################
# compute the multiplicity of the affine plane curve C at the point P

@doc raw"""
    multiplicity(C::AffinePlaneCurve, P::AbsAffineRationalPoint)

Return the multiplicity of `C` at `P`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> C = plane_curve(x^2*(x+y)*(y^3-x^2));

julia> P = C([2,-2])
Rational point
  of scheme(-x^4 - x^3*y + x^2*y^3 + x*y^4)
with coordinates (2, -2)

julia> multiplicity(C, P)
1

```
"""
function multiplicity(C::AffinePlaneCurve, P::AbsAffineRationalPoint)
  P in C || P in ambient_space(C) || error("The point needs to be in a two dimensional space")
  D = curve_map_point_origin(C, P)
  G = defining_equation(D)
  R = parent(G)
  A, _ = grade(R)
  HC = homogeneous_components(A(G))
  L = collect(keys(HC))
  M = sort(L, by=_sort_helper_multiplicity)
  return M[1].coeff[1, 1]
end

################################################################################
# compute the set of tangent lines of the affine plane curve C at the point P
# (linear factors of the homogeneous part of lowest degree of the equation).

@doc raw"""
    tangent_lines(C::AffinePlaneCurve, P::AbsAffineRationalPoint)

Return the tangent lines at `P` to `C` with their multiplicity.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> C = plane_curve(x^2*(x+y)*(y^3-x^2));

julia> P = C([0, 0])
Rational point
  of scheme(-x^4 - x^3*y + x^2*y^3 + x*y^4)
with coordinates (0, 0)

julia> tangent_lines(C, P)
Dict{AffinePlaneCurve{QQField, MPolyQuoRing{QQMPolyRingElem}}, Int64} with 2 entries:
  scheme(x)     => 3
  scheme(x + y) => 1

```
"""
function tangent_lines(C::AffinePlaneCurve, P::AbsAffineRationalPoint)
  P in C || error("The point needs to be in C")
  D = curve_map_point_origin(C, P)
  G = defining_equation(D)
  R = parent(G)
  V = gens(R)
  A, _ = grade(R)
  HC = homogeneous_components(A(G))
  L = collect(keys(HC))
  M = sort(L, by=_sort_helper_multiplicity)
  Gm = HC[M[1]]
  Z = factor(Gm.f)
  D = Dict{typeof(C), Int}()
  X = V[1] - P[1]
  Y = V[2] - P[2]
  for p in keys(Z.fac)
     if total_degree(p) == 1
        push!(D, AffinePlaneCurve(evaluate(p, [X, Y])) => Z.fac[p])
     end
  end
  return D
end

################################################################################
# Compute the intersection multiplicity of the two affine plane curves at the
# given point.

@doc raw"""
    intersection_multiplicity(C::AffinePlaneCurve, D::AffinePlaneCurve, P::AbsAffineRationalPoint)

Return the intersection multiplicity of `C` and `D` at `P`.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y]);

julia> C = plane_curve((x^2+y^2)*(x^2 + y^2 + 2*y))
Affine plane curve
  defined by 0 = x^4 + 2*x^2*y^2 + 2*x^2*y + y^4 + 2*y^3

julia> D = plane_curve((x^2+y^2)*(y^3*x^6 - y^6*x^2))
Affine plane curve
  defined by 0 = x^7*y + x^5*y^3 - x^3*y^4 - x*y^6

julia> Q = D([0, -2]);

julia> intersection_multiplicity(C, D, Q)
1

```
"""
function intersection_multiplicity(C::AffinePlaneCurve, D::AffinePlaneCurve, P::AbsAffineRationalPoint)
  P in ambient_space(C) || error("The point needs to be in the same affine plane as C")
  P in C && P in D || return 0
  S,_ = stalk(intersect(C,D), P)
  return vector_space_dim(S)
end

################################################################################
# Check if the two curves intersect transversally at the given point.

@doc raw"""
    is_transverse_intersection(C::AffinePlaneCurve, D::AffinePlaneCurve, P::AbsAffineRationalPoint)

Return `true` if `C` and `D` intersect transversally at `P` and `false` otherwise.

# Examples
```jldoctest
julia> R, (x, y) = polynomial_ring(QQ, [:x, :y])
(Multivariate polynomial ring in 2 variables over QQ, QQMPolyRingElem[x, y])

julia> C = plane_curve(x*(x+y))
Affine plane curve
  defined by 0 = x^2 + x*y

julia> D = plane_curve((x-y)*(x-2))
Affine plane curve
  defined by 0 = x^2 - x*y - 2*x + 2*y

julia> P = C([QQ(0), QQ(0)])
Rational point
  of scheme(x^2 + x*y)
with coordinates (0, 0)

julia> Q = C([QQ(2), QQ(-2)])
Rational point
  of scheme(x^2 + x*y)
with coordinates (2, -2)

julia> is_transverse_intersection(C, D, P)
false

julia> is_transverse_intersection(C, D, Q)
true

```
"""
function is_transverse_intersection(C::AffinePlaneCurve, D::AffinePlaneCurve, P::AbsAffineRationalPoint)
  P in C && P in D || return false
  any(P in i for i in common_components(C,D)) && return false
  return intersection_multiplicity(C, D, P) == 1
end

@doc raw"""
    projective_closure(C::AffinePlaneCurve) -> ProjectivePlaneCurve

Return the projective closure of `C`.
"""
function projective_closure(C::AffinePlaneCurve)
  F = defining_equation(C)
  H = homogenizer(parent(F), "z")
  return plane_curve(H(F))
end

geometric_genus(C::AffinePlaneCurve) = geometric_genus(projective_closure(C))


