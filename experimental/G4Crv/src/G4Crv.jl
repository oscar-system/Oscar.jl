################################################################################
#
#          CrvG4/CrvG4.jl : Non-hyperelliptic genus 4 curves
#
# (C) 2026
#
################################################################################

################################################################################
#
#  Types
#
################################################################################

mutable struct G4Crv{BaseRingType<:Field, RingType<:Ring, T<:FieldElem} <: AbsProjectiveScheme{BaseRingType, RingType}
  base_field::Ring
  conic::MPolyRingElem{T}
  cubic::MPolyRingElem{T}
  disc::T
  scheme::ProjectiveScheme{BaseRingType, RingType}

  function G4Crv(conic::MPolyRingElem{T}, cubic::MPolyRingElem{T}, check::Bool = true) where {T}
    R = base_ring(conic)

    @req total_degree(conic) == 2 && total_degree(cubic) == 3 "Nonhyperelliptic genus 4 curve is defined by the intersection of a conic and a cubic."
    @req is_homogeneous(conic) && is_homogeneous(cubic) "Input needs to consist of homogeneous polynomials."
    #Compute discriminant d
    #if d != 0 || check == false

      C = new{parent_type(T), MPolyQuoRing{MPolyDecRingElem{T, typeof(conic)}}, T}()

      C.conic = conic
      C.cubic = cubic
      #C.disc = d
      C.base_field = R


    # else
    #  error("Discriminant is zero")
    # end
    return C
  end
end

#=
mutable struct G4CrvPt{T}
  coordx::T
  coordy::T
  coordz::T
  coordw::T
  is_infinite::Bool
  parent::G4Crv{T}

  function G4CrvPt{T}(C::G4Crv{T}, coords::Vector{T}, check::Bool = true) where {T}
    K = base_field(C)
    P = new{T}()
    if check
      if !is_on_curve(C, coords)
        error("Point is not on the curve")
      end
    end

    P.parent = C
    if coords[4] == 0
      P.coordx = coords[1]
      P.coordy = coords[2]
      P.coordz = coords[3]
      P.coordw = coords[4]
      P.is_infinite = true
    else
      P.is_infinite = false

      #Don't have numerators, denominators and gcd over finite fields
      if T <: FieldElem

        scalar = inv(coords[4])

        P.coordx = coords[1]*scalar
        P.coordy = coords[2]*scalar
        P.coordz = coords[3]*scalar
        P.coordw = coords[4]*scalar
      else

        #Eliminate denominators
        D = prod(map(denominator, coords))
        map!(x -> D*x, coords)
        c = gcd(coords)
        map!(x -> divexact(x, c), coords)
        #Eliminate gcd
        
        P.coordx = coords[1]
        P.coordy = coords[2]
        P.coordz = coords[3]
        P.coordw = coords[4]
      end
    end
    return P
  end
end

function Base.getindex(P::G4CrvPt, i::Int)
  @req 1 <= i <= 4 "Index must be 1, 2 or 3"

  if i == 1
    return P.coordx
  elseif i == 2
    return P.coordy
  elseif i == 3
    return P.coordz
  elseif i == 4
    return P.coordw
  end
end
=#

################################################################################
#
#  Constructors for Non-hyperelliptic Genus 4 Curve
#
################################################################################

@doc raw"""
    g4_curve(conic::PolyRingElem, cubic::PolyRingElem; check::Bool = true) -> G4Crv

Return the non-hyperelliptic genus 4 curve defined by the intersection of the 
conic and the cubic.
"""
function g4_curve(conic::MPolyRingElem{T}, cubic::MPolyRingElem{T}; check::Bool = true) where T <: FieldElem
  @req parent(conic) == parent(cubic) "Conic and cubic need to have the same parent."
  @req total_degree(conic) == 2 && total_degree(cubic) == 3 "First argument needs to be of degree 2
   and second argument needs to be of degree 3."
  L = parent(conic)
  if length(gens(L)) == 4
    return G4Crv(conic, cubic, check)
  elseif length(gens(L)) == 3
    R = base_ring(conic)
    Rxyzw, (x, y, z, w) = polynomial_ring(R, [:x, :y, :z, :w])
    conic_hom = sum([mon*w^(2-total_degree(mon)) for mon in monomials(conic)];init = zero(Rxyzw))
    conic_cubic = sum([mon*w^(3-total_degree(mon)) for mon in monomials(conic)];init = zero(Rxyzw))
    return G4Crv(conic_hom, cubic_hom, check)
  else
    error("Can only construct a genus 4 curve as the intersection of a conic and a cubic in P^3")
  end

  
end


################################################################################
#
#  Underlying scheme
#
################################################################################


function underlying_scheme(C::G4Crv{T}) where T
  if !isdefined(C, :scheme)
    conic, cubic = equations(C)
    R = parent(conic)
    if !is_z_graded(R)
      R, _ = grade(R)
      conic = R(conic)
      cubic = R(cubic)
    end

    I = ideal(R, [conic, cubic])
    C_sch = projective_scheme(I)
    C.scheme = C_sch
  
  end

  return C.scheme
end


################################################################################
#
#  Field access
#
################################################################################

@doc raw"""
    base_field(C::G4Crv) -> Field

Return the base field over which `C` is defined.
"""
function base_field(C::G4Crv{T}) where T
  return C.base_field::parent_type(T)
end

################################################################################
#
#  Base Change
#
################################################################################

@doc raw"""
    base_change(K::Field, C::G4Crv) -> G4Curve

Return the base change of the conic curve $C$ over K if coercion is
possible.
"""
function base_change(K::Field, C::G4Crv)
  conic, cubic = equations(C)
  conic_new = change_coefficient_ring(K, conic)
  cubic_new = change_coefficient_ring(K, cubic)
  return hyperelliptic_curve(conic_new, cubic_new)
end


################################################################################
#
#  Equations
#
################################################################################

@doc raw"""
    equations(C::G4Crv) -> Poly, Poly

Return the conic and cubic defining the canonical embedding of the curve C.
"""
function equations(C::G4Crv)
  return C.conic, C.cubic
end


################################################################################
#
#  Points on Nonhyperelliptic Genus 4 Curves
#
################################################################################

function (C::G4Crv{T})(coords::Vector{S}; check::Bool = true) where {S, T}
  if !(3 <= length(coords) <= 4)
    error("Points need to be given in either affine coordinates (x, y, z) or projective coordinates (x, y, z, w).")
  end

  if length(coords) == 3
    push!(coords, 1)
  end
  if S === T
    parent(coords[1]) != base_field(C) &&
        error("Objects must be defined over same field.")
    return G4CrvPt{T}(C, coords, check)
  else
    return G4CrvPt{T}(C, map(base_field(C), coords)::Vector{T}, check)
  end
end


################################################################################
#
#  Test for inclusion
#
################################################################################
#=
@doc raw"""
    is_on_curve(C::G4Crv{T}, coords::Vector{T}) -> Bool

Return true if `coords` defines a point on $C$ and false otherwise. The array
`coords` must have length 4.
"""
function is_on_curve(C::G4Crv{T}, coords::Vector{T}) where T
  length(coords) != 4 && error("Array must be of length 4.")
  coords
  x, y, z, w = coords

  K = parent(x)

  if all(i -> i == zero(K), coords)
    error("(0 : 0 : 0) is not a point in projective space.")
  end

  conic, cubic = equations(C)
  if conic(coords...) == 0 && cubic(coords...) == 0
    return true
  else
    return false
  end
end
=#
################################################################################
#
#  Parent
#
################################################################################
#=
function parent(P::G4CrvPt)
  return P.parent
end
=#

################################################################################
#
#  ElemType
#
################################################################################

function elem_type(::Type{G4Crv{T}}) where T
  return G4CrvPt{T}
end


################################################################################
#
#  String I/O
#
################################################################################

function show(io::IO, C::G4Crv)
  conic, cubic = equations(C)
  print(io, "Nonhyperelliptic curve of genus 4 given by the intersection of \n $(conic) and \n $(cubic)")
end
#=
function show(io::IO, P::G4CrvPt)
   print(io, "Point  ($(P[1]) : $(P[2]) : $(P[3]) : $(P[4]))  of $(P.parent)")
end

@doc raw"""
    ==(P::G4CurvePoint, Q::G4CurvePoint) -> Bool

Return true if $P$ and $Q$ are equal and live over the same nonhyperelliptic
genus 4 curve $C$.
"""
function ==(P::G4CrvPt{T}, Q::G4CrvPt{T}) where T
  if parent(P) != parent(Q)
    return false
  end
  # Compare coordinates
  if P[1] == Q[1] && P[2] == Q[2] && P[3] == Q[3] && P[4] == Q[4]
    return true
  else
    return false
  end
end

function Base.hash(P::G4CrvPt, h::UInt)
  h = hash(parent(P), h)
  h = hash(P[1], h)
  h = hash(P[2], h)
  h = hash(P[3], h)
  h = hash(P[4], h)
  return h
end
=#
include("Invariants.jl")
include("Reconstruction.jl")
include("Minimization.jl")
