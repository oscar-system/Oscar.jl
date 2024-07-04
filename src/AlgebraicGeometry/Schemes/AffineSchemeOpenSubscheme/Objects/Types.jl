
########################################################################
# Type for Zariski-open subsets of affine schemes                      #
########################################################################
@doc raw"""
    AffineSchemeOpenSubscheme{AffineSchemeType, BRT} <: Scheme{BRT}

Zariski open subset ``U`` of an affine scheme ``X = Spec(R)``.
This stores a list of generators ``f₁,…,fᵣ`` of an ideal
``I`` defining the complement ``Z = X ∖ U``.
The scheme ``X`` is referred to as the *ambient scheme* and
the list ``f₁,…,fᵣ`` as the *generators* for ``U``.
"""
@attributes mutable struct AffineSchemeOpenSubscheme{AffineSchemeType, BRT} <: Scheme{BRT}
  X::AffineSchemeType # the ambient scheme
  gens::Vector # a list of functions defining the complement of the open subset

  # fields used for caching
  name::String
  patches::Vector{AbsAffineScheme}
  intersections::Dict{Tuple{Int, Int}, AbsAffineScheme}
  complement::AbsAffineScheme
  complement_ideal::Ideal
  ring_of_functions::Ring

  function AffineSchemeOpenSubscheme(
      X::AffineSchemeType,
      f::Vector{RET};
      name::String="",
      check::Bool=true
    ) where {AffineSchemeType<:AbsAffineScheme, RET<:RingElem}
    for a in f
      parent(a) == ambient_coordinate_ring(X) || error("element does not belong to the correct ring")
      @check !(!isempty(X) && iszero(OO(X)(a))) "generators must not be zero"
    end
    U = new{AffineSchemeType, typeof(base_ring(X))}(X, f)
    U.intersections = Dict{Tuple{Int, Int}, AbsAffineScheme}()
    length(name) > 0 && set_name!(U, name)
    return U
  end
  ### Conversion from PrincipalOpenSubsets
  function AffineSchemeOpenSubscheme(U::PrincipalOpenSubset)
    X = ambient_scheme(U)
    h = complement_equation(U)
    V = new{typeof(X), typeof(base_ring(X))}(X, [lifted_numerator(h)])
    V.intersections = Dict{Tuple{Int, Int}, AbsAffineScheme}()
    V.patches = [U]
    return V
  end
end

########################################################################
# Common type fo subsets of affine space                               #
########################################################################

const AffineSchemeSubset = Union{AffineSchemeOpenSubscheme,AbsAffineScheme,PrincipalOpenSubset}
