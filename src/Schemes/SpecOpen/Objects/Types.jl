export SpecOpen

########################################################################
# Type for Zariski-open subsets of affine schemes                      #
########################################################################
@Markdown.doc """
    SpecOpen{SpecType, BRT} <: Scheme{BRT}

Zariski open subset ``U`` of an affine scheme ``X = Spec(R)``.
This stores a list of generators ``f₁,…,fᵣ`` of an ideal
``I`` defining the complement ``Z = X ∖ U``.
The scheme ``X`` is referred to as the *ambient scheme* and
the list ``f₁,…,fᵣ`` as the *generators* for ``U``.
"""
@attributes mutable struct SpecOpen{SpecType, BRT} <: Scheme{BRT}
  X::SpecType # the ambient scheme
  gens::Vector # a list of functions defining the complement of the open subset

  # fields used for caching
  name::String
  patches::Vector{AbsSpec}
  intersections::Dict{Tuple{Int, Int}, AbsSpec}
  complement::AbsSpec
  complement_ideal::Ideal
  ring_of_functions::Ring

  function SpecOpen(
      X::SpecType,
      f::Vector{RET};
      name::String="",
      check::Bool=true
    ) where {SpecType<:AbsSpec, RET<:RingElem}
    for a in f
      parent(a) == ambient_coordinate_ring(X) || error("element does not belong to the correct ring")
      if check
        !isempty(X) && iszero(OO(X)(a)) && error("generators must not be zero")
      end
    end
    U = new{SpecType, typeof(base_ring(X))}(X, f)
    U.intersections = Dict{Tuple{Int, Int}, AbsSpec}()
    length(name) > 0 && set_name!(U, name)
    return U
  end
  ### Conversion from PrincipalOpenSubsets
  function SpecOpen(U::PrincipalOpenSubset)
    X = ambient_scheme(U)
    h = complement_equation(U)
    V = new{typeof(X), typeof(base_ring(X))}(X, [lifted_numerator(h)])
    V.intersections = Dict{Tuple{Int, Int}, AbsSpec}()
    V.patches = [U]
    return V
  end
end

########################################################################
# Common type fo subsets of affine space                               #
########################################################################

SpecSubset = Union{SpecOpen,AbsSpec,PrincipalOpenSubset}
