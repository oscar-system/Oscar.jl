export AbsPrePreSheaf
export space, restriction_map
export PreSheafOnScheme
export StructureSheafOfRings
export underlying_presheaf

########################################################################
# The AbsPreSheaf interface                                               #
########################################################################
@Markdown.doc """
    space(F::AbsPreSheaf) 

For a sheaf ``ℱ`` on a space ``X`` return ``X``.
"""
function space(F::AbsPreSheaf) 
  return space(underlying_presheaf(F))
end

@Markdown.doc """
    (F::AbsPreSheaf)(U; cached=true) 

For a sheaf ``ℱ`` on a space ``X`` and an (admissible) open set 
``U ⊂ X`` check whether ``U`` is open in ``X`` and return ``ℱ(U)``.
"""
function (F::AbsPreSheaf{<:Any, OpenType, OutputType})(U::T; cached::Bool=true) where {OpenType, OutputType, T<:OpenType}
  return (underlying_presheaf(F))(U, cached=cached)::OutputType
end

@Markdown.doc """
    restriction_map(F::AbsPreSheaf, U, V)

For a sheaf ``ℱ`` on a space ``X`` and an (admissible) pair of 
open sets ``U, V ⊂ X`` check whether ``U ⊂ V ⊂ X`` are open and 
return the restriction map ``ℱ(V) → ℱ(U)``.
"""
function restriction_map(F::AbsPreSheaf{<:Any, OpenType, OutputType, RestrictionType},
    U::Type1, V::Type2
  ) where {OpenType, OutputType, RestrictionType, Type1<:OpenType, Type2<:OpenType}
  return restriction_map(underlying_presheaf(F), U, V)::RestrictionType
end

# An alias for shorter notation
function (F::AbsPreSheaf{<:Any, OpenType, OutputType, RestrictionType})(
    U::Type1, V::Type2
  ) where {OpenType, OutputType, RestrictionType, Type1<:OpenType, Type2<:OpenType}
  return restriction_map(F, U, V)::RestrictionType
end

@Markdown.doc """
    is_open_func(F::AbsPreSheaf)

For a sheaf ``ℱ`` on a space ``X`` return a function `f` on two 
arguments such that `f(U, V)` returns `true` whenever ``U ⊂ V ⊂ X``
are open inclusions and `false` otherwise.

***Note:*** This function is expected to accept ``X`` as the second argument 
to check whether ``U ⊂ X`` is open in ``X``!
"""
function is_open_func(F::AbsPreSheaf)
  return is_open_func(underlying_presheaf(F))
end

########################################################################
# Implementation of PreSheafOnScheme                                      #
########################################################################

### implementing the essential functionality
space(F::PreSheafOnScheme) = F.X
object_cache(F::PreSheafOnScheme) = F.obj_cache
restriction_cache(F::PreSheafOnScheme) = F.res_cache
is_open_func(F::PreSheafOnScheme) = F.is_open_func
production_func(F::PreSheafOnScheme) = F.production_func
restriction_func(F::PreSheafOnScheme) = F.restriction_func

function (F::PreSheafOnScheme{<:Any, OpenType, OutputType})(U::T; cached::Bool=true, check::Bool=true) where {OpenType, OutputType, T<:OpenType}
  haskey(object_cache(F), U) && return (object_cache(F)[U])::OutputType #We can always look whether or not the asked for result has been computed before

  # Testing openness might be expensive, so it can be skipped
  check && is_open_func(F)(U, space(F)) || error("the given set is not open or admissible")
  G = production_func(F)(U, object_cache(F), restriction_cache(F))
  cached && (object_cache(F)[U] = G)
  return G::OutputType
end

function restriction_map(F::PreSheafOnScheme{<:Any, OpenType, OutputType, RestrictionType},
    U::Type1, V::Type2
  ) where {OpenType, OutputType, RestrictionType, Type1<:OpenType, Type2<:OpenType}
  haskey(restriction_cache(F), (U, V)) && return (restriction_cache(F)[(U, V)])::RestrictionType

  is_open_func(F)(V, U) || error("the second argument is not open in the first")
  FV = F(V)
  FU = F(U)
  rho = restriction_func(F)(U, V, object_cache(F), restriction_cache(F))
  restriction_cache(F)[(U, V)] = rho
  return rho::RestrictionType
end

########################################################################
# The implementation of StructureSheafOfRings                         #
########################################################################

underlying_presheaf(S::StructureSheafOfRings) = S.OO

@attr StructureSheafOfRings function OO(X::AbsCoveredScheme)
  return StructureSheafOfRings(X)
end

function Base.show(io::IO, R::StructureSheafOfRings)
  print(io, "𝒪_{$(space(R))}")
end
### Missing methods for compatibility of SimpleGlueings with Glueings
function restrict(
    a::Union{MPolyElem, MPolyQuoElem, 
             MPolyLocalizedRingElem, MPolyQuoLocalizedRingElem}, 
    U::PrincipalOpenSubset)
  parent(a) == OO(ambient_scheme(U)) || return OO(U)(a)
  return OO(U)(a, check=false)
end

function restrict(
    a::Union{MPolyElem, MPolyQuoElem, 
             MPolyLocalizedRingElem, MPolyQuoLocalizedRingElem}, 
    U::SpecOpen)
  parent(a) == OO(ambient_scheme(U)) || return OO(U)(a)
  return OO(U)(a, check=false)
end

########################################################################
# Sections of sheaves                                                  #
########################################################################

abstract type AbsPreSheafSection{SpaceType, 
                                 ParentType<:AbsPreSheaf, 
                                 OpenType, 
                                 ElemType
                                } end

@Markdown.doc """
    (v::AbsPreSheafSection{<:Any, <:AbsPreSheaf, OpenType})(U::OpenType) where {OpenType}

For a section ``v`` in a presheaf ``ℱ`` on ``X`` and an admissible 
open subset ``U ⊂ X`` return an element ``s ∈ ℱ(U)`` representing the section.
"""
function (v::AbsPreSheafSection{<:Any, <:AbsPreSheaf, OpenType})(U::OpenType) where {OpenType}
  error("method not implemented")
end

parent(v::AbsPreSheafSection) = parent(underlying_section(v))
space(v::AbsPreSheafSection) = space(underlying_section(v))
cache_dict(v::AbsPreSheafSection) = cache_dict(underlying_section(v))
production_func(v::AbsPreSheafSection) = production_func(underlying_section(v))

mutable struct PreSheafSection{SpaceType, 
                               ParentType<:AbsPreSheaf, 
                               OpenType,
                               ElemType
                              } <: AbsPreSheafSection{SpaceType, 
                                                      ParentType, 
                                                      OpenType,
                                                      ElemType
                                                     }
  X::SpaceType
  F::ParentType
  D::IdDict{OpenType, ElemType}
  production_func

  function PreSheafSection(F::AbsPreSheaf, 
      D::IdDict{OpenType, ElemType},
      production_func
    ) where {OpenType, ElemType}
    X = space(F)
    return new{typeof(X), typeof(F), OpenType, ElemType}(X, F, D)
  end
end

parent(v::PreSheafSection) = v,F
space(v::PreSheafSection) = v.X
cache_dict(v::PreSheafSection) = v.D
production_func(v::PreSheafSection) = v.production_func

function (v::PreSheafSection{<:Any, <:AbsPreSheaf, OpenType, ElemType})(U::OpenType) where {OpenType, ElemType}
  # First, look up whether this local representative has already been cached
  haskey(cache_dict(v)) && return cache_dict(v)[U]
  # If not, delegate to the production function.
  # These can restrict from suitable bigger sets, assemble from compatible 
  # sections on refinements, extend trivially from other patches...
  s = production_func(v)(U, cache_dict(v)) 
  cache_dict(v)[U] = s
  return s::ElemType
end

function _production_func_for_coherent_sheaves(F::AbsCoherentSheaf)
  X = scheme(F)
  # We assume the following.
  # The section is defined on all affine charts of X.
  # Every admissible open subset is a PrincipalOpenSubset of 
  # some affine chart.
  function prod_func(U::AbsSpec, D::IdDict)
    U in affine_charts(X) || error("open set is not an affine chart of the scheme")
    error("section is not defined in this chart")
  end
  function prod_func(U::PrincipalOpenSubset, D::IdDict)
    V = ambient_scheme(U)
    t = haskey(D, V) ? D[V] : production_func(V, D)
    # Per default we do not cache local sections 
    return F(V, U)(t)
  end
end

mutable struct CoherentSheafSection{SpaceType, ParentType<:AbsPreSheaf, 
                                    OpenType, ElemType
                                   } <: AbsPreSheafSection{SpaceType, ParentType, 
                                                           OpenType, ElemType
                                                          }
  s::PreSheafSection{SpaceType, ParentType, OpenType, ElemType}

  function CoherentSheafSection(
      F::AbsCoherentSheaf,
      D::Dict{AbsSpec, ModuleFPElem}
    )
    X = scheme(F)
    all(U->haskey(D, U), affine_charts(X)) || error("section must be defined on all affine charts")
    s = PreSheafSection(F, D, _production_func_for_coherent_sheaves(F))
    return new{typeof(X), typeof(F), AbsSpec, ModuleFPElem}(s)
  end
end

underlying_section(s::CoherentSheafSection) = s.s

