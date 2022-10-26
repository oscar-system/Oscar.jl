export AbsPrePreSheaf
export space, restriction_map
export PreSheafOnScheme
export StructureSheafOfRings

########################################################################
# The AbsPreSheaf interface                                               #
########################################################################
@Markdown.doc """
    space(F::AbsPreSheaf) 

For a sheaf ``ℱ`` on a space ``X`` return ``X``.
"""
function space(F::AbsPreSheaf) 
  return space(underlying_sheaf(F))
end

@Markdown.doc """
    (F::AbsPreSheaf)(U; cached=true) 

For a sheaf ``ℱ`` on a space ``X`` and an (admissible) open set 
``U ⊂ X`` check whether ``U`` is open in ``X`` and return ``ℱ(U)``.
"""
function (F::AbsPreSheaf{<:Any, OpenType, OutputType})(U::T; cached::Bool=true) where {OpenType, OutputType, T<:OpenType}
  return (underlying_sheaf(F))(U, cached=cached)::OutputType
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
  return restriction_map(underlying_sheaf(F), U, V)::RestrictionType
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
  return is_open_func(underlying_sheaf(F))
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
  G = production_func(F)(U)
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
  rho = restriction_func(F)(U, FU, V, FV)
  restriction_cache(F)[(U, V)] = rho
  return rho::RestrictionType
end

########################################################################
# The implementation of StructureSheafOfRings                         #
########################################################################

underlying_sheaf(S::StructureSheafOfRings) = S.OO

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
  parent(a) == OO(ambient(U)) || return OO(U)(a)
  return OO(U)(a, check=false)
end


