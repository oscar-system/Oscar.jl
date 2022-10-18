export AbsSheaf
export space, restriction_map
export SheafOnScheme
export RingOfRegularFunctions

########################################################################
# The AbsSheaf interface                                               #
########################################################################
@Markdown.doc """
    space(F::AbsSheaf) 

For a sheaf ``ℱ`` on a space ``X`` return ``X``.
"""
function space(F::AbsSheaf) 
  return space(underlying_sheaf(F))
end

@Markdown.doc """
    (F::AbsSheaf{<:Any, OpenType})(U::T) where {T<:OpenType}

For a sheaf ``ℱ`` on a space ``X`` and an (admissible) open set 
``U ⊂ X`` check whether ``U`` is open in ``X`` and return ``ℱ(U)``.
"""
function (F::AbsSheaf{<:Any, OpenType, OutputType})(U::T) where {OpenType, OutputType, T<:OpenType}
  return (underlying_sheaf(F))(U)::OutputType
end

@Markdown.doc """
    restriction_map(F::AbsSheaf{<:Any, OpenType, OutputType},
                    U::Type1, V::Type2
                   ) where {
                     OpenType, OutputType, 
                     Type1<:OpenType, Type2<:OpenType
                   }

For a sheaf ``ℱ`` on a space ``X`` and an (admissible) pair of 
open sets ``U, V ⊂ X`` check whether ``U ⊂ V ⊂ X`` are open and 
return the restriction map ``ℱ(V) → ℱ(U)``.
"""
function restriction_map(F::AbsSheaf{<:Any, OpenType, OutputType, RestrictionType},
    U::Type1, V::Type2
  ) where {OpenType, OutputType, RestrictionType, Type1<:OpenType, Type2<:OpenType}
  return restriction_map(underlying_sheaf(F), U, V)::RestrictionType
end

# An alias for shorter notation
function (F::AbsSheaf{<:Any, OpenType, OutputType, RestrictionType})(
    U::Type1, V::Type2
  ) where {OpenType, OutputType, RestrictionType, Type1<:OpenType, Type2<:OpenType}
  return restriction_map(F, U, V)::RestrictionType
end

@Markdown.doc """
    is_open_func(F::AbsSheaf)

For a sheaf ``ℱ`` on a space ``X`` return a function `f` on two 
arguments such that `f(U, V)` returns `true` whenever ``U ⊂ V ⊂ X``
are open inclusions and `false` otherwise.

***Note:*** This function must accept ``X`` as the second argument 
to check whether ``U ⊂ X`` is open in ``X``!
"""
function is_open_func(F::AbsSheaf)
  return is_open_func(underlying_sheaf(F))
end

########################################################################
# Implementation of SheafOnScheme                                      #
########################################################################

### implementing the essential functionality
space(F::SheafOnScheme) = F.X
object_cache(F::SheafOnScheme) = F.obj_cache
restriction_cache(F::SheafOnScheme) = F.res_cache
is_open_func(F::SheafOnScheme) = F.is_open_func
production_func(F::SheafOnScheme) = F.production_func
restriction_func(F::SheafOnScheme) = F.restriction_func

function (F::SheafOnScheme{<:Any, OpenType, OutputType})(U::T) where {OpenType, OutputType, T<:OpenType}
  haskey(object_cache(F), U) && return (object_cache(F)[U])::OutputType

  is_open_func(F)(U, space(F)) || error("the given set is not open or admissible")
  G = production_func(F)(U)
  object_cache(F)[U] = G
  return G::OutputType
end

function restriction_map(F::SheafOnScheme{<:Any, OpenType, OutputType, RestrictionType},
    U::Type1, V::Type2
  ) where {OpenType, OutputType, RestrictionType, Type1<:OpenType, Type2<:OpenType}
  haskey(restriction_cache(F), (U, V)) && return (restriction_cache(F)[(U, V)])::RestrictionType

  is_open_func(F)(V, U) || error("the second argument is not open in the first")
  rho = restriction_func(F)(U, V)
  restriction_cache(F)[(U, V)] = rho
  return rho::RestrictionType
end

########################################################################
# The implementation of RingOfRegularFunctions                         #
########################################################################

underlying_sheaf(S::RingOfRegularFunctions) = S.OO


