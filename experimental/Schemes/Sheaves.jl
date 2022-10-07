export AbsSheaf
export space, restriction_map
export SheafOnSpec
@Markdown.doc """
    AbsSheaf{SpaceType, OpenType, OutputType, RestrictionType}

Abstract type for a sheaf ℱ on a space X.

 * `SpaceType` is a parameter for the type of the space ``X`` on which ``ℱ`` is defined. 

 * `OpenType` is a type (most probably abstract!) for the open sets ``U ⊂ X`` which are admissible as input for ``ℱ(U)``.

 * `OutputType` is a type (most probably abstract!) for the values that ``ℱ`` takes on admissible open sets ``U``.

 * `RestrictionType` is a parameter for the type of the restriction maps ``ℱ(V) → ℱ(U)`` for ``U ⊂ V ⊂ X`` open.
"""
abstract type AbsSheaf{SpaceType, OpenType, OutputType, RestrictionType} end

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
                    U::Type1, V::Type1
                   ) where {
                     OpenType, OutputType, 
                     Type1<:OpenType, Type2<:OpenType
                   }

For a sheaf ``ℱ`` on a space ``X`` and an (admissible) pair of 
open sets ``U, V ⊂ X`` check whether ``U ⊂ V ⊂ X`` are open and 
return the restriction map ``ℱ(V) → ℱ(U)``.
"""
function restriction_map(F::AbsSheaf{<:Any, OpenType, OutputType, RestrictionType},
    U::Type1, V::Type1
  ) where {OpenType, OutputType, RestrictionType, Type1<:OpenType, Type2<:OpenType}
  return restriction_map(underlying_sheaf(F), U, V)::RestrictionType
end

# An alias for shorter notation
function (F::AbsSheaf{<:Any, OpenType, OutputType, RestrictionType})(
    U::Type1, V::Type1
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
# A minimal implementation of the sheaf interface on an affine scheme  #
########################################################################
@attributes mutable struct SheafOnSpec{SpaceType, OpenType, OutputType, RestrictionType, 
                                       IsOpenFuncType, ProductionFuncType,
                                       RestrictionFuncType
                                      } <: AbsSheaf{SpaceType, OpenType, OutputType, RestrictionType}
  X::SpaceType

  # caches
  obj_cache::IdDict{<:OpenType, <:OutputType} # To cache values that have already been computed
  res_cache::IdDict{<:Tuple{<:OpenType, <:OpenType}, <:RestrictionType} # To cache already computed restrictions

  # production functions for new objects
  is_open_func::IsOpenFuncType # To check whether one set is open in the other
  production_func::ProductionFuncType # To produce ℱ(U) for U ⊂ X
  restriction_func::RestrictionFuncType

  function SheafOnSpec(X::AbsSpec, production_func::Any, restriction_func::Any;
      OpenType=AbsSpec, OutputType=Any, RestrictionType=Any,
      is_open_func::Any=is_open_embedding
    )
    return new{typeof(X), OpenType, OutputType, RestrictionType, 
               typeof(is_open_func), typeof(production_func), typeof(restriction_func)
              }(X, IdDict{OpenType, OutputType}(), 
                IdDict{Tuple{OpenType, OpenType}, RestrictionType}(),
                is_open_func, production_func, restriction_func
               )
  end
end

### implementing the essential functionality
space(F::SheafOnSpec) = F.X
object_cache(F::SheafOnSpec) = F.obj_cache
restriction_cache(F::SheafOnSpec) = F.res_cache
is_open_func(F::SheafOnSpec) = F.is_open_func
production_func(F::SheafOnSpec) = F.production_func
restriction_func(F::SheafOnSpec) = F.restriction_func

function (F::SheafOnSpec{<:Any, OpenType, OutputType})(U::T) where {OpenType, OutputType, T<:OpenType}
  haskey(object_cache(F), U) && return object_cache(F)[U]

  is_open_func(F)(U, space(F)) || error("the given set is not open or admissible")
  G = production_func(F)(U)
  object_cache(F)[U] = G
  return G
end

function restriction_map(F::SheafOnSpec{<:Any, OpenType, OutputType, RestrictionType},
    U::Type1, V::Type1
  ) where {OpenType, OutputType, RestrictionType, Type1<:OpenType, Type2<:OpenType}
  haskey(restriction_cache(F), (U, V)) && return object_cache(F)[(U, V)]

  is_open_func(F)(V, U) || error("the second argument is not open in the first")
  rho = restriction_func(F)(U, V)
  restriction_cache(F)[(U, V)] = rho
  return rho
end

