export AbsSheaf
export space, restriction_map
export SheafOnScheme
export SheafOO

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
# A minimal implementation of the sheaf interface on an affine scheme  #
########################################################################
@attributes mutable struct SheafOnScheme{SpaceType, OpenType, OutputType, RestrictionType, 
                                       IsOpenFuncType, ProductionFuncType,
                                       RestrictionFuncType
                                      } <: AbsSheaf{
                                       SpaceType, OpenType, 
                                       OutputType, RestrictionType
                                      }
  X::SpaceType

  # caches
  obj_cache::IdDict{<:OpenType, <:OutputType} # To cache values that have already been computed
  res_cache::IdDict{<:Tuple{<:OpenType, <:OpenType}, <:RestrictionType} # To cache already computed restrictions

  # production functions for new objects
  is_open_func::IsOpenFuncType # To check whether one set is open in the other
  production_func::ProductionFuncType # To produce ℱ(U) for U ⊂ X
  restriction_func::RestrictionFuncType

  function SheafOnScheme(X::Scheme, production_func::Any, restriction_func::Any;
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
# The structure sheaf of affine schemes                                #
########################################################################
@attributes mutable struct SheafOO{SpaceType, OpenType, OutputType,
                                          RestrictionType, ProductionFuncType,
                                          RestrictionFuncType,
                                          SheafType
                                         } <: AbsSheaf{
                                          SpaceType, OpenType, 
                                          OutputType, RestrictionType
                                         }
  OO::SheafType

  ### Structure sheaf on affine schemes
  function SheafOO(X::AbsSpec)
    function is_open_func(U::AbsSpec, V::AbsSpec)
      return is_subset(V, X) && is_open_embedding(U, V) # Note the restriction to subsets of X
    end
    function production_func(U::AbsSpec)
      return OO(U)
    end
    function restriction_func(V::AbsSpec, U::AbsSpec)
      return hom(OO(V), OO(U), gens(OO(U)), check=false) # check=false assures quicker computation
    end

    R = SheafOnScheme(X, production_func, restriction_func, 
                    OpenType=AbsSpec, OutputType=Ring, 
                    RestrictionType=Hecke.Map,
                    is_open_func=is_open_func
                   )
    return new{typeof(X), AbsSpec, Ring, Hecke.Map, 
               typeof(production_func), typeof(restriction_func), 
               typeof(R)}(R)
  end

  ### Structure sheaf on covered schemes
  function SheafOO(X::AbsCoveredScheme)

    ### Checks for open containment. 
    #
    # We allow the following cases:
    #
    #  * U::PrincipalOpenSubset in W===ambient_scheme(U) in the basic charts of X
    #  * U::PrincipalOpenSubset ⊂ V::PrincipalOpenSubset with ambient_scheme(U) === ambient_scheme(V) in the basic charts of X
    #  * U::PrincipalOpenSubset ⊂ V::PrincipalOpenSubset with ambient_scheme(U) != ambient_scheme(V) both in the basic charts of X
    #    and U and V contained in the glueing domains of their ambient schemes
    #  * U::AbsSpec ⊂ U::AbsSpec in the basic charts of X
    #  * U::AbsSpec ⊂ X for U in the basic charts
    #  * U::PrincipalOpenSubset ⊂ X with ambient_scheme(U) in the basic charts of X
    function is_open_func(U::PrincipalOpenSubset, V::PrincipalOpenSubset)
      C = default_covering(X)
      A = ambient_scheme(U) 
      A in C || return false
      B = ambient_scheme(V) 
      B in C || return false
      if A == B
        is_subset(U, V) || return false
      else
        G = C[A, B] # Get the glueing
        f, g = glueing_morphisms(G)
        is_subset(U, domain(f)) || return false
        is_subset(V, domain(g)) || return false
        gU = preimage(g, U)
        is_subset(gU, V) || return false
      end
      return true
    end
    function is_open_func(U::PrincipalOpenSubset, Y::AbsCoveredScheme)
      return Y === X && ambient_scheme(U) in default_covering(X)
    end
    function is_open_func(U::AbsSpec, Y::AbsCoveredScheme)
      return Y === X && U in default_covering(X)
    end
    function is_open_func(Z::AbsCoveredScheme, Y::AbsCoveredScheme)
      return X === Y === Z
    end
    function is_open_func(U::AbsSpec, V::AbsSpec)
      return U === V && U in default_covering(X)
    end
    function is_open_func(U::PrincipalOpenSubset, V::AbsSpec)
      V in default_covering(X) || return false
      ambient_scheme(U) === V && return true
      W = ambient_scheme(U)
      W in default_covering(X) || return false
      G = default_covering(X)[W, V]
      return is_subset(U, glueing_domains(G)[1])
    end


    function production_func(U::AbsSpec)
      return OO(U)
    end

    function restriction_func(V::AbsSpec, U::AbsSpec)
      X === U || error("schemes must be the same")
      return identity_map(OO(X))
    end
    function restriction_func(V::AbsSpec, U::PrincipalOpenSubset)
      if ambient_scheme(U) === V
        return hom(OO(V), OO(U), gens(OO(U)), check=false)
      else
        W = ambient_scheme(U)
        G = default_covering(X)[V, W]
        f, g = glueing_morphisms(G)
        function rho_func(x::RingElem)
          parent(x) == OO(V) || error("element does not belong to the correct domain")
          return restrict(pullback(g)(x), U) # should probably be tuned to avoid checks. 
        end
        return MapFromFunc(rho_func, OO(V), OO(U))
      end
      error("arguments are not valid")
    end
    function restriction_func(V::PrincipalOpenSubset, U::PrincipalOpenSubset)
      A = ambient_scheme(V)
      if A === ambient_scheme(U)
        return hom(OO(V), OO(U), gens(OO(U)), check=false)
      else 
        B = ambient_scheme(U)
        G = default_covering(X)[A, B]
        f, g = glueing_morphisms(G)
        gV = preimage(g, V)
        gres = restrict(g, gV, V, check=false)
        h = inclusion_morphism(U, gV)
        function rho_func(x::RingElem)
          parent(x) == OO(V) || error("input not valid")
          return pullback(h)(pullback(gres)(x))
        end
        return MapFromFunc(rho_func, OO(V), OO(U))
      end
      error("arguments are invalid")
    end

    R = SheafOnScheme(X, production_func, restriction_func, 
                    OpenType=AbsSpec, OutputType=Ring, 
                    RestrictionType=Hecke.Map,
                    is_open_func=is_open_func
                   )
    return new{typeof(X), AbsSpec, Ring, Hecke.Map, 
               typeof(production_func), typeof(restriction_func), 
               typeof(R)}(R)
  end
end

underlying_sheaf(S::SheafOO) = S.OO


