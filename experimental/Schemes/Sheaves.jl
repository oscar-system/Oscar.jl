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
# Implementation of PreSheafOnScheme                                   #
########################################################################

### implementing the essential functionality
space(F::PreSheafOnScheme) = F.X
object_cache(F::PreSheafOnScheme) = F.obj_cache # an IdDict caching the values of F on open subsets
#restriction_cache(F::PreSheafOnScheme) = F.res_cache # Caching is now done via attributes in the objects
is_open_func(F::PreSheafOnScheme) = F.is_open_func
production_func(F::PreSheafOnScheme) = F.production_func
restriction_func(F::PreSheafOnScheme) = F.restriction_func

### Production and caching of the values of F on admissible open sets
function (F::PreSheafOnScheme{<:Any, OpenType, OutputType})(U::T; cached::Bool=true, check::Bool=true) where {OpenType, OutputType, T<:OpenType}
  #First look whether or not the asked for result has been computed before
  haskey(object_cache(F), U) && return (object_cache(F)[U])::OutputType 

  # Testing openness might be expensive, so it can be skipped
  check && is_open_func(F)(U, space(F)) || error("the given set is not open or admissible")
  G = production_func(F)(F, U)
  cached && (object_cache(F)[U] = G)
  return G::OutputType
end

### Production and caching of the restriction maps
@Markdown.doc """
    restriction_map(F::PreSheafOnScheme{<:Any, OpenType, OutputType, RestrictionType},
        U::Type1, V::Type2
      ) where {OpenType, OutputType, RestrictionType, Type1<:OpenType, Type2<:OpenType}

For a `F` produce (and cache) the restriction map `F(U) → F(V)`.
"""
function restriction_map(F::PreSheafOnScheme{<:Any, OpenType, OutputType, RestrictionType},
    U::Type1, V::Type2
  ) where {OpenType, OutputType, RestrictionType, Type1<:OpenType, Type2<:OpenType}
  # First, look up whether this restriction had already been asked for previously.
  inc = incoming_restrictions(F, F(V)) 
  !(inc == nothing) && haskey(inc, U) && return (inc[U])::RestrictionType

  # Check whether the given pair is even admissible.
  is_open_func(F)(V, U) || error("the second argument is not open in the first")

  # Hand the production of the restriction over to the internal method 
  rho = restriction_func(F)(F, U, V)

  # Cache the result in the attributes of F(V)
  inc isa IdDict{<:OpenType, <:RestrictionType} && (inc[U] = rho) # It is the restriction coming from U.
  return rho::RestrictionType
end

@Markdown.doc """
  add_incoming_restriction!(F::AbsPreSheaf{<:Any, OpenType, <:Any, RestrictionType}, 
    U::OpenType,
    V::OpenType,
    rho::RestrictionType
  ) where {OpenType, RestrictionType}

**Note:** This is a method for internal use! 

For an `AbsPreSheaf` `F`, a pair of open sets ``U ⊃ V`` and a manually computed 
morphism ``ρ : F(U) → F(V)``, this method stores the map `rho` in the internal caching system as 
the restriction map for `F` from `U` to `V`.
"""
function add_incoming_restriction!(F::AbsPreSheaf{<:Any, OpenType, OutputType, RestrictionType}, 
    U::OpenType,
    M::OutputType,
    rho::RestrictionType
  ) where {OpenType, OutputType, RestrictionType}
  # First, look up the incoming restriction maps for F(V).
  # This will create the dictionary, if necessary.
  incoming_res = incoming_restrictions(F, M)
  incoming_res == nothing && return F # This indicates that no 
  incoming_res::IdDict{<:OpenType, <:RestrictionType}
  incoming_res[U] = rho
  return F
end

@Markdown.doc """
    incoming_restrictions(F::AbsPreSheaf{<:Any, OpenType, OutputType, RestrictionType, M::OutputType) 

Supposing `M` is the value `M = F(U)` of some `AbsPreSheaf` `F` on an admissible open 
set `U`, return an `IdDict` whose keys `V` are those admissible open sets for `F` 
for which a restriction map `ρ : F(V) → F(U)` has already been computed and cached. 
The values of the dictionary are precisely those restriction maps for the respective keys.

**Note:** This 
"""
function incoming_restrictions(
    F::AbsPreSheaf{<:Any, OpenType, OutputType, RestrictionType},
    M::OutputType
  ) where {OpenType, OutputType, RestrictionType}
  hasfield(typeof(M), :__attrs) || return nothing # M has to be attributable to allow for caching!
  if !has_attribute(M, :incoming_restrictions)
    D = IdDict{OpenType, RestrictionType}()
    set_attribute!(M, :incoming_restrictions, D)
    return D
  end
  return get_attribute(M, :incoming_restrictions)::IdDict{OpenType, RestrictionType}
end

########################################################################
# The implementation of StructureSheafOfRings                         #
########################################################################

function _is_open_func_for_schemes(X::AbsCoveredScheme)
  ### Checks for open containment.
  #
  # We allow the following cases:
  #
  #  * U::PrincipalOpenSubset with one ancestor W in the basic charts of X
  #  * U::SimplifiedSpec with one ancestor W in the basic charts of X
  #  * U::PrincipalOpenSubset ⊂ V::PrincipalOpenSubset with ambient_scheme(U) === ambient_scheme(V) in the basic charts of X
  #  * U::PrincipalOpenSubset ⊂ V::PrincipalOpenSubset with ambient_scheme(U) != ambient_scheme(V) both in the basic charts of X
  #    and U and V contained in the glueing domains of their ambient schemes
  #  * U::AbsSpec ⊂ U::AbsSpec in the basic charts of X
  #  * U::AbsSpec ⊂ X for U in the basic charts
  #  * U::PrincipalOpenSubset ⊂ X with ambient_scheme(U) in the basic charts of X
  #  * W::SpecOpen ⊂ X with ambient_scheme(U) in the basic charts of X
  function is_open_func(
      U::Union{<:PrincipalOpenSubset, <:SimplifiedSpec}, 
      V::Union{<:PrincipalOpenSubset, <:SimplifiedSpec}
    )
    inc_U_flat = _flatten_open_subscheme(U, default_covering(X))
    inc_V_flat = _flatten_open_subscheme(V, default_covering(X))
    A = ambient_scheme(codomain(inc_U_flat))
    B = ambient_scheme(codomain(inc_V_flat))
    Udirect = codomain(inc_U_flat)
    Vdirect = codomain(inc_V_flat)
   # some_ancestor(W->any(WW->(WW===W), affine_charts(X)), U) || return false
   # some_ancestor(W->any(WW->(WW===W), affine_charts(X)), V) || return false
   # incU, dU = _find_chart(U, default_covering(X))
   # incV, dU = _find_chart(V, default_covering(X))
   # A = codomain(incU)
   # B = codomain(incV)
   # Udirect = (ambient_scheme(U) === A ? U : PrincipalOpenSubset(codomain(incU), dU))
   # Vdirect = (ambient_scheme(V) === B ? V : PrincipalOpenSubset(codomain(incV), dV))

    if A === B
      is_subset(Udirect, Vdirect) || return false
    else
      G = default_covering(X)[A, B] # Get the glueing
      f, g = glueing_morphisms(G)
      is_subset(Udirect, domain(f)) || return false
      gU = preimage(g, Udirect)
      is_subset(gU, Vdirect) || return false
    end
    return true
  end
  function is_open_func(
      U::Union{<:PrincipalOpenSubset, <:SimplifiedSpec}, 
      Y::AbsCoveredScheme
    )
    return Y === X && some_ancestor(W->(W in affine_charts(X)), U)
  end
  function is_open_func(U::AbsSpec, Y::AbsCoveredScheme)
    return Y === X && some_ancestor(W->any(WW->(WW===W), affine_charts(X)), U)
  end
  function is_open_func(Z::AbsCoveredScheme, Y::AbsCoveredScheme)
    return X === Y === Z
  end
  function is_open_func(U::AbsSpec, V::AbsSpec)
    U in affine_charts(X) || return false
    V in affine_charts(X) || return false
    G = default_covering(X)[U, V]
    return issubset(U, glueing_domains(G)[1])
  end
  function is_open_func(
      U::Union{<:PrincipalOpenSubset, <:SimplifiedSpec}, 
      V::AbsSpec
    )
    V in affine_charts(X) || return false
    inc_U_flat = _flatten_open_subscheme(U, default_covering(X))
    A = ambient_scheme(codomain(inc_U_flat))
    Udirect = codomain(inc_U_flat)
    W = ambient_scheme(Udirect)
    # some_ancestor(W->(W===V), U) && return true
    # incU, dU = _find_chart(U, default_covering(X))
    # W = codomain(incU)
    haskey(glueings(default_covering(X)), (W, V)) || return false # In this case, they are not glued
    G = default_covering(X)[W, V]
    #Udirect = (ambient_scheme(U) === W ? U : PrincipalOpenSubset(W, dU))
    return is_subset(Udirect, glueing_domains(G)[1])
  end
  function is_open_func(W::SpecOpen, Y::AbsCoveredScheme)
    return Y === X && ambient_scheme(W) in default_covering(X)
  end

  function is_open_func(W::SpecOpen, V::AbsSpec)
    V in default_covering(X) || return false
    ambient_scheme(W) === V && return true
    U = ambient_scheme(W)
    U in default_covering(X) || return false
    G = default_covering(X)[U, V]
    return is_subset(W, glueing_domains(G)[1])
  end
  function is_open_func(
      W::SpecOpen, 
      V::Union{<:PrincipalOpenSubset, <:SimplifiedSpec}
    )
    PW = ambient_scheme(W)
    inc_V_flat = _flatten_open_subscheme(V, default_covering(X))
    Vdirect = codomain(inc_V_flat)
    PV = ambient_scheme(Vdirect)
    #incV, dV = _find_chart(V, default_covering(X))
    #Vdirect = PrincipalOpenSubset(codomain(incV), dV)
    #PV = codomain(incV)
    PW in default_covering(X) || return false
    PV in default_covering(X) || return false
    if PW === PV
      return issubset(W, V)
      #return all(x->(issubset(x, V)), affine_patches(W))
    else
      haskey(glueings(default_covering(X)), (PW, PV)) || return false
      G = default_covering(X)[PW, PV]
      preV = preimage(glueing_morphisms(G)[1], Vdirect)
      return issubset(W, preV)
    end
  end
  function is_open_func(W::SpecOpen, V::SpecOpen)
    PW = ambient_scheme(W)
    PV = ambient_scheme(V)
    PW in default_covering(X) || return false
    PV in default_covering(X) || return false
    if PW === PV
      return issubset(W, V)
      #return all(x->(issubset(x, V)), affine_patches(W))
    else
      G = default_covering(X)[PW, PV]
      preV = preimage(glueing_morphisms(G)[1], V)
      return issubset(W, preV)
    end
  end
  function is_open_func(U::AbsSpec, W::SpecOpen)
    U in default_covering(X) || return false
    if U === ambient_scheme(W)
      # in this case W must be equal to U
      return issubset(W, U)
      #return one(OO(U)) in complement_ideal(W)
    else
      G = default_covering(X)[ambient_scheme(W), U]
      issubset(U, glueing_domains(G)[2]) || return false
      preU = preimage(glueing_morphisms(G)[1], U)
      return issubset(preU, W)
    end
  end
  function is_open_func(
      U::Union{<:PrincipalOpenSubset, <:SimplifiedSpec}, 
      W::SpecOpen
    )
    inc_U_flat = _flatten_open_subscheme(U, default_covering(X))
    A = ambient_scheme(codomain(inc_U_flat))
    Udirect = codomain(inc_U_flat)
    #some_ancestor(W->any(WW->(WW===W), affine_charts(X)), U)
    #incU, dU = _find_chart(U, default_covering(X))
    #U_direct = PrincipalOpenSubset(codomain(incU), dU)
    #ambient_scheme(U) in default_covering(X) || return false
    U_flat = codomain(inc_U_flat)
    PU = ambient_scheme(U_flat)
    if PU === ambient_scheme(W)
      # in this case W must be equal to U
      return issubset(W, U_flat)
      #return one(OO(U)) in complement_ideal(W)
    else
      G = default_covering(X)[ambient_scheme(W), PU]
      issubset(U_flat, glueing_domains(G)[2]) || return false
      preU = preimage(glueing_morphisms(G)[1], U_flat)
      return issubset(preU, W)
    end
  end

end

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
  return OO(U)(a)
end

########################################################################
# Sections of sheaves                                                  #
########################################################################

abstract type AbsPreSheafSection{SpaceType, 
                                 ParentType<:AbsPreSheaf, 
                                 OpenType, 
                                 ElemType
                                } end

### The interface for sections in presheaves

# Calling for a representative of the section on some open subset
@Markdown.doc """
    (v::AbsPreSheafSection{<:Any, <:AbsPreSheaf, OpenType})(U::OpenType) where {OpenType}

For a section ``v`` in a presheaf ``ℱ`` on ``X`` and an admissible 
open subset ``U ⊂ X`` return an element ``s ∈ ℱ(U)`` representing the section.
"""
function (v::AbsPreSheafSection{<:Any, <:AbsPreSheaf, OpenType})(U::OpenType) where {OpenType}
  return (underlying_section(v))(U)
end

### Further functionality

# The sheaf in which this is a section
parent(v::AbsPreSheafSection) = parent(underlying_section(v))

# The space on which the section is defined 
space(v::AbsPreSheafSection) = space(underlying_section(v))

# The dictionary with cached representatives of that section;
# for internal use only!
cache_dict(v::AbsPreSheafSection) = cache_dict(underlying_section(v))

# A function to produce a representative of that section on an 
# admissible open subset which is not yet cached; 
# for internal use only!
#
# The signature of this function should be 
#
#   production_func(w::PreSheafSection, U::OpenType) where {OpenType}
# 
# where OpenType is any admissible type for the open sets considered. 
# In the internal calls, the first argument w will always be the internal 
# representative of the section v itself. 
#
# This function should then take care of producing a representative 
# of v on U from the data already provided by w.
# 
production_func(v::AbsPreSheafSection) = production_func(underlying_section(v))


########################################################################
# Minimal concrete type for sections in presheaves                     #
########################################################################

# This type is for internal use in the user-facing concrete instances 
# of AbsPreSheafSection; to be returned via underlying_section().
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
    return new{typeof(X), typeof(F), OpenType, ElemType}(X, F, D, production_func)
  end
end

parent(v::PreSheafSection) = v.F
space(v::PreSheafSection) = v.X
cache_dict(v::PreSheafSection) = v.D
production_func(v::PreSheafSection) = v.production_func

function (v::PreSheafSection{<:Any, <:AbsPreSheaf, OpenType, ElemType})(U::OpenType) where {OpenType, ElemType}
  # First, look up whether this local representative has already been cached
  haskey(cache_dict(v), U) && return cache_dict(v)[U]
  # If not, delegate to the production function.
  # These can restrict from suitable bigger sets, assemble from compatible 
  # sections on refinements, extend trivially from other patches...
  s = production_func(v)(v, U)
  #cache_dict(v)[U] = s #TODO: Really cache all local representatives?
  return s::ElemType
end

