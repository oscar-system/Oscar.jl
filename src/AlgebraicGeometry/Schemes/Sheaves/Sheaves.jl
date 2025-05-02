
########################################################################
# The AbsPreSheaf interface                                            
#
# Architecture for an `AbsPreSheaf` `F` on a scheme `X`
#
# User facing calls are `F(U)` for the object and `F(U, V)` for the 
# restriction maps for admissible open subsets `U` and `V` of `X`.
#
# The objects and restriction maps must be cached. This can be 
# automated by using a concrete instance of `PreSheafOnScheme` in the 
# internals and making it available as `underlying_presheaf` of `F`. 
#
# If an object `M` for `U` is not found in the cache, then it is 
# produced by a call to 
#
#   produce_object(F, U)
#
# and then cached. The method of this function must thus be overwritten 
# for the concrete type of `F`. 
#
# Similarly for the production of the restriction maps: In this case 
# the user needs to implement a method of 
#
#   produce_restriction_map(F, U, V)
#
# for the concrete type of sheaf. 
########################################################################
@doc raw"""
    space(F::AbsPreSheaf)

For a sheaf ``ℱ`` on a space ``X`` return ``X``.
"""
function space(F::AbsPreSheaf)
  return space(underlying_presheaf(F))
end

function scheme(F::AbsPreSheaf)
  return scheme(underlying_presheaf(F))
end

scheme(F::PreSheafOnScheme) = space(F)

@doc raw"""
    (F::AbsPreSheaf)(U; cached=true)

For a sheaf ``ℱ`` on a space ``X`` and an (admissible) open set
``U ⊂ X`` check whether ``U`` is open in ``X`` and return ``ℱ(U)``.
"""
function (F::AbsPreSheaf{<:Any, OpenType, OutputType})(U::T; cached::Bool=true) where {OpenType, OutputType, T<:OpenType}
  cached && haskey(object_cache(F), U) && return object_cache(F)[U]
  result = produce_object(F, U)::OutputType
  cached && (object_cache(F)[U] = result)
  return result
end

function produce_object(F::AbsPreSheaf, U)
  error("method for `produce_object` must be overwritten for sheaves of type $(typeof(F))")
end

@doc raw"""
    restriction_map(F::AbsPreSheaf, U, V)

For a sheaf ``ℱ`` on a space ``X`` and an (admissible) pair of
open sets ``U, V ⊂ X`` check whether ``U ⊂ V ⊂ X`` are open and
return the restriction map ``ℱ(V) → ℱ(U)``.
"""
function restriction_map(F::AbsPreSheaf{<:Any, OpenType, OutputType, RestrictionType},
    U::Type1, V::Type2;
    check::Bool=false
  ) where {OpenType, OutputType, RestrictionType, Type1<:OpenType, Type2<:OpenType}
  inc = incoming_restrictions(F, F(V))
  inc !== nothing && haskey(inc, U) && return (inc[U])::RestrictionType

  # Check whether the given pair is even admissible.
  check && (is_open_func(F)(V, U) || error("the second argument is not open in the first"))

  # Hand the production of the restriction over to the internal method
  rho = produce_restriction_map(F, U, V)

  # Sanity checks
  # disabled for the moment because of the ideal sheaves: They use the ring maps.
  #domain(rho) === F(U) || error("domain of the produced restrition is not correct")
  #codomain(rho) === F(V) || error("codomain of the produced restrition is not correct")

  # Cache the result in the attributes of F(V)
  inc isa IdDict{<:OpenType, <:RestrictionType} && (inc[U] = rho) # It is the restriction coming from U.
  return rho::RestrictionType

  return restriction_map(underlying_presheaf(F), U, V)::RestrictionType
end

function produce_restriction_map(F::AbsPreSheaf, U, V)
  error("method for `produce_restriction_map` must be overwritten for sheaves of type $(typeof(F))")
end

### Temporary workaround
produce_restriction_map(F::AbsPreSheaf, U::AbsAffineScheme, V::AbsAffineScheme) = restriction_func(F)(F, U, V)

# An alias for shorter notation
function (F::AbsPreSheaf{<:Any, OpenType, OutputType, RestrictionType})(
    U::Type1, V::Type2
  ) where {OpenType, OutputType, RestrictionType, Type1<:OpenType, Type2<:OpenType}
  return restriction_map(F, U, V)::RestrictionType
end

@doc raw"""
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

function object_cache(F::AbsPreSheaf)
  return object_cache(underlying_presheaf(F))
end

########################################################################
# Implementation of PreSheafOnScheme                                   #
########################################################################

### implementing the essential functionality
space(F::PreSheafOnScheme) = F.X
object_cache(F::PreSheafOnScheme) = F.obj_cache # an IdDict caching the values of F on open subsets
#restriction_cache(F::PreSheafOnScheme) = F.res_cache # Caching is now done via attributes in the objects
is_open_func(F::PreSheafOnScheme) = F.is_open_func

### Production and caching of the values of F on admissible open sets
function (F::PreSheafOnScheme{<:Any, OpenType, OutputType})(U::T; cached::Bool=true, check::Bool=false) where {OpenType, OutputType, T<:OpenType}
  error("execution should never get here; please implement `produce_object` for your type of sheaf!")
end

### Production and caching of the restriction maps
@doc raw"""
    restriction_map(F::PreSheafOnScheme{<:Any, OpenType, OutputType, RestrictionType},
        U::Type1, V::Type2
      ) where {OpenType, OutputType, RestrictionType, Type1<:OpenType, Type2<:OpenType}

For a `F` produce (and cache) the restriction map `F(U) → F(V)`.
"""
function restriction_map(F::PreSheafOnScheme{<:Any, OpenType, OutputType, RestrictionType},
    U::Type1, V::Type2;
    check::Bool=false
  ) where {OpenType, OutputType, RestrictionType, Type1<:OpenType, Type2<:OpenType}
  # First, look up whether this restriction had already been asked for previously.
  inc = incoming_restrictions(F, F(V))
  inc !== nothing && haskey(inc, U) && return (inc[U])::RestrictionType

  # Check whether the given pair is even admissible.
  check && (is_open_func(F)(V, U) || error("the second argument is not open in the first"))

  # Hand the production of the restriction over to the internal method
  rho = produce_restriction_map(F, U, V)

  # Sanity checks
  # disabled for the moment because of the ideal sheaves: They use the ring maps.
  #domain(rho) === F(U) || error("domain of the produced restrition is not correct")
  #codomain(rho) === F(V) || error("codomain of the produced restrition is not correct")

  # Cache the result in the attributes of F(V)
  inc isa IdDict{<:OpenType, <:RestrictionType} && (inc[U] = rho) # It is the restriction coming from U.
  return rho::RestrictionType
end

@doc raw"""
    add_incoming_restriction!(F::AbsPreSheaf{<:Any, OpenType, <:Any, RestrictionType},
      U::OpenType,
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
  incoming_res === nothing && return F # This indicates that no
  incoming_res::IdDict{<:OpenType, <:RestrictionType}
  # sanity checks
  domain(rho) === F(U) || error("domain is not correct")
  codomain(rho) === M || error("codomain is not correct")
  incoming_res[U] = rho
  return F
end

@doc raw"""
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
  #  * U::SimplifiedAffineScheme with one ancestor W in the basic charts of X
  #  * U::PrincipalOpenSubset ⊂ V::PrincipalOpenSubset with ambient_scheme(U) === ambient_scheme(V) in the basic charts of X
  #  * U::PrincipalOpenSubset ⊂ V::PrincipalOpenSubset with ambient_scheme(U) != ambient_scheme(V) both in the basic charts of X
  #    and U and V contained in the gluing domains of their ambient schemes
  #  * U::AbsAffineScheme ⊂  V::AbsAffineScheme in the basic charts of X
  #  * U::AbsAffineScheme ⊂ X for U in the basic charts
  #  * U::PrincipalOpenSubset ⊂ X with ambient_scheme(U) in the basic charts of X
  #  * W::AffineSchemeOpenSubscheme ⊂ X with ambient_scheme(U) in the basic charts of X
  function is_open_func(
      U::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme},
      V::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme}
    )
    inc_U_flat = _flatten_open_subscheme(U, default_covering(X))
    inc_V_flat = _flatten_open_subscheme(V, default_covering(X))
    A = ambient_scheme(codomain(inc_U_flat))
    B = ambient_scheme(codomain(inc_V_flat))
    Udirect = codomain(inc_U_flat)
    Vdirect = codomain(inc_V_flat)

    if A === B
      is_subset(Udirect, Vdirect) || return false
    else
      G = default_covering(X)[A, B] # Get the gluing
      f, g = gluing_morphisms(G)
      is_subset(Udirect, domain(f)) || return false
      gU = preimage(g, Udirect)
      is_subset(gU, Vdirect) || return false
    end
    return true
  end
  function is_open_func(
      U::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme},
      Y::AbsCoveredScheme
    )
    return Y === X && has_ancestor(W->(any(WW->(WW === W), affine_charts(X))), U)
  end
  function is_open_func(U::AbsAffineScheme, Y::AbsCoveredScheme)
    return Y === X && has_ancestor(W->any(WW->(WW===W), affine_charts(X)), U)
  end
  # The following is implemented for the sake of completeness for boundary cases.
  function is_open_func(Z::AbsCoveredScheme, Y::AbsCoveredScheme)
    return X === Y === Z
  end
  function is_open_func(U::AbsAffineScheme, V::AbsAffineScheme)
    any(x->x===U, affine_charts(X)) || return false
    any(x->x===U, affine_charts(X)) || return false
    G = default_covering(X)[U, V]
    return is_subscheme(U, gluing_domains(G)[1])
  end
  function is_open_func(
      U::AbsAffineScheme,
      V::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme}
    )
    is_subscheme(U, V) && return true
    any(x->x===U, affine_charts(X)) || return false
    inc_V_flat = _flatten_open_subscheme(V, default_covering(X))
    A = ambient_scheme(codomain(inc_V_flat))
    Vdirect = codomain(inc_V_flat)
    W = ambient_scheme(Vdirect)
    haskey(gluings(default_covering(X)), (W, U)) || return false # In this case, they are not glued
    G = default_covering(X)[W, U]
    f, g = gluing_morphisms(G)
    pre_V = preimage(g, V)
    return is_subset(U, pre_V)
  end
  function is_open_func(
      U::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme},
      V::AbsAffineScheme
    )
    any(x->x===V, affine_charts(X)) || return false
    inc_U_flat = _flatten_open_subscheme(U, default_covering(X))
    A = ambient_scheme(codomain(inc_U_flat))
    Udirect = codomain(inc_U_flat)
    W = ambient_scheme(Udirect)
    haskey(gluings(default_covering(X)), (W, V)) || return false # In this case, they are not glued
    G = default_covering(X)[W, V]
    return is_subset(Udirect, gluing_domains(G)[1])
  end
  function is_open_func(W::AffineSchemeOpenSubscheme, Y::AbsCoveredScheme)
    return Y === X && ambient_scheme(W) in default_covering(X)
  end

  function is_open_func(W::AffineSchemeOpenSubscheme, V::AbsAffineScheme)
    V in default_covering(X) || return false
    ambient_scheme(W) === V && return true
    U = ambient_scheme(W)
    U in default_covering(X) || return false
    G = default_covering(X)[U, V]
    return is_subset(W, gluing_domains(G)[1])
  end
  function is_open_func(
      W::AffineSchemeOpenSubscheme,
      V::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme}
    )
    PW = ambient_scheme(W)
    inc_V_flat = _flatten_open_subscheme(V, default_covering(X))
    Vdirect = codomain(inc_V_flat)
    PV = ambient_scheme(Vdirect)
    PW in default_covering(X) || return false
    PV in default_covering(X) || return false
    if PW === PV
      return is_subscheme(W, V)
    else
      haskey(gluings(default_covering(X)), (PW, PV)) || return false
      G = default_covering(X)[PW, PV]
      preV = preimage(gluing_morphisms(G)[1], Vdirect)
      return is_subscheme(W, preV)
    end
  end
  function is_open_func(W::AffineSchemeOpenSubscheme, V::AffineSchemeOpenSubscheme)
    PW = ambient_scheme(W)
    PV = ambient_scheme(V)
    PW in default_covering(X) || return false
    PV in default_covering(X) || return false
    if PW === PV
      return is_subscheme(W, V)
    else
      G = default_covering(X)[PW, PV]
      preV = preimage(gluing_morphisms(G)[1], V)
      return is_subscheme(W, preV)
    end
  end
  function is_open_func(U::AbsAffineScheme, W::AffineSchemeOpenSubscheme)
    U in default_covering(X) || return false
    if U === ambient_scheme(W)
      # in this case W must be equal to U
      return is_subscheme(W, U)
      #return one(OO(U)) in complement_ideal(W)
    else
      G = default_covering(X)[ambient_scheme(W), U]
      is_subscheme(U, gluing_domains(G)[2]) || return false
      preU = preimage(gluing_morphisms(G)[1], U)
      return is_subscheme(preU, W)
    end
  end
  function is_open_func(
      U::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme},
      W::AffineSchemeOpenSubscheme
    )
    inc_U_flat = _flatten_open_subscheme(U, default_covering(X))
    A = ambient_scheme(codomain(inc_U_flat))
    Udirect = codomain(inc_U_flat)
    U_flat = codomain(inc_U_flat)
    PU = ambient_scheme(U_flat)
    if PU === ambient_scheme(W)
      # in this case W must be equal to U
      return is_subscheme(W, U_flat)
    else
      G = default_covering(X)[ambient_scheme(W), PU]
      is_subscheme(U_flat, gluing_domains(G)[2]) || return false
      preU = preimage(gluing_morphisms(G)[1], U_flat)
      return is_subscheme(preU, W)
    end
  end

  return is_open_func
end

function _is_open_func_for_schemes_without_affine_scheme_open_subscheme(X::AbsCoveredScheme)
  ### Checks for open containment.
  #
  # We allow the following cases:
  #
  #  * U::PrincipalOpenSubset with one ancestor W in the basic charts of X
  #  * U::SimplifiedAffineScheme with one ancestor W in the basic charts of X
  #  * U::PrincipalOpenSubset ⊂ V::PrincipalOpenSubset with ambient_scheme(U) === ambient_scheme(V) in the basic charts of X
  #  * U::PrincipalOpenSubset ⊂ V::PrincipalOpenSubset with ambient_scheme(U) != ambient_scheme(V) both in the basic charts of X
  #    and U and V contained in the gluing domains of their ambient schemes
  #  * U::AbsAffineScheme ⊂  V::AbsAffineScheme in the basic charts of X
  #  * U::AbsAffineScheme ⊂ X for U in the basic charts
  #  * U::PrincipalOpenSubset ⊂ X with ambient_scheme(U) in the basic charts of X
  function is_open_func(
      U::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme},
      V::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme}
    )
    inc_U_flat = _flatten_open_subscheme(U, default_covering(X))
    inc_V_flat = _flatten_open_subscheme(V, default_covering(X))
    A = ambient_scheme(codomain(inc_U_flat))
    B = ambient_scheme(codomain(inc_V_flat))
    Udirect = codomain(inc_U_flat)
    Vdirect = codomain(inc_V_flat)

    if A === B
      is_subset(Udirect, Vdirect) || return false
    else
      G = default_covering(X)[A, B] # Get the gluing
      f, g = gluing_morphisms(G)
      is_subscheme(Udirect, domain(f)) || return false
      gU = preimage(g, Udirect)
      is_subscheme(gU, Vdirect) || return false
    end
    return true
  end
  function is_open_func(
      U::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme},
      Y::AbsCoveredScheme
    )
    return Y === X && has_ancestor(W->(any(WW->(WW === W), affine_charts(X))), U)
  end
  function is_open_func(U::AbsAffineScheme, Y::AbsCoveredScheme)
    return Y === X && has_ancestor(W->any(WW->(WW===W), affine_charts(X)), U)
  end
  # The following is implemented for the sake of completeness for boundary cases.
  function is_open_func(Z::AbsCoveredScheme, Y::AbsCoveredScheme)
    return X === Y === Z
  end
  function is_open_func(U::AbsAffineScheme, V::AbsAffineScheme)
    U in affine_charts(X) || return false
    V in affine_charts(X) || return false
    G = default_covering(X)[U, V]
    return is_subscheme(U, gluing_domains(G)[1])
  end
  function is_open_func(
      U::Union{<:PrincipalOpenSubset, <:SimplifiedAffineScheme},
      V::AbsAffineScheme
    )
    V in affine_charts(X) || return false
    inc_U_flat = _flatten_open_subscheme(U, default_covering(X))
    A = ambient_scheme(codomain(inc_U_flat))
    Udirect = codomain(inc_U_flat)
    W = ambient_scheme(Udirect)
    haskey(gluings(default_covering(X)), (W, V)) || return false # In this case, they are not glued
    G = default_covering(X)[W, V]
    return is_subset(Udirect, gluing_domains(G)[1])
  end
  return is_open_func
end

underlying_presheaf(S::StructureSheafOfRings) = S.OO

@doc raw"""
    OO(X::AbsCoveredScheme) -> StructureSheafOfRings

Given a covered scheme `X`, return the structure sheaf of rings `𝒪(X)` of `X`.

# Examples
```jldoctest
julia> P, (x, y, z) = graded_polynomial_ring(QQ, [:x, :y, :z]);

julia> I = ideal([x^3-y^2*z]);

julia> Y = proj(P, I);

julia> Ycov = covered_scheme(Y)
Scheme
  over rational field
with default covering
  described by patches
    1: scheme(-(y//x)^2*(z//x) + 1)
    2: scheme((x//y)^3 - (z//y))
    3: scheme((x//z)^3 - (y//z)^2)
  in the coordinate(s)
    1: [(y//x), (z//x)]
    2: [(x//y), (z//y)]
    3: [(x//z), (y//z)]

julia> OO(Ycov)
Structure sheaf of rings of regular functions
  on scheme over QQ covered with 3 patches
    1: [(y//x), (z//x)]   scheme(-(y//x)^2*(z//x) + 1)
    2: [(x//y), (z//y)]   scheme((x//y)^3 - (z//y))
    3: [(x//z), (y//z)]   scheme((x//z)^3 - (y//z)^2)
```
"""
@attr StructureSheafOfRings function OO(X::AbsCoveredScheme)
  return StructureSheafOfRings(X)
end

function Base.show(io::IO, R::StructureSheafOfRings)
  io = pretty(io)
  if is_terse(io)
    print(io, "Structure sheaf of rings")
  else
    if is_unicode_allowed()
      print(io, "𝒪_{")
      print(terse(io), space(R), "}")
    else
      print(io, "Structure sheaf of ", Lowercase(), space(R))
    end
  end
end

function Base.show(io::IO, ::MIME"text/plain", R::StructureSheafOfRings)
  io = pretty(io)
  println(io, "Structure sheaf of rings of regular functions")
  print(io, Indent(), "on ", Lowercase())
  show(IOContext(io, :show_semi_compact => true), space(R))
  print(io, Dedent())
end

### Missing methods for compatibility of SimpleGluings with Gluings
function restrict(
    a::MPolyRingElem,
    U::PrincipalOpenSubset;
    check::Bool=true
  )
  parent(a) === OO(ambient_scheme(U)) || return OO(U)(a)
  return OO(U)(a, check=false)
end

function restrict(
    a::MPolyLocRingElem,
    U::PrincipalOpenSubset;
    check::Bool=true
  )
  parent(a) === OO(ambient_scheme(U)) || return OO(U)(a)
  return OO(U)(a, check=check)
end

function restrict(
    a::MPolyQuoRingElem,
    U::PrincipalOpenSubset;
    check::Bool=true
  )
  parent(a) === OO(ambient_scheme(U)) || return OO(U)(lift(a))
  return OO(U)(a, check=check)
end

function restrict(
    a::MPolyQuoLocRingElem,
    U::PrincipalOpenSubset;
    check::Bool=true
  )
  #parent(a) === OO(ambient_scheme(U)) || return OO(U)(lift(a))
  parent(a) === OO(ambient_scheme(U)) || return convert(OO(U), fraction(a))
  return OO(U)(a, check=check)
end

function restrict(
    a::Union{MPolyRingElem, MPolyQuoRingElem,
             MPolyLocRingElem, MPolyQuoLocRingElem},
    U::AffineSchemeOpenSubscheme;
    check::Bool=true
  )
  parent(a) === OO(ambient_scheme(U)) || return OO(U)(a)
  return OO(U)(a, check=check)
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
@doc raw"""
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

