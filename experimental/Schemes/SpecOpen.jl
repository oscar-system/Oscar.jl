export SpecOpen, ambient, gens, complement, npatches, affine_patches, intersections, name, intersect, issubset, closure, find_non_zero_divisor, is_non_zero_divisor, is_dense, open_subset_type, ambient_type, is_canonically_isomorphic

export SpecOpenRing, scheme, domain, OO, structure_sheaf_ring_type, isdomain_type, isexact_type

export SpecOpenRingElem, domain, restrictions, patches, restrict, npatches, structure_sheaf_elem_type

export SpecOpenMor, maps_on_patches, restriction, identity_map, preimage, generic_fractions, pullback, maximal_extension

export adjoint

@Markdown.doc """
    SpecOpen{SpecType, BRT, BRET} <: Scheme{BRT, BRET}

Zariski open subset ``U`` of an affine scheme ``X = Spec(R)``. 
This stores a list of generators ``f₁,…,fᵣ`` of an ideal 
``I`` defining the complement ``Z = X ∖ U``. 
The scheme ``X`` is referred to as the *ambient scheme* and 
the list ``f₁,…,fᵣ`` as the *generators* for ``U``.

The type parameters stand for the following: The ring 
``R = (𝕜[x₁,…,xₙ]/I)[S⁻¹]`` is a localization of the quotient 
of a polynomial ring and 

 * `SpecType` is the type of the affine scheme ``X`` of which 
this is an open subset;
 * `BRT` is the type of the coefficient ring ``𝕜``;
 * `BRET` is the type of the elements of ``𝕜``.
"""
@attributes mutable struct SpecOpen{SpecType, BRT, BRET} <: Scheme{BRT, BRET}
  X::SpecType # the ambient scheme
  gens::Vector # a list of functions defining the complement of the open subset

  # fields used for caching
  name::String
  patches::Vector{SpecType}
  intersections::Dict{Tuple{Int, Int}, SpecType}
  complement::SpecType

  function SpecOpen(
      X::SpecType, 
      f::Vector{RET}; 
      name::String="", 
      check::Bool=true
    ) where {SpecType<:Spec, RET<:RingElem}
    for a in f
      parent(a) == base_ring(OO(X)) || error("element does not belong to the correct ring")
      if check
        !isempty(X) && iszero(OO(X)(a)) && error("generators must not be zero")
      end
    end
    U = new{SpecType, typeof(base_ring(X)), elem_type(base_ring(X))}(X, f)
    length(name) > 0 && set_name!(U, name)
    return U
  end
end

open_subset_type(X::Spec) = SpecOpen{typeof(X), typeof(coefficient_ring(base_ring(OO(X)))), elem_type(coefficient_ring(base_ring(OO(X))))}
open_subset_type(::Type{Spec{BRT, BRET, RT, RET, MST}}) where {BRT, BRET, RT, RET, MST} = SpecOpen{Spec{BRT, BRET, RT, RET, MST}, BRT, BRET}

ambient_type(U::SpecOpen{SpecType, BRT, BRET}) where {SpecType<:Spec, BRT, BRET} = SpecType
ambient_type(::Type{SpecOpen{SpecType, BRT, BRET}}) where {SpecType<:Spec, BRT, BRET} = SpecType

poly_type(::Type{SpecOpenType}) where {SpecOpenType<:SpecOpen} = poly_type(ambient_type(SpecOpenType))
poly_type(U::SpecOpen) = poly_type(typeof(U))

@Markdown.doc """
    ambient(U::SpecOpen)

Return the ambient scheme ``X`` of a Zariski open subset ``U ⊂ X``.
"""
ambient(U::SpecOpen) = U.X

@Markdown.doc """
    npatches(U::SpecOpen)

Return the number of generators stored for describing the complement of ``U``.
"""
npatches(U::SpecOpen) = length(U.gens)

@Markdown.doc """
    gens(U::SpecOpen)

Return the generators ``[f₁,…,fᵣ]`` stored for the description 
of the complement of ``U``.
"""
gens(U::SpecOpen) = U.gens::Vector{elem_type(base_ring(OO(ambient(U))))}

@Markdown.doc """
    affine_patch(U::SpecOpen, i::Int)

Return the hypersurface complement of ``fᵢ`` in the 
ambient scheme ``X`` of ``U`` where ``f₁,…,fᵣ`` are 
the generators stored for the description of the complement 
of ``U``. This function can also be called using the 
`getindex` method or simply via `U[i]`.
"""
affine_patch(U::SpecOpen, i::Int) = affine_patches(U)[i]
getindex(U::SpecOpen, i::Int) = affine_patches(U)[i]

function getindex(U::SpecOpen, X::Spec) 
  for i in 1:npatches(U)
    X == U[i] && return i
  end
  error("scheme $X not found among the open patches in $U")
end

function Base.show(io::IO, U::SpecOpen)
  if isdefined(U, :name) 
    print(io, name(U))
    return
  end
  print(io, "complement of zero locus of $(gens(U)) in $(ambient(U))")
end

@Markdown.doc """
    function SpecOpen(X::Spec, I::MPolyLocalizedIdeal)

Return the complement of the zero locus of ``I`` in ``X``.
"""
function SpecOpen(X::Spec, I::MPolyLocalizedIdeal)
  base_ring(I) === localized_ring(OO(X)) || error("Ideal does not belong to the correct ring")
  f = [reduce(f, groebner_basis(localized_modulus(OO(X)))) for f in gens(I)]
  g = [numerator(a) for a in f if !iszero(numerator(a))]
  return SpecOpen(X, g)
end

function SpecOpen(X::Spec, I::MPolyIdeal)
  return SpecOpen(X, localized_ring(OO(X))(I))
end

function complement(X::T, Z::T) where {T<:Spec}
  if !issubset(Z, X) 
    Z = intersect(X, Z)
  end
  return SpecOpen(X, modulus(OO(Z)))
end

SpecOpen(X::Spec) = SpecOpen(X, [one(base_ring(OO(X)))], check=false)

function complement(U::SpecOpen) 
  if !isdefined(U, :complement)
    #I = radical(saturated_ideal(ideal(localized_ring(OO(ambient(U))), gens(U))))
    #U.complement = subscheme(ambient(U), I)
    U.complement = subscheme(ambient(U), gens(U))
  end
  return U.complement
end

@Markdown.doc """
    affine_patches(U::SpecOpen)

Return a list of principal affine open subschemes covering ``U``.
"""
function affine_patches(U::SpecOpen)
  if !isdefined(U, :patches)
    X = ambient(U)
    U.patches = [hypersurface_complement(X, f) for f in gens(U)]
  end
  return U.patches
end

@Markdown.doc """
    intersections(U::SpecOpen)

Return a list of pairwise intersections of the 
principal open subschemes covering ``U``.
"""
function intersections(U::SpecOpen)
  if !isdefined(U, :intersections)
    X = ambient(U)
    V = affine_patches(U)
    for i in 2:length(V)
      for j in 1:i-1
        U.intersections[(i,j)] = U.intersections[(j,i)] = intersect(V[i], V[j])
      end
    end
  end
  return U.intersections
end

function name(U::SpecOpen) 
  if isdefined(U, :name)
    return U.name
  end
  return "open subset of $(ambient(U))"
end

function intersect(
    Y::Spec, 
    U::SpecOpen;
    check::Bool=true
  )
  X = ambient(U)
  base_ring(OO(X)) === base_ring(OO(Y)) || error("Schemes can not be compared")
  if !issubset(Y, X)
    Y = intersect(Y, X)
  end
  return SpecOpen(Y, gens(U), check=check)
end

function intersect(
    U::SpecOpen,
    Y::Spec
  )
  return intersect(Y, U)
end

function intersect(
    U::SpecOpen,
    V::SpecOpen
  )
  X = ambient(U) 
  is_canonically_isomorphic(X, ambient(V)) || error("ambient schemes do not coincide")
  return SpecOpen(X, [a*b for a in gens(U) for b in gens(V)])
end

function Base.union(U::T, V::T) where {T<:SpecOpen}
  is_canonically_isomorphic(ambient(U), ambient(V)) || error("the two open sets do not lay in the same ambient scheme")
  return SpecOpen(ambient(U), vcat(gens(U), gens(V)))
end

function issubset(
    Y::Spec, 
    U::SpecOpen
  )
  return one(localized_ring(OO(Y))) in ideal(OO(Y), gens(U))
end

function issubset(
    U::SpecOpen,
    Y::Spec
  ) 
  for V in affine_patches(U)
    issubset(V, Y) || return false
  end
  return true
end

function issubset(U::T, V::T) where {T<:SpecOpen}
  base_ring(OO(ambient(U))) === base_ring(OO(ambient(V))) || return false
  W = V
  issubset(W, ambient(U)) || (W = intersect(V, ambient(U)))
  Z = complement(W)
  # perform an implicit radical membership test (Rabinowitsch) that is way more 
  # efficient than computing radicals.
  for g in gens(U)
    isempty(hypersurface_complement(Z, g)) || return false
  end
  return true
  #return issubset(complement(intersect(V, ambient(U))), complement(U))
end

function is_canonically_isomorphic(U::T, V::T) where {T<:SpecOpen}
  return issubset(U, V) && issubset(V, U)
end

function is_canonically_isomorphic(
    U::SpecOpen,
    Y::Spec
  )
  return issubset(U, Y) && issubset(Y, U)
end

function is_canonically_isomorphic(
    Y::Spec,
    U::SpecOpen
  )
  return is_canonically_isomorphic(U, Y)
end

@Markdown.doc """
    closure(U::SpecOpen)

Compute the Zariski closure of an open set ``U ⊂ X`` 
where ``X`` is the affine ambient scheme of ``U``.
"""
function closure(U::SpecOpen)
  X = ambient(U)
  R = base_ring(OO(X))
  I = saturated_ideal(localized_modulus(OO(X)))
  I = saturation(I, ideal(R, gens(U)))
  return subscheme(X, I)
end

@Markdown.doc """
    closure(U::SpecOpen, Y::Spec)

Compute the closure of ``U ⊂ Y``.
"""
function closure(
    U::SpecOpen,
    Y::Spec 
  )
  issubset(U, Y) || error("the first set is not contained in the second")
  X = closure(U)
  return intersect(X, Y)
end


@Markdown.doc """
    SpecOpenRing{SpecType, OpenType}

The ring of regular functions ``𝒪(X, U)`` on an open subset ``U`` of an 
affine scheme ``X``.

 * `SpecType` is the type of the affine scheme ``X`` on which 
this sheaf is defined;
 * `OpenType` is the type of the (Zariski) open subsets of ``U``.
"""
mutable struct SpecOpenRing{SpecType, OpenType} <: Ring
  scheme::SpecType
  domain::OpenType

  function SpecOpenRing(
      X::SpecType, 
      U::OpenType
    ) where {SpecType<:Spec, OpenType<:SpecOpen}
    issubset(U, X) || error("open set does not lay in the scheme")
    return new{SpecType, OpenType}(X, U)
  end
end

SpecOpenRing(U::SpecOpen) = SpecOpenRing(ambient(U), U)

spec_open_ring_type(::Type{T}) where {T<:Spec} = SpecOpenRing{T, open_subset_type(T)}
spec_open_ring_type(X::Spec) = spec_open_ring_type(typeof(X))

@Markdown.doc """
    scheme(R::SpecOpenRing)

The ring ``R = 𝒪(X, U)`` belongs to a sheaf of rings ``𝒪(X, -)`` and this returns 
the scheme ``X`` on which ``𝒪`` is defined.
"""
scheme(R::SpecOpenRing) = R.scheme

@Markdown.doc """
    domain(R::SpecOpenRing)

For a ring ``R = 𝒪(X, U)``, return ``U``.
"""
domain(R::SpecOpenRing) = R.domain

OO(U::SpecOpen) = SpecOpenRing(ambient(U), U)
OO(X::Spec, U::SpecOpen) = SpecOpenRing(X, U)

function is_canonically_isomorphic(R::T, S::T) where {T<:SpecOpenRing}
  return is_canonically_isomorphic(scheme(R), scheme(S)) && is_canonically_isomorphic(domain(R), domain(S))
end

@Markdown.doc """
    SpecOpenRingElem{SpecOpenType, RestrictionType}

An element ``f ∈ 𝒪(X, U)`` of the ring of regular functions on 
an open set ``U`` of an affine scheme ``X``.

 * `SpecOpenType` is the type of the open sets ``U`` of ``X``;
 * `RestrictionType` is the type of the restrictions of ``f`` to
the affine patches of ``U``.
"""
mutable struct SpecOpenRingElem{
      SpecOpenRingType<:SpecOpenRing, 
      RestrictionType<:MPolyQuoLocalizedRingElem
    } <: RingElem
  parent::SpecOpenRingType
  restrictions::Vector{RestrictionType}

  function SpecOpenRingElem(
      R::SpecOpenRingType,
      f::Vector{RestrictionType};
      check::Bool=true
    ) where {
        SpecOpenRingType<:SpecOpenRing, 
        RestrictionType<:MPolyQuoLocalizedRingElem
    }
    n = length(f)
    U = domain(R)
    n == length(affine_patches(U)) || error("the number of restrictions does not coincide with the number of affine patches")
    g = RestrictionType[OO(U[i])(f[i]) for i in 1:n] # will throw if conversion is not possible
    if check
      for i in 1:n-1
        for j in i+1:n
          W = intersections(U)[(i,j)]
          OO(W)(f[i], check=false) == OO(W)(f[j], check=false) || error("elements are not compatible on overlap")
        end
      end
    end
    return new{SpecOpenRingType, RestrictionType}(R, g)
  end
end

### type getters
elem_type(::Type{SpecOpenRing{S, T}}) where {S, T} = SpecOpenRingElem{SpecOpenRing{S, T}, elem_type(ring_type(S))}

elem_type(R::SpecOpenRing) = elem_type(typeof(R))

parent_type(::Type{SpecOpenRingElem{S, T}}) where {S, T} = S
parent_type(f::SpecOpenRingElem) = parent_type(typeof(f))

parent(f::SpecOpenRingElem) = f.parent
scheme(f::SpecOpenRingElem) = scheme(parent(f))
domain(f::SpecOpenRingElem) = domain(parent(f))
restrictions(f::SpecOpenRingElem) = f.restrictions
affine_patches(f::SpecOpenRingElem) = affine_patches(domain(f))
npatches(f::SpecOpenRingElem) = length(restrictions(f))
getindex(f::SpecOpenRingElem, i::Int) = getindex(restrictions(f), i)
getindex(f::SpecOpenRingElem, U::Spec) = restrictions(f)[domain(f)[U]]

### copying
function Base.deepcopy_internal(f::SpecOpenRingElem, dict::IdDict)
  return SpecOpenRingElem(parent(f), copy(restrictions(f)), check=false)
end

function restrict(
    f::SpecOpenRingElem, 
    V::Spec
  )
  isempty(V) && return zero(OO(V))
  for i in 1:length(restrictions(f))
    if V == affine_patches(domain(f))[i]
      return restrictions(f)[i]
    end
  end
  issubset(V, domain(f)) || error("the set is not contained in the domain of definition of the function")
  VU = [intersect(V, U) for U in affine_patches(domain(f))]
  g = [OO(VU[i])(f[i]) for i in 1:length(VU)]
  l = write_as_linear_combination(one(OO(V)), OO(V).(lifted_denominator.(g)))
  return dot(l, OO(V).(lifted_numerator.(g)))
end

(R::MPolyQuoLocalizedRing)(f::SpecOpenRingElem) = restrict(f, Spec(R))

(R::SpecOpenRing)(f::RingElem) = SpecOpenRingElem(R, [OO(U)(f) for U in affine_patches(domain(R))])
(R::SpecOpenRing)(f::Vector{T}) where {T<:RingElem} = SpecOpenRingElem(R, [OO(domain(R)[i])(f[i]) for i in 1:length(f)])

function (R::SpecOpenRing)(f::SpecOpenRingElem) 
  parent(f) === R && return f
  return SpecOpenRingElem(R, [restrict(f, U) for U in affine_patches(domain(R))])
end

function +(a::T, b::T) where {T<:SpecOpenRingElem}
  parent(a) === parent(b) || return a + (parent(a)(b))
  return SpecOpenRingElem(parent(a), [a[i] + b[i] for i in 1:length(restrictions(a))], check=false)
end

function -(a::T, b::T) where {T<:SpecOpenRingElem}
  parent(a) === parent(b) || return a - (parent(a)(b))
  return SpecOpenRingElem(parent(a), [a[i] - b[i] for i in 1:length(restrictions(a))], check=false)
end

function -(a::T) where {T<:SpecOpenRingElem}
  return SpecOpenRingElem(parent(a), [-a[i] for i in 1:length(restrictions(a))], check=false)
end

function *(a::T, b::T) where {T<:SpecOpenRingElem}
  parent(a) === parent(b) || return a * (parent(a)(b))
  return SpecOpenRingElem(parent(a), [a[i] * b[i] for i in 1:length(restrictions(a))], check=false)
end

function *(a::RingElem, b::T) where {T<:SpecOpenRingElem}
  return b*(parent(b)(a))
end

function *(a::Integer, b::T) where {T<:SpecOpenRingElem}
  return b*(parent(b)(a))
end

function *(b::T, a::RingElem) where {T<:SpecOpenRingElem}
  return a*b
end

function ==(a::T, b::T) where {T<:SpecOpenRingElem}
  parent(a) === parent(b) || return a == (parent(a)(b))
  for i in 1:length(restrictions(a))
    a[i] == b[i] || return false
  end
  return true
end

function ^(a::SpecOpenRingElem, i::Int64)
  return SpecOpenRingElem(parent(a), [a[k]^i for k in 1:length(restrictions(a))])
end
function ^(a::SpecOpenRingElem, i::Integer)
  return SpecOpenRingElem(parent(a), [a[k]^i for k in 1:length(restrictions(a))])
end
function ^(a::SpecOpenRingElem, i::fmpz)
  return SpecOpenRingElem(parent(a), [a[k]^i for k in 1:length(restrictions(a))])
end

function divexact(a::T, b::T; check::Bool=false) where {T<:SpecOpenRingElem} 
  parent(a) === parent(b) || return divexact(a, (parent(a)(b)))
  return SpecOpenRingElem(parent(a), [divexact(a[i], b[i]) for i in 1:length(restrictions(a))])
end

function isunit(a::SpecOpenRingElem) 
  return all(x->isunit(x), restrictions(a))
end

inv(a::SpecOpenRingElem) = SpecOpenRingElem(parent(a), [inv(f) for f in restrictions(a)], check=false)

one(R::SpecOpenRing) = SpecOpenRingElem(R, [one(OO(U)) for U in affine_patches(domain(R))], check=false)
zero(R::SpecOpenRing) = SpecOpenRingElem(R, [zero(OO(U)) for U in affine_patches(domain(R))], check=false)
(R::SpecOpenRing)() = zero(R)
(R::SpecOpenRing)(a::Integer) = SpecOpenRingElem(R, [OO(U)(a) for U in affine_patches(domain(R))], check=false)
(R::SpecOpenRing)(a::Int64) = SpecOpenRingElem(R, [OO(U)(a) for U in affine_patches(domain(R))], check=false)
(R::SpecOpenRing)(a::fmpz) = SpecOpenRingElem(R, [OO(U)(a) for U in affine_patches(domain(R))], check=false)

isdomain_type(::Type{T}) where {T<:SpecOpenRingElem} = true
isdomain_type(a::SpecOpenRingElem) = isdomain_type(typeof(a))
isexact_type(::Type{T}) where {T<:SpecOpenRingElem} = true
isexact_type(a::SpecOpenRingElem) = isexact_type(typeof(a))
isdomain_type(::Type{T}) where {T<:SpecOpenRing} = true
isdomain_type(R::SpecOpenRing) = isdomain_type(typeof(R))
isexact_type(::Type{T}) where {T<:SpecOpenRing} = true
isexact_type(R::SpecOpenRing) = isexact_type(typeof(R))

AbstractAlgebra.promote_rule(::Type{T}, ::Type{RET}) where {T<:SpecOpenRingElem, RET<:Integer} = T
AbstractAlgebra.promote_rule(::Type{RET}, ::Type{T}) where {T<:SpecOpenRingElem, RET<:Integer} = T

### TODO: Rethink this. For instance, restrictions can happen both from and to Specs.
function AbstractAlgebra.promote_rule(::Type{T}, ::Type{RET}) where {T<:SpecOpenRingElem, RET<:RingElem} 
  return T
end


@Markdown.doc """
    maximal_extension(X::Spec, f::AbstractAlgebra.Generic.Frac)

Return the maximal extension of the restriction of ``f`` 
to a rational function on ``X`` on a maximal domain of 
definition ``U ⊂ X``. 

**Note:** When ``X = Spec(R)`` with ``R = (𝕜[x₁,…,xₙ]/I)[S⁻¹]``, 
the numerator and denominator of ``f`` have to be elements of 
the ring ``𝕜[x₁,…,xₙ]``.
"""
function maximal_extension(
    X::Spec, 
    f::AbstractAlgebra.Generic.Frac{RET}
  ) where {RET<:RingElem}
  a = numerator(f)
  b = denominator(f)
  W = localized_ring(OO(X))
  I = quotient(ideal(W, b) + localized_modulus(OO(X)), ideal(W, a))
  U = SpecOpen(X, I)
  g = [OO(V)(f) for V in affine_patches(U)]
  R = SpecOpenRing(X, U)
  return R(g)
end

@Markdown.doc """
    maximal_extension(X::Spec, f::Vector{AbstractAlgebra.Generic.Frac})

Return the extension of the restriction of the ``fᵢ`` as a 
set of rational functions on ``X`` as *regular* functions to a 
common maximal domain of definition ``U ⊂ X``.

**Note:** When ``X = Spec(R)`` with ``R = (𝕜[x₁,…,xₙ]/I)[S⁻¹]``, 
the numerators and denominators of the entries of ``f`` have to 
be elements of the ring ``𝕜[x₁,…,xₙ]``.
"""
function maximal_extension(
    X::Spec, 
    f::Vector{AbstractAlgebra.Generic.Frac{RET}}
  ) where {RET<:RingElem}
  if length(f) == 0
    return SpecOpen(X), Vector{structure_sheaf_elem_type(X)}()
  end
  R = base_ring(parent(f[1]))
  for a in f
    R == base_ring(parent(a)) || error("fractions do not belong to the same ring")
  end
  R == base_ring(OO(X)) || error("fractions do not belong to the base ring of the scheme")
  a = numerator.(f)
  b = denominator.(f)
  W = localized_ring(OO(X))
  I = ideal(W, one(W))
  for p in f
    I = intersect(quotient(ideal(W, denominator(p)) + localized_modulus(OO(X)), ideal(W, numerator(p))), I)
  end
  U = SpecOpen(X, I)
  S = SpecOpenRing(X, U)
  # TODO: For some reason, the type of the inner vector is not inferred if it has no entries. 
  # Investigate why? Type instability?
  return U, [SpecOpenRingElem(S, (elem_type(OO(X))[OO(V)(a) for V in affine_patches(U)])) for a in f]
end

#TODO: implement the catchall versions of the above functions.

@Markdown.doc """
    SpecOpenMor{DomainType<:SpecOpen, CodomainType<:SpecOpen, SpecMorType<:SpecMor}

Morphisms ``f : U → V`` of open sets ``U ⊂ X`` and ``V ⊂ Y`` of affine schemes.
These are stored as morphisms ``fᵢ: Uᵢ→ Y`` on the affine patches 
``Uᵢ`` of ``U``.

The type parameters stand for the following: When ``X = Spec(R)`` and 
``Y = Spec(S)`` with ``R = (𝕜[x₁,…,xₘ]/I)[A⁻¹]`` and ``S = (𝕜[y₁,…,yₙ]/J)[B⁻¹]``
then 

 * `DomainType` is the type of the domain;
 * `CodomainType` is the type of the codomain;
 * `SpecMorType` is the type of the restriction of the morphism to the
affine patches of the domain to the affine ambient scheme of the codomain. 
"""
mutable struct SpecOpenMor{DomainType<:SpecOpen, CodomainType<:SpecOpen, SpecMorType<:SpecMor}
  domain::DomainType
  codomain::CodomainType
  maps_on_patches::Vector{SpecMorType}

  # fields used for caching
  inverse::SpecOpenMor

  function SpecOpenMor(
      U::DomainType,
      V::CodomainType,
      f::Vector{SpecMorType};
      check::Bool=true
    ) where {DomainType<:SpecOpen, CodomainType<:SpecOpen, SpecMorType<:SpecMor}
    Y = ambient(V)
    n = length(f)
    n == length(affine_patches(U)) || error("number of patches does not coincide with the number of maps")
    if check
      for i in 1:n
        domain(f[i]) == affine_patches(U)[i] || error("domain of definition of the map does not coincide with the patch")
        codomain(f[i]) == Y || error("codomain is not compatible")
      end
      for i in 1:n-1
	for j in i+1:n
	  A = intersect(domain(f[i]), domain(f[j]))
	  restrict(f[i], A, Y) == restrict(f[j], A, Y) || error("maps don't glue")
	end
      end
      I = ideal(localized_ring(OO(Y)), gens(V))
      for g in f
	one(localized_ring(OO(domain(g)))) in pullback(g)(I) + localized_modulus(OO(domain(g))) || error("image is not contained in the codomain")
      end
    end
    return new{DomainType, CodomainType, SpecMorType}(U, V, f)
  end
end

morphism_type(U::S, V::T) where {S<:SpecOpen, T<:SpecOpen} = SpecOpenMor{S, T, morphism_type(ambient(U), ambient(V))}
morphism_type(::Type{DomainType}, ::Type{CodomainType}) where {DomainType<:SpecOpen, CodomainType<:SpecOpen} = SpecOpenMor{DomainType, CodomainType, morphism_type(ambient_type(DomainType), ambient_type(CodomainType))}

domain(f::SpecOpenMor) = f.domain
codomain(f::SpecOpenMor) = f.codomain
maps_on_patches(f::SpecOpenMor) = f.maps_on_patches
getindex(f::SpecOpenMor, i::Int) = maps_on_patches(f)[i]

function getindex(f::SpecOpenMor, D::Spec)
  U = affine_patches(domain(f))
  D in U || error("affine patch not found in the domain")
  for i = 1:length(U)
    D == U[i] && return f[i]
  end
end

function Base.show(io::IO, f::SpecOpenMor) 
  print(io, "Morphism from $(domain(f)) to $(codomain(f)) given by the rational map $(generic_fractions(f))")
end

function SpecOpenMor(U::T, V::T, f::Vector; check::Bool=true) where {T<:SpecOpen}
  Y = ambient(V)
  return SpecOpenMor(U, V, [SpecMor(W, Y, f) for W in affine_patches(U)], check=check) 
end

function inclusion_morphism(U::T, V::T; check::Bool=true) where {T<:SpecOpen}
  X = ambient(U)
  if check 
    issubset(U, V) || error("method not implemented")
  end
  return SpecOpenMor(U, V, gens(base_ring(OO(X))), check=false)
end

function SpecOpenMor(X::SpecType, d::RET, Y::SpecType, e::RET, f::Vector{RET}; check::Bool=true) where {SpecType<:Spec, RET<:RingElem}
  U = SpecOpen(X, [d], check=check)
  V = SpecOpen(Y, [e], check=check)
  return SpecOpenMor(U, V, [SpecMor(U[1], Y, OO(U[1]).(f), check=check)], check=check)
end




@Markdown.doc """
    compose(f::T, g::T) where {T<:SpecOpenMor}

Compute the composition of two morphisms 

      f    g
    U →  V →  W
    ∩    ∩    ∩
    X    Y    Z

of open sets of affine varieties.
"""
function compose(f::T, g::T; check::Bool=true) where {T<:SpecOpenMor}
  U = domain(f)
  Cf = codomain(f)
  V = domain(g)
  if check
    issubset(Cf, V) || error("maps are not compatible")
  end
  W = codomain(g)
  X = ambient(U)
  Y = ambient(V)
  Z = ambient(W)
  g_maps = maps_on_patches(g)
  f_maps = maps_on_patches(f)
  #####################################################################
  # The method proceeds as follows.
  # We extract the affine open sets 
  #   Uᵢⱼ = Uᵢ∩ f⁻¹(Vⱼ)  
  # where Uᵢ are the patches of U and Vⱼ those of V.
  # Then we can pass to the restrictions 
  #   fᵢⱼ : Uᵢⱼ → Vⱼ 
  # and the compositions 
  #   gᵢⱼ :=  gⱼ ∘ fᵢⱼ: Uᵢⱼ → Z.
  # For any variable z for Z we can write the pullbacks gᵢⱼ*z ∈ 𝒪(Uᵢⱼ). 
  # These functions necessarily coincide on the overlaps of the Uᵢⱼ in 
  # Uᵢ, so we can extend them to a global function on Uᵢ. 
  # From these global functions, we can then assemble a morphism 
  # on each one of the affine schemes Uᵢ which glue together to a 
  # SpecOpenMor over these patches.
  #####################################################################
  m = length(gens(U))
  n = length(gens(V))
  result_maps = Vector{morphism_type(X, Y)}()
  for i in 1:m
    U_i = affine_patches(U)[i]
    f_i = f_maps[i]
    z_i = Vector{elem_type(OO(U_i))}() # the pullbacks of the coordinate functions of Z on Uᵢ
    z_ij = Vector{Vector{elem_type(OO(U_i))}}()
    for j in 1:n
      V_j = affine_patches(V)[j]
      g_j = g_maps[j]
      U_ij = intersect(U_i, preimage(f_i, V_j))
      g_ij = compose(restrict(f_i, U_ij, V_j), g_j)
      push!(z_ij, [pullback(g_ij)(z) for z in gens(OO(Z))])
      ### This is using the well known trick that if aᵢ// dᵢ coincide pairwise
      # wherever dᵢ and dⱼ are defined, and 1 = ∑ᵢ λᵢ⋅ dᵢ, then 
      # ∑ᵢλᵢ⋅ aᵢ = aᵢ// dᵢ on all open sets {dᵢ≠ 0}. 
    end
    for k in 1:length(gens(OO(Z)))
      l = write_as_linear_combination(one(OO(U_i)), lifted_denominator.([z_ij[j][k] for j in 1:n]))
      push!(z_i, dot(l, OO(U_i).(lifted_numerator.([z_ij[j][k] for j in 1:n]))))
    end
    push!(result_maps, SpecMor(U_i, Z, MPolyQuoLocalizedRingHom(OO(Z), OO(U_i), z_i)))
    #push!(result_maps, SpecMor(U_i, Z, z_i))
  end
  return SpecOpenMor(U, W, result_maps)
end

function pullback(f::SpecOpenMor, a::RET) where {RET<:RingElem}
  U = domain(f)
  X = ambient(U)
  V = codomain(f)
  Y = ambient(V)
  R = base_ring(OO(Y))
  parent(a) == R || error("element does not belong to the correct ring")
  pb_a = elem_type(OO(X))[pullback(f[i])(a) for i in 1:npatches(U)]
  return SpecOpenRingElem(SpecOpenRing(X, U), pb_a)
end

function pullback(f::SpecOpenMor, a::SpecOpenRingElem)
  error("not implemented")
end

@Markdown.doc """
    maximal_extension(X::Spec, Y::Spec, f::AbstractAlgebra.Generic.Frac)

Given a rational map ``ϕ : X ---> Y ⊂ Spec 𝕜[y₁,…,yₙ]`` of affine schemes 
determined by ``ϕ*(yⱼ) = fⱼ = aⱼ/bⱼ``, find the maximal open subset ``U⊂ X`` 
to which ``ϕ`` can be extended to a regular map ``g : U → Y`` and return ``g``.
"""
function maximal_extension(
    X::SpecType,
    Y::SpecType,
    f::Vector{AbstractAlgebra.Generic.Frac{RET}}
  ) where {SpecType<:Spec, RET<:RingElem}
  U, g = maximal_extension(X, f)
  n = length(affine_patches(U))
  maps = Vector{morphism_type(X, Y)}()
  for i in 1:n
    push!(maps, SpecMor(affine_patches(U)[i], Y, [restrictions(a)[i] for a in g]))
  end
  return SpecOpenMor(U, SpecOpen(Y), maps)
end

function maximal_extension(
    X::T, 
    Y::T, 
    f::Vector
  ) where {T<:Spec}
  return maximal_extension(X, Y, FractionField(base_ring(OO(X))).(f))
end

### the restriction of a morphism to open subsets in domain and codomain
function restriction(
    f::SpecOpenMor,
    U::SpecOpen,
    V::SpecOpen;
    check::Bool=true
  )
  if check
    issubset(U, domain(f)) || error("the given open is not an open subset of the domain of the map")
    issubset(V, codomain(f)) || error("the given open is not an open subset of the codomain of the map")
  end
  inc = inclusion_morphism(U, domain(f), check=check)
  help_map = compose(inc, f, check=check)
  return SpecOpenMor(U, V, maps_on_patches(help_map), check=check)
end

### the restriction of a morphism to closed subsets in domain and codomain
function restriction(
    f::SpecOpenMor,
    X::SpecType,
    Y::SpecType;
    check::Bool=true
  ) where {SpecType<:Spec}
  U = intersect(X, domain(f), check=check)
  V = intersect(Y, codomain(f), check=check)

  new_maps_on_patches = [restrict(f[i], U[i], Y, check=check) for i in 1:npatches(U)]

  return SpecOpenMor(U, V, new_maps_on_patches, check=check)
end

identity_map(U::SpecOpen) = SpecOpenMor(U, U, [SpecMor(V, ambient(U), gens(OO(V)), check=false) for V in affine_patches(U)], check=false)

function ==(f::T, g::T) where {T<:SpecOpenMor} 
  is_canonically_isomorphic(domain(f), domain(g)) || return false
  is_canonically_isomorphic(codomain(f), codomain(g)) || return false
  Y = ambient(codomain(f))
  m = length(affine_patches(domain(f)))
  n = length(affine_patches(domain(g)))
  for i in 1:m
    for j in 1:n
      U = intersect(domain(f)[i], domain(g)[i])
      restrict(f[i], U, Y) == restrict(g[i], U, Y) || return false
    end
  end
  return true
end

function preimage(f::SpecOpenMor, Z::Spec)
  U = domain(f) 
  X = ambient(U)
  n = length(affine_patches(U))
  W = localized_ring(OO(X))
  I = ideal(W, one(W))
  for i in 1:n 
    I = intersect(I, localized_modulus(OO(closure(preimage(f[i], Z), X))))
  end
  fZbar = subscheme(X, I)
  return SpecOpen(fZbar, [g for g in gens(U) if !iszero(OO(fZbar)(g))])
end

function preimage(f::SpecOpenMor, V::SpecOpen)
  U = domain(f)
  X = ambient(U)
  R = base_ring(OO(X))
  I = ideal(R, one(R))
  for i in 1:npatches(U)
    I = intersect(I, saturated_ideal(ideal(OO(U[i]), pullback(f[i]).(gens(V)))))
  end
  return intersect(U, SpecOpen(X, I))
end

function is_non_zero_divisor(f::RET, U::SpecOpen) where {RET<:RingElem}
  for V in affine_patches(U)
    zero_ideal = ideal(OO(V), [zero(OO(V))])
    zero_ideal == quotient(zero_ideal, ideal(OO(V), [f])) || return false
  end
  return true
end

function find_non_zero_divisor(U::SpecOpen)
  n = length(gens(U))
  X = ambient(U)
  R = base_ring(OO(X))
  n == 0 && return zero(R)
  kk = coefficient_ring(R)
  coeff = elem_type(kk)[rand(kk, 0:100) for i in 1:n]
  d = sum([coeff[i]*gens(U)[i] for i in 1:n])
  while !is_non_zero_divisor(d, U)
    d = dot([rand(kk, 0:100) for i in 1:n], gens(U))
  end
  return d
end

@Markdown.doc """
    generic_fractions(f::SpecOpenMor)

Given a morphism ``f : U → V`` of Zariski open subsets ``U ⊂ X ⊂ 𝔸ᵐ`` and ``V ⊂ Y ⊂ 𝔸ⁿ``, 
produce a tuple of fractions ``[a₁/b₁,…,aₙ/bₙ]`` such that ``f`` can be recovered 
as the maximal extension of the rational map given by 
```
   U ⊃ U' → 𝔸ⁿ,  x ↦ [a₁(x)/b₁(x),…,aₙ(x)/bₙ(x)]
```
where ``U'`` is the complement of the zero loci of the denominators ``bᵢ`` in ``U``.
In particular, this requires ``U'`` to be dense in ``U`` and this subset 
is chosen at random.
"""
function generic_fractions(f::SpecOpenMor)
  U = domain(f)
  X = ambient(U)
  V = codomain(f)
  Y = ambient(V)
  d = find_non_zero_divisor(U)
  V = hypersurface_complement(X, d)
  result = fraction.([restrict(pullback(f, y), V) for y in gens(base_ring(OO(Y)))])
end

function is_dense(U::SpecOpen)
  X = ambient(U)
  I = [localized_modulus(OO(closure(V, X))) for V in affine_patches(U)]
  J = ideal(OO(X), [one(OO(X))])
  for i in I
    J = intersect(J, i)
  end
  return J == localized_modulus(OO(X))
end

#function Base.adjoint(M::AbstractAlgebra.Generic.MatSpaceElem{T}) where {T} 
function Base.adjoint(M::MatElem)
  n = ncols(M)
  n == nrows(M) || error("matrix is not square")
  N = zero(M)
  rows = collect(1:n)
  cols = collect(1:n)
  row_sign = 1
  for i in 1:n
    column_sign = row_sign
    for j in 1:n
      N[j, i] = column_sign* det(M[deleteat!(copy(rows), i), deleteat!(copy(cols), j)])
      column_sign = column_sign*(-1)
    end
    row_sign = row_sign*(-1)
  end
  return N
end

Base.inv(M::MatElem) = inv(det(M))*adjoint(M)

function preimage(f::SpecMor, V::SpecOpen; check::Bool=true)
  Z = preimage(f, ambient(V))
  new_gens = pullback(f).(gens(V))
  return SpecOpen(Z, lifted_numerator.(new_gens), check=check)
end
