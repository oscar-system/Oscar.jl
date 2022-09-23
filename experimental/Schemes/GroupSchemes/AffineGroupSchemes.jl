export AbsAffineGroupScheme
export multiplication_map, product_over_ground_field, inverse_map, first_inclusion, second_inclusion, diagonal_embedding, first_projection, second_projection

export AffineGroupScheme

########################################################################
#
# The following documents the interface for the functionality of an 
# abstract affine group scheme.
#
# Below you will find a concrete type that can be used for a default 
# implementation of these methods by forwarding it via the 
# `underlying_group_scheme` function.
#
########################################################################

abstract type AbsAffineGroupScheme{BRT, BRET} <: AbsSpec{BRT, BRET} end

@Markdown.doc """
    multiplication_map(G::AbsAffineGroupScheme)

Returns a morphism of affine schemes `G × G → G` 
where `G × G` is the `product_over_ground_field` of `G`.
"""
@attr SpecMor function multiplication_map(G::AbsAffineGroupScheme)
  return SpecMor(product_over_ground_field(G), G, 
                 pullback(multiplication_map(underlying_group_scheme(G)))
                )
end

@Markdown.doc """
    product_over_ground_field(G::AbsAffineGroupScheme)

Returns the fiber product `G × G` over the ground field 𝕜 of `G`, 
together with its two canonical projections. 
"""
function product_over_ground_field(G::AbsAffineGroupScheme)
  return product_over_ground_field(underlying_group_scheme(G))
end

@Markdown.doc """
    first_projection(G::AbsAffineGroupScheme)

For a group scheme ``G`` this returns the projection ``G × G → G, (g, h) ↦ g`` 
where ``G×G`` is the `product_over_ground_field` of ``G``.
"""
@attr SpecMor function first_projection(G::AbsAffineGroupScheme)
  return SpecMor(product_over_ground_field(G), G, 
                 pullback(first_projection(underlying_group_scheme(G)))
                )
end

@Markdown.doc """
    second_projection(G::AbsAffineGroupScheme)

For a group scheme ``G`` this returns the projection ``G × G → G, (g, h) ↦ h`` 
where ``G×G`` is the `product_over_ground_field` of ``G``.
"""
@attr SpecMor function second_projection(G::AbsAffineGroupScheme)
  return SpecMor(product_over_ground_field(G), G, 
                 pullback(second_projection(underlying_group_scheme(G)))
                )
end

@Markdown.doc """
   inverse_map(G::AbsAffineGroupScheme)

Returns the map `G → G` assigning to each point ``g ∈ G`` its inverse 
``g⁻¹`` with respect to the group law on ``G``.
"""
@attr SpecMor function inverse_map(G::AbsAffineGroupScheme)
  return SpecMor(G, G, pullback(inverse_map(underlying_group_scheme(G))))
end

@Markdown.doc """
    first_inclusion(G::AbsAffineGroupScheme)

For a group ``G`` with neutral element ``e ∈ G`` this returns the morphism 
``G → G × G, g ↦ (g, e)`` where ``G × G`` is the `product_over_ground_field` of ``G``.
"""
@attr SpecMor function first_inclusion(G::AbsAffineGroupScheme)
  return SpecMor(G, product_over_ground_field(G),
                 pullback(first_inclusion(underlying_group_scheme(G)))
                )
end

@Markdown.doc """
    second_inclusion(G::AbsAffineGroupScheme)

For a group ``G`` with neutral element ``e ∈ G`` this returns the morphism 
``G → G × G, g ↦ (e, g)`` where ``G × G`` is the `product_over_ground_field` of ``G``.
"""
@attr SpecMor function second_inclusion(G::AbsAffineGroupScheme)
  return SpecMor(G, product_over_ground_field(G),
                 pullback(second_inclusion(underlying_group_scheme(G)))
                )
end

@Markdown.doc """
    diagonal_embedding(G::AbsAffineGroupScheme)

For a group ``G`` this returns the morphism ``G → G × G, g ↦ (g, g)`` 
where ``G × G`` is the `product_over_ground_field` of ``G``.
"""
@attr SpecMor function diagonal_embedding(G::AbsAffineGroupScheme)
  return SpecMor(G, product_over_ground_field(G), 
                 pullback(diagonal_embedding(underlying_group_scheme(G)))
                )
end

@Markdown.doc """
    neutral_element_coordinates(G::AbsAffineGroupScheme)

For an affine group scheme ``G ⊂ 𝔸ⁿ`` this returns the 
coordinates ``(x₁,…,xₙ) ∈ 𝕜ⁿ`` of the neutral element ``e ∈ G``
for the given embedding.
"""
function neutral_element_coordinates(G::AbsAffineGroupScheme)
  return neutral_element_coordinates(underlying_group_scheme(G))
end

### This method can be implemented to forward an existing 
# implementation of the above interface
function underlying_group_scheme(G::AbsAffineGroupScheme)
  error("function `underlying_group_scheme` not implemented for schemes of type $(typeof(G))")
end

### Forwarding of the scheme interface
underlying_scheme(G::AbsAffineGroupScheme) = underlying_scheme(underlying_group_scheme(G))
underlying_scheme_type(::Type{T}) where {T<:AbsAffineGroupScheme} = underlying_scheme_type(underlying_group_scheme_type(T))


########################################################################
#
# A minimal implementation of the interface of an affine group scheme. 
#
# This will usually be used as an internal field to more specific 
# types of affine group schemes, see the example `_kk_star` below.
#
########################################################################
@attributes mutable struct AffineGroupScheme{BRT, BRET, SpecType} <: AbsAffineGroupScheme{BRT, BRET}
  X::SpecType
  product_over_ground_field::AbsSpec
  diagonal_embedding::SpecMor
  first_projection::SpecMor
  second_projection::SpecMor
  first_inclusion::SpecMor
  second_inclusion::SpecMor
  multiplication_map::SpecMor
  inverse_map::SpecMor
  neutral_element::Vector{BRET}

  function AffineGroupScheme(
      X::AbsSpec, XxX::AbsSpec, 
      diag::SpecMor,
      p1::SpecMor, p2::SpecMor, 
      i1::SpecMor, i2::SpecMor, 
      mult_map::SpecMor, inv::SpecMor, 
      neutral_element::Vector{T};
      check::Bool=true
    ) where {T<:FieldElem}
    all(x->(parent(x) == coefficient_ring(base_ring(OO(X)))), neutral_element) || error("coordinates of the neutral element do not belong to the correct field")
    domain(p1) == XxX || error("domain of first projection not compatible")
    domain(p2) == XxX || error("domain of second projection not compatible")
    codomain(p1) == X || error("codomain of first projection not compatible")
    codomain(p2) == X || error("codomain of second projection not compatible")
    domain(diag) == X || error("domain of diagonal embedding is not compatible")
    codomain(diag) == XxX || error("codomain of diagonal embedding is not compatible")
    domain(i1) == X || error("domain of first inclusion not compatible")
    domain(i2) == X || error("domain of second inclusion not compatible")
    codomain(i1) == XxX || error("codomain of first inclusion not compatible")
    codomain(i2) == XxX || error("codomain of second inclusion not compatible")
    domain(inv) == codomain(inv) == X || error("domain or codomain of the inverse map is not compatible")

    if check
      is_identity_map(compose(i1, p1)) || error("composition of the first inclusion and projection is not the identity")
      is_identity_map(compose(i2, p2)) || error("composition of the second inclusion and projection is not the identity")
      is_isomorphism(inv) || error("the inverse map is not an isomorphism")
      is_identity_map(compose(inv, inv)) || error("the composition of the inverse map with itself is not the identity")
      is_identity_map(compose(diag, p2)) || error("composition of the diagonal embedding and the second projection is not the identity")
      is_identity_map(compose(diag, p1)) || error("composition of the diagonal embedding and the first projection is not the identity")
      is_identity_map(compose(i1, mult_map)) || error("composition of the first inclusion and the multiplication map is not the identity")
      is_identity_map(compose(i2, mult_map)) || error("composition of the second inclusion and the multiplication map is not the identity")
      # TODO: Add some further checks about the neutral element.
    end

    G = new{
            base_ring_type(X), 
            ring_type(X), 
            typeof(X)
           }(
             X, XxX
             )
    # We need to manually promote the maps so that they have access 
    # to G itself as (co-)domains.
    G.diagonal_embedding = SpecMor(G, XxX, pullback(diag))
    G.first_projection = SpecMor(XxX, G, pullback(p1))
    G.second_projection = SpecMor(XxX, G, pullback(p2))
    G.first_inclusion = SpecMor(G, XxX, pullback(i1))
    G.second_inclusion = SpecMor(G, XxX, pullback(i2))
    G.multiplication_map = SpecMor(XxX, G, pullback(mult_map))
    G.inverse_map = SpecMor(G, G, pullback(inv))
    G.neutral_element = neutral_element
    return G
  end
end

### essential internal getters for forwarding of the scheme functionality
underlying_scheme(G::AffineGroupScheme) = G.X

### type getters
function underlying_scheme_type(
    ::Type{AffineGroupSchemeType}
  ) where {SpecType, AffineGroupSchemeType<:AffineGroupScheme{<:Any, <:Any, SpecType}}
  return SpecType
end

### user-facing essential getters
product_over_ground_field(G::AffineGroupScheme) = G.product_over_ground_field
diagonal_embedding(G::AffineGroupScheme) = G.diagonal_embedding
first_projection(G::AffineGroupScheme) = G.first_projection
second_projection(G::AffineGroupScheme) = G.second_projection
first_inclusion(G::AffineGroupScheme) = G.first_inclusion
second_inclusion(G::AffineGroupScheme) = G.second_inclusion
multiplication_map(G::AffineGroupScheme) = G.multiplication_map
inverse_map(G::AffineGroupScheme) = G.inverse_map
neutral_element_coordinates(G::AffineGroupScheme) = G.neutral_element


########################################################################
# 
# The following is an example of a concrete implementation of an 
# affine group scheme of type `<:AbsAffineGroupScheme` using the 
# concrete minimal type `AffineGroupScheme` in the background. 
#
########################################################################

@attributes mutable struct _kk_star{BRT, BRET, GroupSchemeType} <: AbsAffineGroupScheme{BRT, BRET}
  X::GroupSchemeType

  function _kk_star(kk::AbstractAlgebra.Field; var_name::String="x")
    P, (x,) = PolynomialRing(kk, [var_name])
    X = hypersurface_complement(Spec(P), x)
    XxX, p1, p2 = product(X, X)
    y = pullback(p1).(gens(OO(X)))[1]
    z = pullback(p2).(gens(OO(X)))[1]
    mult_map = SpecMor(XxX, X, hom(OO(X), OO(XxX), [y*z]))
    e = [one(kk)]
    inv_map = SpecMor(X, X, hom(OO(X), OO(X), [1//x]))
    i1 = SpecMor(X, XxX, hom(OO(XxX), OO(X), [x, one(x)]))
    i2 = SpecMor(X, XxX, hom(OO(XxX), OO(X), [one(x), x]))
    diag = SpecMor(X, XxX, hom(OO(XxX), OO(X), [x, x]))

    G = AffineGroupScheme(X, XxX, diag, p1, p2, i1, i2, mult_map, inv_map, e)

    return new{typeof(kk), elem_type(kk), typeof(G)}(G)
  end
end

### forwarding of the essential getters is achieved with this one line
underlying_group_scheme(G::_kk_star) = G.X

### necessary forwarding of the type getters
underlying_group_scheme_type(::Type{T}) where {GST, T<:_kk_star{<:Any, <:Any, GST}} = GST

### additional methods for this particular type can be added at will
function variable_name(G::_kk_star)
  return String(symbols(base_ring(OO(G)))[1])
end
    
