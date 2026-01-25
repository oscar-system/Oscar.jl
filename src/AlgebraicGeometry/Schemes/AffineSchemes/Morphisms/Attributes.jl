


########################################################################
# (1) Properties for AbsAffineSchemeMor
########################################################################

# In an attempt to mimic inheritance, any new concrete instance of AbsAffineSchemeMor
# can internally store a AffineSchemeMor instance and must implement the method
# underlying_morphism to access it. Then, all functionality for AffineSchemeMor is
# automatically forwarded.
underlying_morphism(f::AbsAffineSchemeMor) = error("`underlying_morphism(f)` not implemented for `f` of type $(typeof(f))")


@doc raw"""
    domain(f::AbsAffineSchemeMor)

On a morphism ``f : X → Y`` of affine schemes, this returns ``X``.

# Examples
```jldoctest
julia> Y = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates [x1, x2, x3]

julia> R = OO(Y)
Multivariate polynomial ring in 3 variables x1, x2, x3
  over rational field

julia> (x1,x2,x3) = gens(R)
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> X = subscheme(Y, x1)
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    by ideal (x1)

julia> f = inclusion_morphism(X, Y)
Affine scheme morphism
  from [x1, x2, x3]  scheme(x1)
  to   [x1, x2, x3]  affine 3-space over QQ
given by the pullback function
  x1 -> x1
  x2 -> x2
  x3 -> x3

julia> domain(f)
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    by ideal (x1)
```
"""
domain(f::AbsAffineSchemeMor) = domain(underlying_morphism(f))


@doc raw"""
    codomain(f::AbsAffineSchemeMor)

On a morphism ``f : X → Y`` of affine schemes, this returns ``Y``.

# Examples
```jldoctest
julia> Y = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates [x1, x2, x3]

julia> R = OO(Y)
Multivariate polynomial ring in 3 variables x1, x2, x3
  over rational field

julia> (x1,x2,x3) = gens(R)
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> X = subscheme(Y, x1)
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    by ideal (x1)

julia> f = inclusion_morphism(X, Y)
Affine scheme morphism
  from [x1, x2, x3]  scheme(x1)
  to   [x1, x2, x3]  affine 3-space over QQ
given by the pullback function
  x1 -> x1
  x2 -> x2
  x3 -> x3

julia> codomain(f)
Affine space of dimension 3
  over rational field
with coordinates [x1, x2, x3]
```
"""
codomain(f::AbsAffineSchemeMor) = codomain(underlying_morphism(f))


@doc raw"""
    pullback(f::AbsAffineSchemeMor)

On a morphism ``f : X → Y`` of affine schemes ``X = Spec(S)`` and
``Y = Spec(R)``, this returns the ring homomorphism ``f^* : R → S``.

# Examples
```jldoctest
julia> Y = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates [x1, x2, x3]

julia> R = OO(Y)
Multivariate polynomial ring in 3 variables x1, x2, x3
  over rational field

julia> (x1,x2,x3) = gens(R)
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> X = subscheme(Y, x1)
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    by ideal (x1)

julia> pullback(inclusion_morphism(X, Y))
Ring homomorphism
  from multivariate polynomial ring in 3 variables over QQ
  to quotient of multivariate polynomial ring by ideal (x1)
defined by
  x1 -> x1
  x2 -> x2
  x3 -> x3
```
"""
pullback(f::AbsAffineSchemeMor) = pullback(underlying_morphism(f))



########################################################################
# (2) Properties for AffineSchemeMor
########################################################################

pullback(phi::AffineSchemeMor) = phi.pullback
domain(phi::AffineSchemeMor) = phi.domain
codomain(phi::AffineSchemeMor) = phi.codomain



########################################################################
# (3) Properties for OpenInclusion
########################################################################

underlying_morphism(f::OpenInclusion) = f.inc
complement_ideal(f::OpenInclusion) = f.I
complement_scheme(f::OpenInclusion) = f.Z


########################################################################
# (4) Preimage
########################################################################

function preimage(
    phi::AbsAffineSchemeMor,
    Z::AbsAffineScheme{<:Ring, <:MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any,
                                               <:MPolyPowersOfElement}};
    check::Bool=true
  )
  X = domain(phi)
  Y = codomain(phi)
  check && (is_subscheme(Z, Y) || (Z = intersect(Y, Z)))
  IZ = modulus(underlying_quotient(OO(Z)))
  a = denominators(inverted_set(OO(Z)))
  R = ambient_coordinate_ring(X)
  f = pullback(phi)
  new_units = elem_type(R)[lifted_numerator(f(d)) for d in a]
  new_gens = Vector{elem_type(R)}(lifted_numerator.(f.(gens(IZ))))
  return hypersurface_complement(subscheme(X, new_gens), new_units)
end

function preimage(f::AbsAffineSchemeMor, Z::AbsAffineScheme{<:Ring, <:MPolyRing}; check::Bool=true)
  OO(Z) == ambient_coordinate_ring(codomain(f)) || error("schemes can not be compared")
  return subscheme(domain(f), ideal(OO(domain(f)), [zero(OO(domain(f)))]))
end

function preimage(f::AbsAffineSchemeMor,
    Z::AbsAffineScheme{<:Ring, <:MPolyLocRing{<:Any, <:Any, <:Any, <:Any,
                                            <:MPolyPowersOfElement}};
    check::Bool=true)
  return hypersurface_complement(domain(f), pullback(f).(denominators(inverted_set(OO(Z)))))
end

function preimage(f::AbsAffineSchemeMor, Z::AbsAffineScheme; check::Bool=true)
  pbf = pullback(f)
  R = OO(codomain(f))
  S = OO(domain(f))
  I = ideal(R, gens(modulus(OO(Z))))
  J = ideal(S, pbf.(gens(I)))
  return subscheme(domain(f), J)
end



##############################################
# (5) The graph of a morphism
##############################################

@doc raw"""
    graph(f::AbsAffineSchemeMor)

Return the graph of ``f : X → Y`` as a subscheme of ``X×Y`` as well as the two projections to ``X`` and ``Y``.

# Examples
```jldoctest
julia> Y = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates [x1, x2, x3]

julia> R = OO(Y)
Multivariate polynomial ring in 3 variables x1, x2, x3
  over rational field

julia> (x1,x2,x3) = gens(R)
3-element Vector{QQMPolyRingElem}:
 x1
 x2
 x3

julia> X = subscheme(Y, x1)
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables x1, x2, x3
      over rational field
    by ideal (x1)

julia> f = inclusion_morphism(X, Y)
Affine scheme morphism
  from [x1, x2, x3]  scheme(x1)
  to   [x1, x2, x3]  affine 3-space over QQ
given by the pullback function
  x1 -> x1
  x2 -> x2
  x3 -> x3

julia> graph(f)
(scheme(x1, -x1, x2 - x2, x3 - x3), Hom: scheme(x1, -x1, x2 - x2, x3 - x3) -> scheme(x1), Hom: scheme(x1, -x1, x2 - x2, x3 - x3) -> affine 3-space over QQ with coordinates [x1, x2, x3])
```
"""
function graph(f::AbsAffineSchemeMor{<:AbsAffineScheme{BRT}, <:AbsAffineScheme{BRT}}) where {BRT}
  X = domain(f)
  Y = codomain(f)
  XxY, prX, prY = product(X, Y)
  pb_X = pullback(prX)
  pb_Y = pullback(prY)
  pb_f = pullback(f)
  I = ideal(OO(XxY), pb_X.(pb_f.(gens(OO(Y)))) - pb_Y.(gens(OO(Y))))
  G = subscheme(XxY, I)
  return G, restrict(prX, G, X), restrict(prY, G, Y)
end



##############################################
# (6) The inverse of a morphism
##############################################

@attr AbsAffineSchemeMor function inverse(f::AbsAffineSchemeMor)
  is_isomorphism(f) || error("the given morphism is not an isomorphism")
  return get_attribute(f, :inverse)::morphism_type(codomain(f), domain(f))
end

########################################################################
# (7) Type getters
########################################################################

### Type getters
pullback_type(::Type{T}) where {DomType, CodType, PbType, T<:AbsAffineSchemeMor{DomType, CodType, PbType}} = PbType
pullback_type(f::AbsAffineSchemeMor) = pullback_type(typeof(f))

function morphism_type(::Type{AffineSchemeType1}, ::Type{AffineSchemeType2}) where {AffineSchemeType1<:AbsAffineScheme, AffineSchemeType2<:AbsAffineScheme}
  return AffineSchemeMor{AffineSchemeType1, AffineSchemeType2, morphism_type(ring_type(AffineSchemeType2), ring_type(AffineSchemeType1))}
end

morphism_type(X::AbsAffineScheme, Y::AbsAffineScheme) = morphism_type(typeof(X), typeof(Y))


