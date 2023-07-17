


########################################################################
# (1) Properties for AbsSpecMor
########################################################################

# In an attempt to mimic inheritance, any new concrete instance of AbsSpecMor
# can internally store a SpecMor instance and must implement the method
# underlying_morphism to access it. Then, all functionality for SpecMor is
# automatically forwarded.
underlying_morphism(f::AbsSpecMor) = error("`underlying_morphism(f)` not implemented for `f` of type $(typeof(f))")


@doc raw"""
    domain(f::AbsSpecMor)

On a morphism ``f : X → Y`` of affine schemes, this returns ``X``.

# Examples
```jldoctest
julia> Y = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates x1, x2, x3

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
    of multivariate polynomial ring in 3 variables over QQ
    by ideal(x1)

julia> f = inclusion_morphism(X, Y)
Morphism
  from [x1, x2, x3]  spec of quotient of multivariate polynomial ring
  to   [x1, x2, x3]  affine 3-space over QQ
given by the pullback function
  x1 -> 0
  x2 -> x2
  x3 -> x3

julia> domain(f)
Spectrum
  of quotient
    of multivariate polynomial ring in 3 variables over QQ
    by ideal(x1)
```
"""
domain(f::AbsSpecMor) = domain(underlying_morphism(f))


@doc raw"""
    codomain(f::AbsSpecMor)

On a morphism ``f : X → Y`` of affine schemes, this returns ``Y``.

# Examples
```jldoctest
julia> Y = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates x1, x2, x3

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
    of multivariate polynomial ring in 3 variables over QQ
    by ideal(x1)

julia> f = inclusion_morphism(X, Y)
Morphism
  from [x1, x2, x3]  spec of quotient of multivariate polynomial ring
  to   [x1, x2, x3]  affine 3-space over QQ
given by the pullback function
  x1 -> 0
  x2 -> x2
  x3 -> x3

julia> codomain(f)
Affine space of dimension 3
  over rational field
with coordinates x1, x2, x3
```
"""
codomain(f::AbsSpecMor) = codomain(underlying_morphism(f))


@doc raw"""
    pullback(f::AbsSpecMor)

On a morphism ``f : X → Y`` of affine schemes ``X = Spec(S)`` and
``Y = Spec(R)``, this returns the ring homomorphism ``f^* : R → S``.

# Examples
```jldoctest
julia> Y = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates x1, x2, x3

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
    of multivariate polynomial ring in 3 variables over QQ
    by ideal(x1)

julia> pullback(inclusion_morphism(X, Y))
Map with following data
Domain:
=======
Multivariate polynomial ring in 3 variables over QQ
Codomain:
=========
Quotient of multivariate polynomial ring by ideal with 1 generator
```
"""
pullback(f::AbsSpecMor) = pullback(underlying_morphism(f))



########################################################################
# (2) Properties for SpecMor
########################################################################

pullback(phi::SpecMor) = phi.pullback
domain(phi::SpecMor) = phi.domain
codomain(phi::SpecMor) = phi.codomain



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
    phi::AbsSpecMor,
    Z::AbsSpec{<:Ring, <:MPolyQuoLocRing{<:Any, <:Any, <:Any, <:Any,
                                               <:MPolyPowersOfElement}};
    check::Bool=true
  )
  X = domain(phi)
  Y = codomain(phi)
  check && (issubset(Z, Y) || (Z = intersect(Y, Z)))
  IZ = modulus(underlying_quotient(OO(Z)))
  a = denominators(inverted_set(OO(Z)))
  R = ambient_coordinate_ring(X)
  f = pullback(phi)
  new_units = elem_type(R)[lifted_numerator(f(d)) for d in a]
  new_gens = Vector{elem_type(R)}(lifted_numerator.(f.(gens(IZ))))
  return hypersurface_complement(subscheme(X, new_gens), new_units)
end

function preimage(f::AbsSpecMor, Z::AbsSpec{<:Ring, <:MPolyRing}; check::Bool=true)
  OO(Z) == ambient_coordinate_ring(codomain(f)) || error("schemes can not be compared")
  return subscheme(domain(f), ideal(OO(domain(f)), [zero(OO(domain(f)))]))
end

function preimage(f::AbsSpecMor,
    Z::AbsSpec{<:Ring, <:MPolyLocRing{<:Any, <:Any, <:Any, <:Any,
                                            <:MPolyPowersOfElement}};
    check::Bool=true)
  return hypersurface_complement(domain(f), pullback(f).(denominators(inverted_set(OO(Z)))))
end

function preimage(f::AbsSpecMor, Z::AbsSpec; check::Bool=true)
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
    graph(f::AbsSpecMor)

Return the graph of ``f : X → Y`` as a subscheme of ``X×Y`` as well as the two projections to ``X`` and ``Y``.

# Examples
```jldoctest
julia> Y = affine_space(QQ,3)
Affine space of dimension 3
  over rational field
with coordinates x1, x2, x3

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
    of multivariate polynomial ring in 3 variables over QQ
    by ideal(x1)

julia> f = inclusion_morphism(X, Y)
Morphism
  from [x1, x2, x3]  spec of quotient of multivariate polynomial ring
  to   [x1, x2, x3]  affine 3-space over QQ
given by the pullback function
  x1 -> 0
  x2 -> x2
  x3 -> x3

julia> graph(f)
(Spec of quotient of multivariate polynomial ring, Morphism: spec of quotient of multivariate polynomial ring -> spec of quotient of multivariate polynomial ring, Morphism: spec of quotient of multivariate polynomial ring -> affine 3-space over QQ with coordinates x1, x2, x3)
```
"""
function graph(f::AbsSpecMor{<:AbsSpec{BRT}, <:AbsSpec{BRT}}) where {BRT}
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

@attr AbsSpecMor function inverse(f::AbsSpecMor)
  is_isomorphism(f) || error("the given morphism is not an isomorphism")
  return get_attribute(f, :inverse)::morphism_type(codomain(f), domain(f))
end

########################################################################
# (7) Type getters
########################################################################

### Type getters
pullback_type(::Type{T}) where {DomType, CodType, PbType, T<:AbsSpecMor{DomType, CodType, PbType}} = PbType
pullback_type(f::AbsSpecMor) = pullback_type(typeof(f))
domain_type(::Type{T}) where {DomType, CodType, PbType, T<:AbsSpecMor{DomType, CodType, PbType}} = DomType
domain_type(f::AbsSpecMor) = domain_type(typeof(f))
codomain_type(::Type{T}) where {DomType, CodType, PbType, T<:AbsSpecMor{DomType, CodType, PbType}} = CodType
codomain_type(f::AbsSpecMor) = codomain_type(typeof(f))

function morphism_type(::Type{SpecType1}, ::Type{SpecType2}) where {SpecType1<:AbsSpec, SpecType2<:AbsSpec}
  return SpecMor{SpecType1, SpecType2, morphism_type(ring_type(SpecType2), ring_type(SpecType1))}
end

morphism_type(X::AbsSpec, Y::AbsSpec) = morphism_type(typeof(X), typeof(Y))

@doc raw"""
    isomorphism_on_open_subsets(f::AbsSpecMor)

For a birational morphism ``f : X → Y`` of `AbsSpec`s this 
returns an isomorphism of affine schemes ``f' : U → V`` which is 
the restriction of ``f`` to two dense open subsets ``U ⊂ X`` and 
``V ⊂ Y``.
"""
function isomorphism_on_open_subsets(f::AbsSpecMor)
  if !has_attribute(f, :iso_on_open_subset)
    is_birational(f) # Should compute and store the attribute
  end
  return get_attribute(f, :iso_on_open_subset)::AbsSpecMor
end

