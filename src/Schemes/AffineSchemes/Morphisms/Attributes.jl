export underlying_morphism, domain, codomain, pullback, inverse, preimage
export pullback_type, domain_type, codomain_type, morphism_type



########################################################################
# (1) Properties for AbsSpecMor
########################################################################

# In an attempt to mimic inheritance, any new concrete instance of AbsSpecMor
# can internally store a SpecMor instance and must implement the method
# underlying_morphism to access it. Then, all functionality for SpecMor is
# automatically forwarded.
underlying_morphism(f::AbsSpecMor) = error("`underlying_morphism(f)` not implemented for `f` of type $(typeof(f))")


@Markdown.doc """
    domain(f::AbsSpecMor)

On a morphism ``f : X → Y`` of affine schemes, this returns ``X``.
"""
domain(f::AbsSpecMor) = domain(underlying_morphism(f))


@Markdown.doc """
    codomain(f::AbsSpecMor)

On a morphism ``f : X → Y`` of affine schemes, this returns ``Y``.
"""
codomain(f::AbsSpecMor) = codomain(underlying_morphism(f))


@Markdown.doc """
    pullback(f::AbsSpecMor)

On a morphism ``f : X → Y`` of affine schemes ``X = Spec(S)`` and
``Y = Spec(R)``, this returns the ring homomorphism ``f^* : R → S``.
"""
pullback(f::AbsSpecMor) = pullback(underlying_morphism(f))



########################################################################
# (2) Properties for SpecMor
########################################################################

pullback(phi::SpecMor) = phi.pullback
domain(phi::SpecMor) = phi.domain
codomain(phi::SpecMor) = phi.codomain



########################################################################
# (3) Properties for OpenInclusion                                            #
########################################################################

underlying_morphism(f::OpenInclusion) = f.inc
complement_ideal(f::OpenInclusion) = f.I
complement_scheme(f::OpenInclusion) = f.Z


########################################################################
# (4) Preimage
########################################################################

function preimage(
    phi::AbsSpecMor,
    Z::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any,
                                               <:MPolyPowersOfElement}};
    check::Bool=true
  )
  X = domain(phi)
  Y = codomain(phi)
  check && (issubset(Z, Y) || (Z = intersect(Y, Z)))
  IZ = modulus(quotient_ring(OO(Z)))
  a = denominators(inverted_set(OO(Z)))
  R = ambient_ring(X)
  f = pullback(phi)
  new_units = [lifted_numerator(f(d)) for d in a]
  new_gens = lifted_numerator.(f.(gens(IZ)))
  return hypersurface_complement(subscheme(X, new_gens), new_units)
end

function preimage(f::AbsSpecMor, Z::AbsSpec{<:Ring, <:MPolyRing}; check::Bool=true)
  OO(Z) == ambient_ring(codomain(f)) || error("schemes can not be compared")
  return subscheme(domain(f), ideal(OO(domain(f)), [zero(OO(domain(f)))]))
end

function preimage(f::AbsSpecMor,
    Z::AbsSpec{<:Ring, <:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any,
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

function graph(f::AbsSpecMor)
  X = standard_spec(domain(f))
  Y = standard_spec(codomain(f))
  fres = restrict(f, X, Y)
  G, prX, prY = graph(fres)
  return G, compose(prX, SpecMor(X, domain(f), gens(OO(X)))), compose(prY, SpecMor(Y, codomain(f), gens(OO(Y))))
end

function graph(f::AbsSpecMor{SpecType, SpecType}) where {SpecType<:StdSpec}
  X = domain(f)
  Y = codomain(f)
  XxY, prX, prY = product(X, Y)
  pb_X = pullback(prX)
  pb_Y = pullback(prY)
  pb_f = pullback(f)
  I = ideal(localized_ring(OO(XxY)), lift.(pb_X.(images(pb_f)) - pb_Y.(gens(OO(Y)))))
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
