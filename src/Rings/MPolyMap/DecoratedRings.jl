#const _DecoratedRingTypes = Union{MPolyDecRing, MPolyQuoRing{<:MPolyDecRingElem}}
#const _DecoratedRingElemTypes = Union{MPolyDecRingElem, MPolyQuoRingElem{<:MPolyDecRingElem}}

# Since the interface regarding the decorations for decorated rings and their quotients is not 
# yet coherent, we only allow free polynomial rings for the moment. 
#
# For instance `forget_grading` does not work for quotient rings.
const _DecoratedRingTypes = Union{MPolyDecRing}
const _DecoratedRingElemTypes = Union{MPolyDecRingElem}

@attributes mutable struct DecoratedRingMorphism{
    DomainType<:_DecoratedRingTypes,
    CodomainType<:_DecoratedRingTypes,
    RingMapType<:MPolyAnyMap,
    GradingGroupMorphismType,
    ImgGensType
  } <: Map{
    DomainType, CodomainType, Hecke.Hecke.Map, DecoratedRingMorphism
  }

  domain::DomainType
  codomain::CodomainType
  ring_map::RingMapType
  decoration_map::GradingGroupMorphismType
  img_gens::Vector{ImgGensType}


  function DecoratedRingMorphism(
      R::_DecoratedRingTypes, # the domain
      S::_DecoratedRingTypes, # the codomain
      g::Vector{<:_DecoratedRingElemTypes}, # the images of the generators
      phi::Map, # the associated homomorphism of the grading groups
      coeff_map::Any # the map on the coefficient rings
    )
    @assert all(x->ishomogeneous(x), g) "images of generators must be homogeneous"
    @assert all(i->degree(g[i]) == phi(degree(gen(R, i))), 1:ngens(R)) "morphism on grading groups is not compatible with the images of the generators"

    Phi = hom(forget_grading(R), forget_grading(S), coeff_map, forget_grading.(g))
    return new{typeof(R), typeof(S), typeof(Phi), typeof(phi), elem_type(S)}(R, S, Phi, phi, g)
  end

  function DecoratedRingMorphism(
      R::_DecoratedRingTypes, # the domain
      S::_DecoratedRingTypes, # the codomain
      g::Vector{<:_DecoratedRingElemTypes}, # the images of the generators
      phi::Map # the associated homomorphism of the grading groups
    )
    @assert length(g) == ngens(R) "number of images does not coincide with the number of variables"
    @assert all(x->parent(x)===S, g) "images are not elements of the codomain"
    @assert all(x->ishomogeneous(x), g) "images of generators must be homogeneous"
    @assert all(i->degree(g[i]) == phi(degree(gen(R, i))), 1:ngens(R)) "morphism on grading groups is not compatible with the images of the generators"

    Phi = hom(forget_grading(R), forget_grading(S), forget_grading.(g))
    return new{typeof(R), typeof(S), typeof(Phi), typeof(phi), elem_type(S)}(R, S, Phi, phi, g)
  end
  function DecoratedRingMorphism(
      R::_DecoratedRingTypes, # the domain
      S::_DecoratedRingTypes, # the codomain
      Phi::Map, # the homomorphism of the rings without decoration
      phi::Map # the associated homomorphism of the grading groups
    )
    @assert domain(Phi) === R "domains are incompatible"
    @assert codomain(Phi) === S "codomains are incompatible"
    g = S.(Phi.(gens(domain(Phi))))
    @assert all(i->degree(g[i]) == phi(degree(gen(R, i))), 1:ngens(R)) "morphism on grading groups is not compatible with the images of the generators"
    return new{typeof(R), typeof(S), typeof(Phi), typeof(phi), elem_type(S)}(R, S, Phi, phi, g)
  end
end
      
domain(psi::DecoratedRingMorphism) = psi.domain
codomain(psi::DecoratedRingMorphism) = psi.codomain
forget_grading(psi::DecoratedRingMorphism) = psi.ring_map
decoration_map(psi::DecoratedRingMorphism) = psi.decoration_map

# Note: The codomain of the map below does not necessarily coincide with the codomain of 
# psi, but elements can be coerced in the latter. 
coefficient_map(psi::DecoratedRingMorphism) = coefficient_map(forget_grading(psi))

########################################################################
# User facing constructors
########################################################################

function hom(R::_DecoratedRingTypes, S::_DecoratedRingTypes, g::Vector;
    decoration_map=begin
      G = grading_group(R)
      H = grading_group(S)
      # try to infer the morphism of grading groups from g
      gg_imgs = elem_type(H)[]
      @assert all(x->parent(x)===S, g) "images of generators do not have the correct parent"
      for v in gens(G)
        i = findfirst(x->degree(x)==v, gens(R))
        i === nothing && error("morphism of grading groups can not be inferred from the ring map")
        push!(gg_imgs, degree(g[i]))
      end
      hom(G, H, gg_imgs)
    end, 
    coefficient_map=nothing
  )
  # Deflect to coercion if necessary
  all(x->parent(x) === S, g) || return hom(R, S, S.(g))

  # Since we have to use keyword arguments for the coefficient_map, we can not do anything better here:
  if coefficient_map === nothing
    return DecoratedRingMorphism(R, S, g, decoration_map)
  else 
    return DecoratedRingMorphism(R, S, g, decoration_map, coefficient_map)
  end
end

########################################################################
# Composition of morphisms of decorated rings
#
# In particular, we need to take care of the different constellations 
# with (non-)trivial coefficient maps, so there are four cases to cover
########################################################################
const _trivial_coeff_map = MPolyAnyMap{<:_DomainTypes, <:NCRing, Nothing, <:Any}
const _non_trivial_coeff_map = MPolyAnyMap{<:_DomainTypes, <:NCRing, <:Any, <:Any}
const _hecke_coeff_map = MPolyAnyMap{<:_DomainTypes, <:NCRing, <:Map, <:Any}

function compose(
    psi1::DecoratedRingMorphism{
              <:_DecoratedRingTypes, 
              <:_DecoratedRingTypes, 
              <:_trivial_coeff_map
            },
    psi2::DecoratedRingMorphism{
              <:_DecoratedRingTypes, 
              <:_DecoratedRingTypes, 
              <:_trivial_coeff_map
            }
  )
  g = psi2.(psi1.(gens(domain(psi1))))
  coeff_ring_map = coefficient_map(psi1)
  dec_map = compose(decoration_map(psi1), decoration_map(psi2))
  return DecoratedRingMorphism(domain(psi1), codomain(psi2), g, dec_map)
end

function compose(
    psi1::DecoratedRingMorphism{
              <:_DecoratedRingTypes, 
              <:_DecoratedRingTypes, 
              <:_trivial_coeff_map
            },
    psi2::DecoratedRingMorphism{
              <:_DecoratedRingTypes, 
              <:_DecoratedRingTypes, 
              <:_non_trivial_coeff_map
            }
  )
  g = psi2.(psi1.(gens(domain(psi1))))
  dec_map = compose(decoration_map(psi1), decoration_map(psi2))
  return DecoratedRingMorphism(domain(psi1), codomain(psi2), g, dec_map, x->coefficient_map(psi2)(x))
end

function compose(
    psi1::DecoratedRingMorphism{
              <:_DecoratedRingTypes, 
              <:_DecoratedRingTypes, 
              <:_non_trivial_coeff_map
            },
    psi2::DecoratedRingMorphism{
              <:_DecoratedRingTypes, 
              <:_DecoratedRingTypes, 
              <:_trivial_coeff_map
            }
  )
  g = psi2.(psi1.(gens(domain(psi1))))
  coeff_ring_map = coefficient_map(psi1)
  dec_map = compose(decoration_map(psi1), decoration_map(psi2))
  return DecoratedRingMorphism(domain(psi1), codomain(psi2), g, dec_map, coefficient_map(psi1))
end

function compose(
    psi1::DecoratedRingMorphism{
              <:_DecoratedRingTypes, 
              <:_DecoratedRingTypes, 
              <:_non_trivial_coeff_map
            },
    psi2::DecoratedRingMorphism{
              <:_DecoratedRingTypes, 
              <:_DecoratedRingTypes, 
              <:_non_trivial_coeff_map
            }
  )
  g = psi2.(psi1.(gens(domain(psi1))))
  coeff_ring_map = MapFromFunc(coefficient_ring(domain(psi1)), codomain(psi2), 
                               x->psi2(domain(psi2)(coefficient_map(psi1)(x))))
  dec_map = compose(decoration_map(psi1), decoration_map(psi2))
  return DecoratedRingMorphism(domain(psi1), codomain(psi2), g, dec_map, coeff_ring_map)
end


### Mapping of elements

function (psi::DecoratedRingMorphism)(a::_DecoratedRingElemTypes)
  parent(a) === domain(psi) || return psi(domain(psi)(a))

  return codomain(psi)(forget_grading(psi)(forget_grading(a)))
end

### Equality test
function ==(psi1::DecoratedRingMorphism, psi2::DecoratedRingMorphism)
  domain(psi1) === domain(psi2) || return false
  codomain(psi1) === codomain(psi2) || return false
  decoration_map(psi1) == decoration_map(psi2) || return false
  psi1.(gens(domain(psi1))) == psi2.(gens(domain(psi1))) || return false
  return true
end
