export SpecMor, identity_map, inclusion_map, inclusion_morphism, restrict, compose



########################################################################
# (1) General constructors
########################################################################

function SpecMor(
      X::AbsSpec,
      Y::AbsSpec{<:Ring, <:MPolyRing},
      f::Vector{<:RingElem};
      check::Bool=true
  )
  return SpecMor(X, Y, hom(OO(Y), OO(X), OO(X).(f)), check=check)
end

function SpecMor(
      X::AbsSpec,
      Y::AbsSpec,
      f::Vector{<:RingElem};
      check::Bool=true
  )
  return SpecMor(X, Y, hom(OO(Y), OO(X), OO(X).(f), check=check), check=check)
end

function SpecMor(
      X::AbsSpec,
      Y::AbsSpec,
      f::Vector;
      check::Bool=true
  )
  return SpecMor(X, Y, OO(X).(f), check=check)
end



########################################################################
# (2) Special constructors
########################################################################

# (2.1) Identity maps
identity_map(X::AbsSpec{<:Any, <:MPolyQuoLocalizedRing}) = SpecMor(X, X, hom(OO(X), OO(X), gens(ambient_ring(X)), check=false))
identity_map(X::AbsSpec{<:Any, <:MPolyLocalizedRing}) = SpecMor(X, X, hom(OO(X), OO(X), gens(ambient_ring(X)), check=false))
identity_map(X::AbsSpec{<:Any, <:MPolyRing}) = SpecMor(X, X, hom(OO(X), OO(X), gens(OO(X))))
identity_map(X::AbsSpec{<:Any, <:MPolyQuo}) = SpecMor(X, X, hom(OO(X), OO(X), gens(ambient_ring(X))))

# (2.2) Inclusion morphisms
inclusion_map(X::AbsSpec, Y::AbsSpec) = SpecMor(X, Y, gens(ambient_ring(Y)))  # TODO: Remove
inclusion_morphism(X::AbsSpec, Y::AbsSpec; check::Bool=true) = SpecMor(X, Y, gens(ambient_ring(Y)), check=check)

# (2.3) Composition
function compose(f::AbsSpecMor, g::AbsSpecMor; check::Bool=true)
  codomain(f) == domain(g) || error("Morphisms can not be composed")
  return SpecMor(domain(f), codomain(g), compose(pullback(g), pullback(f)), check=check)
end

# (2.4) Restriction of a morphism
function restrict(f::SpecMor, U::AbsSpec, V::AbsSpec; check::Bool=true)
  if check
    issubset(U, domain(f)) || error("second argument does not lay in the domain of the map")
    issubset(V, codomain(f)) || error("third argument does not lay in the codomain of the map")
    issubset(U, preimage(f, V)) || error("the image of the restriction is not contained in the restricted codomain")
  end
  return SpecMor(U, V, OO(U).(pullback(f).(gens(domain(pullback(f))))), check=check)
end
