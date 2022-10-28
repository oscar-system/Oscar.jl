export is_empty, issubset, is_open_embedding, is_closed_embedding



####################################################################################
# (1) Check if a scheme is empty
####################################################################################

Base.isempty(X::AbsSpec) = iszero(one(OO(X)))
is_empty(X::EmptyScheme) = true



####################################################################################
# (2) IsSubset for all other schemes
####################################################################################


# (2.0) For empty schemes and whenever issubset cannot be implemented
# (2.0) For empty schemes and whenever issubset cannot be implemented

# TODO: Does not comply with naming convention. Need is_subset etc.

@Markdown.doc """
    issubset(X::AbsSpec, Y::AbsSpec)

Checks whether ``X`` is a subset of ``Y`` based on the comparison of their coordinate rings.
"""
function issubset(X::AbsSpec, Y::AbsSpec)
  error("method `issubset(X, Y)` not implemented for `X` of type $(typeof(X)) and `Y` of type $(typeof(Y))")
end

issubset(X::EmptyScheme{BRT}, Y::Scheme{BRT}) where {BRT} = true

function issubset(Y::AbsSpec{BRT, <:Any}, X::EmptyScheme{BRT}) where {BRT}
  return iszero(one(OO(Y)))
end


# (2.1) MPolyRing in first argument
# (2.1) MPolyRing in first argument

function issubset(
    X::AbsSpec{BRT, RT},
    Y::AbsSpec{BRT, RT}
  ) where {BRT, RT<:MPolyRing}
  return OO(X) === OO(Y)
end


function issubset(
    X::AbsSpec{BRT, <:MPolyRing}, 
    Y::AbsSpec{BRT, <:MPolyQuo}
  ) where {BRT}
  R = OO(X)
  R === ambient_ring(Y) || return false
  return iszero(modulus(OO(Y)))
end


function issubset(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyLocalizedRing}
  ) where {BRT}
  R = OO(X)
  R === ambient_ring(Y) || return false
  return issubset(inverted_set(OO(Y)), units_of(R))
end


function issubset(
    X::AbsSpec{BRT, <:MPolyRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing}
  ) where {BRT}
  R = OO(X)
  R === ambient_ring(Y) || return false
  return issubset(inverted_set(OO(Y)), units_of(R)) && iszero(defining_ideal(Y))
end


# (2.2) MPolyQuo in first argument
# (2.2) MPolyQuo in first argument

function issubset(
    X::AbsSpec{BRT, <:MPolyQuo},
    Y::AbsSpec{BRT, <:MPolyRing}
  ) where {BRT}
  R = ambient_ring(Y)
  R === ambient_ring(X) || return false
  return true
end


function issubset(
    X::AbsSpec{BRT, RT},
    Y::AbsSpec{BRT, RT}
  ) where {BRT, RT<:MPolyQuo}
  R = ambient_ring(X)
  R === ambient_ring(Y) || return false
  return issubset(defining_ideal(Y), defining_ideal(X))
end


function issubset(
    X::AbsSpec{BRT, <:MPolyQuo},
    Y::AbsSpec{BRT, <:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  ) where {BRT}
  R = ambient_ring(X)
  R === ambient_ring(Y) || return false
  all(x->isunit(OO(X)(x)), denominators(inverted_set(OO(Y)))) || return false
  return iszero(localized_ring(OO(Y))(modulus(OO(X))))
end


function issubset(
    X::AbsSpec{BRT, <:MPolyQuo},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  ) where {BRT}
  R = ambient_ring(X)
  R === ambient_ring(Y) || return false
  all(x->isunit(OO(X)(x)), denominators(inverted_set(OO(Y)))) || return false
  return issubset(modulus(OO(Y)), localized_ring(OO(Y))(modulus(OO(X))))
end


# (2.3) MPolyLocalizedRing in first argument
# (2.3) MPolyLocalizedRing in first argument

function issubset(
    X::AbsSpec{BRT, <:MPolyLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyRing}
  ) where {BRT}
  R = OO(Y)
  R == ambient_ring(X) || return false
  return true
end


function issubset(
    X::AbsSpec{BRT, <:MPolyLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyQuo}
  ) where {BRT}
  R = ambient_ring(Y)
  R == ambient_ring(X) || return false
  return all(x->(iszero(OO(X)(x))), gens(modulus(OO(Y))))
end


function issubset(
    X::AbsSpec{BRT, RT},
    Y::AbsSpec{BRT, RT}
  ) where {BRT, RT<:MPolyLocalizedRing}
  R = ambient_ring(X)
  R === ambient_ring(Y) || return false
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  return issubset(UY, UX)
end


function issubset(
    X::AbsSpec{BRT, <:MPolyLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyQuoLocalizedRing}
  ) where {BRT}
  R = ambient_ring(X)
  R === ambient_ring(Y) || return false
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  return issubset(UY, UX) && iszero(modulus(OO(Y)))
end


# (2.4) MPolyQuoLocalizedRing in first argument
# (2.4) MPolyQuoLocalizedRing in first argument

function issubset(
    X::AbsSpec{BRT, <:MPolyQuoLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyRing}
  ) where {BRT}
  R = OO(Y)
  R == ambient_ring(X) || return false
  return true
end


function issubset(
    X::AbsSpec{BRT, <:MPolyQuoLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyQuo}
  ) where {BRT}
  R = ambient_ring(Y)
  R == ambient_ring(X) || return false
  L = localized_ring(OO(X))
  return issubset(L(modulus(OO(Y))), modulus(OO(X)))
end


function issubset(
    X::AbsSpec{BRT, <:MPolyQuoLocalizedRing},
    Y::AbsSpec{BRT, <:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  ) where {BRT}
  R = ambient_ring(X)
  R === ambient_ring(Y) || return false
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  if !issubset(UY, UX)
    # check whether the inverted elements in Y are units anyway
    for a in denominators(UY)
      is_unit(OO(X)(a)) || return false
    end
  end
  return true  # Spec R/I[S^-1] is a closed subscheme of Spec R[S^-1]
end


function issubset(
    X::AbsSpec{BRT, RT},
    Y::AbsSpec{BRT, RT}
  ) where {BRT, RT<:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, <:MPolyPowersOfElement}}
  R = ambient_ring(X)
  R === ambient_ring(Y) || return false
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  if !issubset(UY, UX)
    # check whether the inverted elements in Y are units anyway
    for a in denominators(UY)
      is_unit(OO(X)(a)) || return false
    end
  end
  J = localized_ring(OO(X))(modulus(quotient_ring(OO(Y))))
  return issubset(J, modulus(OO(X)))
end



####################################################################################
# (3) Check if a scheme is openly embeeded in another scheme
####################################################################################

#TODO: Add more cross-type methods as needed.

@Markdown.doc """
    is_open_embedding(X::AbsSpec, Y::AbsSpec)

Checks whether ``X`` is openly embedded in ``Y``.
"""
function is_open_embedding(X::AbsSpec, Y::AbsSpec)
  return is_open_embedding(standard_spec(X), standard_spec(Y))
end


function is_open_embedding(
    X::Spec{BRT, RT},
    Y::Spec{BRT, RT}
  ) where {BRT, RT<:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any,
                                          <:MPolyPowersOfElement}}
  R = ambient_ring(X)
  R === ambient_ring(Y) || return false
  UX = inverted_set(OO(X))
  UY = inverted_set(OO(Y))
  issubset(UY, UX) || return false
  J = localized_ring(OO(X))(modulus(quotient_ring(OO(Y))))
  return modulus(OO(X)) == J 
end


function is_open_embedding(
    X::Spec{BRT, <:MPolyQuoLocalizedRing},
    Y::Spec{BRT, <:MPolyRing}
  ) where {BRT}
  return OO(Y) == ambient_ring(X) && all(iszero, gens(modulus(OO(X))))
end



####################################################################################
# (4) Check if a scheme can be embedded via a closed embeeded in another scheme
####################################################################################

@Markdown.doc """
    is_closed_embedding(X::AbsSpec, Y::AbsSpec)

Checks whether ``X`` is closed embedded in ``Y``.
"""
function is_closed_embedding(X::AbsSpec, Y::AbsSpec)
  error("`is_closed_embedding(X, Y)` not implemented for X of type $(typeof(X)) and Y of type $(typeof(Y))")
end


function is_closed_embedding(
    X::AbsSpec{<:Ring, <:MPolyQuo},
    Y::AbsSpec{<:Ring, <:MPolyRing}
  )
  return ambient_ring(X) === ambient_ring(Y)
end


function is_closed_embedding(
    X::AbsSpec{<:Ring, <:MPolyRing},
    Y::AbsSpec{<:Ring, <:MPolyQuo}
  )
  R = ambient_ring(X)
  R === ambient_ring(Y) || return false
  return iszero(modulus(OO(Y)))
end


function is_closed_embedding(
    X::AbsSpec{<:Ring, <:MPolyQuo},
    Y::AbsSpec{<:Ring, <:MPolyQuo}
  )
  R = ambient_ring(X)
  R === ambient_ring(Y) || return false
  return issubset(modulus(OO(Y)), modulus(OO(X)))
end


function is_closed_embedding(
    X::AbsSpec{<:Ring, <:MPolyRing},
    Y::AbsSpec{<:Ring, <:MPolyRing}
  )
  R = ambient_ring(X)
  R === ambient_ring(Y) || return false
  return true
end


function is_closed_embedding(
    X::AbsSpec{<:Ring, <:MPolyLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                            <:MPolyPowersOfElement}},

    Y::AbsSpec{<:Ring, <:MPolyRing}
  )
  R = ambient_ring(X)
  R === ambient_ring(Y) || return false
  for f in inverted_set(OO(X))
    isunit(OO(Y)(f)) || return false
  end
  return true
end


function is_closed_embedding(
    X::AbsSpec{<:Ring, <:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any, 
                                              <:MPolyPowersOfElement}},

    Y::AbsSpec{<:Ring, <:MPolyQuo}
  )
  R = ambient_ring(X)
  R === ambient_ring(Y) || return false
  for x in inverted_set(OO(X)) 
    isunit(OO(Y)(x)) || return false
  end
  for g in gens(modulus(OO(Y)))
    iszero(OO(X)(g)) || return false
  end
  return true
end


function is_closed_embedding(
    X::AbsSpec{BRT, RT},
    Y::AbsSpec{BRT, RT}
  ) where {BRT, RT<:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any,
                                          <:MPolyPowersOfElement}}
  R = ambient_ring(X)
  R === ambient_ring(Y) || return false
  inverted_set(OO(X)) == inverted_set(OO(Y)) || return false
  J = localized_ring(OO(X))(modulus(quotient_ring(OO(Y))))
  return issubset(J, modulus(OO(X)))
end


function is_closed_embedding(
    X::Spec{BRT, <:MPolyQuo},
    Y::Spec{BRT, <:MPolyRing}
  ) where {BRT}
  OO(Y) === ambient_ring(X) || return false
  return true
end


function is_closed_embedding(
    X::Spec{BRT, <:MPolyQuo},
    Y::Spec{BRT, <:RT}
  ) where {BRT, RT<:MPolyQuoLocalizedRing{<:Any, <:Any, <:Any, <:Any,
                                          <:MPolyPowersOfElement}}
  R = ambient_ring(X)
  R === ambient_ring(Y) || return false
  all(x->(isunit(OO(X)(x))), denominators(inverted_set(OO(Y)))) || return false
  return issubset(modulus(OO(Y)), localized_ring(OO(Y))(modulus(OO(X))))
end
