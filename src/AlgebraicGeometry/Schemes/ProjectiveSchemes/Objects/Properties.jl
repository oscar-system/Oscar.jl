#######################################################################
# Basic properties
#######################################################################
@attr Bool function is_empty(P::AbsProjectiveScheme{<:Field})
  I = defining_ideal(P)
  R = base_ring(I)
  J = ideal(R, gens(R))
  return is_one(saturation(I, J))
end

@attr Int function dim(P::AbsProjectiveScheme{<:Field})
  return dim(defining_ideal(P))-1
end

@attr QQPolyRingElem function hilbert_polynomial(P::AbsProjectiveScheme{<:Field})
  return hilbert_polynomial(homogeneous_coordinate_ring(P))
end

@attr ZZRingElem function degree(P::AbsProjectiveScheme{<:Field})
  return degree(homogeneous_coordinate_ring(P))
end

@attr QQFieldElem function arithmetic_genus(P::AbsProjectiveScheme{<:Field})
  h = hilbert_polynomial(P)
  return (-1)^dim(P) * (first(coefficients(h)) - 1)
end

# TODO: Jacobi criterion
function is_smooth(P::AbsProjectiveScheme; algorithm=:covered)
  get_attribute!(P, :is_smooth) do
    if algorithm == :covered
      return is_smooth(covered_scheme(P))
    end
    if algorithm == :jacobi
      error("jacobi criterion not implemented")
    end
  end
end

@attr Bool function is_irreducible(P::AbsProjectiveScheme)
  is_empty(P) && return false  # the empty set is not irreducible
  I = defining_ideal(P)
  return is_prime(radical(I)) # topological property
end

@attr Bool function is_reduced(P::AbsProjectiveScheme)
  I = defining_ideal(P)
  return is_radical(I)
end

@attr Bool function is_geometrically_reduced(P::AbsProjectiveScheme{<:Field})
  if is_perfect(base_ring(P))
    return is_reduced(P)
  end
  throw(NotImplementedError(:is_geometrically_reduced,"currently we can decide this only over a perfect base field"))
end

@attr Bool function is_geometrically_irreducible(P::AbsProjectiveScheme{<:Field, <:MPolyQuoRing{<:MPolyDecRingElem}})
  is_empty(P) && return false # the empty set is not irreducible
  I = defining_ideal(P)
  is_prime(I) || return false
  AI = absolute_primary_decomposition(I)
  return length(AI)==1
end

@attr Bool function is_geometrically_irreducible(X::AbsProjectiveScheme{<:Field, <:MPolyDecRing})
  return true # projective space is geometrically irreducible
end

@attr Bool function is_integral(X::AbsProjectiveScheme{<:Field, <:MPolyAnyRing})
  !is_empty(X) || return false
  if has_attribute(X,:is_reduced) && has_attribute(X,:is_irreducible)
     return get_attribute(X,:is_reduced) && get_attribute(X,:is_irreducible)
  end
  I = defining_ideal(X)
  return is_prime(I)
end

@attr Bool function is_geometrically_integral(X::AbsProjectiveScheme{<:Field, <:MPolyDecRing})
  return true # projective space is geometrically integral
end

@attr Bool function is_geometrically_integral(X::AbsProjectiveScheme{<:Field})
  is_integral(X) || return false
  throw(NotImplementedError(:is_geometrically_integral, "absolute primary decomposition is currently only available over the rationals"))
end

@attr Bool function is_geometrically_integral(X::AbsProjectiveScheme{<:QQField, <:MPolyQuoRing{<:MPolyDecRingElem}})
  is_integral(X) || return false
  I = defining_ideal(X)
  AI = absolute_primary_decomposition(I)
  @assert length(AI)==1 # it is prime since X is integral
  return AI[1][4]==1
end

